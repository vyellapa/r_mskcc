## Identify and annotate peaks peaks from command line script
##   $ Rscript /Volumes/IrtishaDrive/APA/tagseq/inst/scripts/peaks-identify -v -t 1 /Volumes/IrtishaDrive/APA/sample_run/identify-peaks.yaml
##
## Quantify peaks
##   $ Rscript /Volumes/IrtishaDrive/APA/tagseq/inst/scripts/peaks-quantify -v -t 1 /Volumes/IrtishaDrive/APA/sample_run/quantify-peaks.yaml
##
## The last step will produce a "quantified-peaks.rds" file that has the
## all peaks, their genomic annotation, and quantification across experiments
##

library(TagSeq)
library(rtracklayer)
library(data.table)
library(SummarizedExperiment)
source(system.file("scripts/required_functions.R", package="TagSeq"))

genome <- 'hg19'
annotation.file <- system.file("extdata",
                          sprintf("annotated.genome.%s.rds", genome),
                          package="TagSeq")
annotation <- readRDS(annotation.file)
bsg <- getBsGenome(genome)
a.window.count <- 5L
ip.window <- 6L
ip.distal.downstream <- 12L
ip.internal.downstream <- 12L
min.upstream <- 10L
parallel <- TRUE
threads <- 10
polya.signals <- kPolyA$dna.signal
rescue.from.end <- 40
library(doParallel)
registerDoParallel(cores=threads)

peaks <- readRDS("quantified-peaks.rds")
ip.summary <- internalPrimingSummary(rowRanges(peaks), bsg, quantile.start=NULL,
                                     quantile.end=NULL, quantile.offset=NULL,
                                     .parallel=parallel, internal.flag.by='score',
                                     a.window.count=a.window.count,
                                     ip.window=ip.window,
                                     ip.distal.downstream=ip.distal.downstream,
                                     ip.internal.downstream=ip.internal.downstream,
                                     min.upstream=min.upstream)
saveRDS(ip.summary, 'edge.internal.priming.rds')

## Get some polyA position stats
pa.stats <- collectPolyAPositionStatistics(rowRanges(peaks), bsg,
                                           polya=polya.signals,
                                           .parallel=parallel)
saveRDS(pa.stats, 'polya.signals.rds')

ip.summary <- readRDS('edge.internal.priming.rds')
pa.stats <- readRDS('polya.signals.rds')
ip.info <- combinePrimingInfo(rowRanges(peaks), ip.summary, pa.stats,
                              anno.rescue=annotation,
                              rescue.from.end=rescue.from.end)
if (is(ip.info, 'DataFrame') && nrow(ip.info) == nrow(peaks)) {
  values(peaks) <- cbind(values(peaks), ip.info)
  saveRDS(peaks, 'annotated-peaks-with-ip.rds')
}

################################################################################
## Creating an atlas of robust APA events in the UTR3 and introns
## 1) Remove the peaks in the blacklisted region and ip.axe and change the
## the library size accordingly as these are resulting from artifacts
## 2) Remove peaks that are not associated with any gene
## 3) Group apa events based on the regions where they occur
se.peaks1 <- readRDS('annotated-peaks-with-ip.rds')
se.peaks2 <- removePeaks(se.peaks1, blacklisted = TRUE, ip.axe = TRUE, 
  remove.exon.anno = c("antisense"), change.library.size = TRUE, genome)
se.peaks3 <- se.peaks2[!is.na(values(se.peaks2)$entrez.id)]
se.peaks4 <- assignCleavageUnit(se.peaks3)

################################################################################
## 4) Calculate the TPM
## 5) Calculate usage of UTR3 and intronic events
## 6) Filter for robust UTR3 events
## 7) Filter for robust intronic events
## 8) Filter for genic events
m.tpm <- getTPM(se.peaks4)
assay(se.peaks4, 'tpm') <- m.tpm[, colnames(se.peaks4)]
se.peaks5 <- calculateUTR3AndIntronUsage(se.peaks4, usage.on = "cleavage.unit")


se.peaks6 <- GroupFilterOnUsageAndTPM(se.peaks5, enames = colnames(se.peaks5),
  unit = "utr3", tpm.cutoff = 3, ui.cutoff = 0.1)

##Now filter for intronic
se.peaks7 <- GroupFilterOnUsageAndTPM(se.peaks6, enames = colnames(se.peaks6),
  unit = "intronic", tpm.cutoff = 5, ui.cutoff = 0.1)

se.peaks8 <- GroupFilterOnUsageAndTPM(se.peaks7, enames = colnames(se.peaks7),
  unit = c("genic"), tpm.cutoff = 5, ui.cutoff = 0)

################################################################################
## 9) Cluster the closeby peaks and reannotate them
## 10) Repeat steps 2-5 as the usage of ApA events changes based on the total numbers events in the atlas
se.peaks9 <- clusterPeaks(se.peaks8, max.distance = 200L)
se.peaks10 <- annotateSummarizedExperiment(se.peaks9, annotations = annotation)
se.peaks11 <- se.peaks10[!is.na(values(se.peaks10)$entrez.id)]
se.peaks12 <- assignCleavageUnit(se.peaks11)
se.peaks13 <- calculateUTR3AndIntronUsage(se.peaks12, usage.on = "cleavage.unit")

################################################################################
## 11) Identify the proximal and distal peak in UTR3
## 12) Calculate the cleavage index
se.peaks14 <- prepSE4ApaTesting(se.peaks13)
se.peaks15 <- calculateCI(se.peaks14, 'proximal')
se.peaks16 <- calculateCI(se.peaks15, 'distal')
saveRDS(se.peaks16, 'SE.prep-apa.rds')




