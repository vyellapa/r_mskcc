##' Something happened and manipulating DataFrame objects in the "normal"
##' annotate reads business has been chewing up tons of memory.
tallyTagDistro <- function(expt, annotated.genome, unique.only=TRUE,
                           max.mismatch=NULL,
                           chrs=paste('chr', c(1:22, 'X'), sep="")) {
  tally <- mclapply(chrs, function(chr) {
    cat("===", chr, "===\n")
    x <- getReadsOnChromosome(expt, chr, unique.only=unique.only,
                              max.mismatch=NULL)

    values(x) <- NULL
    x <- resize(x, width=1, fix='center')
    chr.anno <- annotated.genome[seqnames(annotated.genome) == chr]
    anno.sense <- matchToAnnotation(x, chr.anno)
    anno.sense <- anno.sense$exon.anno

    anti.consider <- anno.sense %in% c('intergenic', 'utr5*', 'utr3*')
    tmp <- swapStrand(x[anti.consider])
    anti.anno <- matchToAnnotation(tmp, annotated.genome)$exon.anno
    is.anti <- !(anti.anno %in% c('intergenic', 'utr3*', 'utr5*'))
    anno.sense[which(anti.consider)[is.anti]] <- 'antisense'
    table(anno.sense)
  })

  cats <- c(antisense=0L, cds=0L, intergenic=0L, intron=0L,
            overlap=0L, utr=0L, utr3=0L, "utr3*"=0L, utr5=0L,
            "utr5*"=0L)

  for (xx in tally) {
    xref <- match(names(xx), names(cats))
    cats[xref] <- cats[xref] + xx
  }

  cats
}

##' Read processing report
##'
##' The nut of this function is moved to the command line utility
##'   TagSeq/inst/scripts/tagseq-qc-libstats.R
##'
##' Summarizes number of reads lost/retained at each step of the processing
##' pipeline.
##'
##' (1) no. input reads
##' (2) no. reads after quality trimming
##' (3) no. reads after adapter trimming
##' (4) no. reads after AAA trimming
##' (5) no. reads that align
##' (6) no. of unique reads aligned
##'
##' @param preproc.dir The directory that the preprocessing script was run
##' on the raw reads (ie multiplex-preprocess)
##' @param align.dir The directory the results from bwa2bam are stored in,
##' ie. where the metadata.yaml file is
readProcessingReport <- function(preproc.dir, align.dir,
                                 preproc.log='processing-log.txt',
                                 align.log='metadata.yaml') {
  preproc <- readLines(file.path(preproc.dir, preproc.log))
  meta <- readLines(file.path(align.dir, align.log))

  ## Number of reads from (1) input and (2) adapter trimming
  pre.in <- sapply(preproc[grep('input', preproc, ignore.case=TRUE)],
                   function(line) {
                     as.integer(strsplit(line, ' ')[[1]][[2]])
                   })
  names(pre.in) <- c('quality', 'adapter')
  homo.A.in <- preproc[grep('number of reads processed', preproc,
                            ignore.case=TRUE)]
  homo.A.in <- as.integer(tail(strsplit(homo.A.in, ' ')[[1]], 1))
  homo.A.trash <- preproc[grep('trashed', preproc, ignore.case=TRUE)]
  homo.A.trash <- as.integer(tail(strsplit(homo.A.trash, ' ')[[1]], 1))
  reads.2.align <- homo.A.in - homo.A.trash

  n.aligned <- meta[grep('read.count:', meta, ignore.case=TRUE)]
  n.aligned <- as.integer(tail(strsplit(n.aligned, ' ')[[1]], 1))

  n.unique <- meta[grep('count.unique:', meta, ignore.case=TRUE)]
  n.unique <- as.integer(tail(strsplit(n.unique, ' ')[[1]], 1))

  df <- data.frame(step=c('start', 'after.quality', 'after.adapter', 'after.A',
                     'aligned', 'unique'),
                   count=c(
                     pre.in['quality'],
                     pre.in['adapter'],
                     homo.A.in,
                     reads.2.align,
                     n.aligned,
                     n.unique))
  df
}

