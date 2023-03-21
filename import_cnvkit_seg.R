# Dependencies
suppressMessages(library(dplyr))
suppressMessages(library(stringi))
suppressMessages(library(GenomicRanges))
options(scipen=999)

import_cnvkit_seg <- function(cytoband, path){
    # Import and process segmentation file from CNVkit
    cnvkit <- read.csv(path, sep = "\t", stringsAsFactors = F)
    cnvkit <- cnvkit[-5]
    names(cnvkit) <- c('LEUKID', 'chrom', 'start', 'end', 'seg.mean')
    
    # Annotate cnvkit data with cytoband
    arms_gr <- with(arms, GRanges(chrom, IRanges(start = start, end = end, names = arm)))
    cnvkit_gr <- with(cnvkit, GRanges(chrom, IRanges(start = start, end = end, names = seg.mean)))
    cnvkit_ol <- findOverlaps(query = cnvkit_gr, subject = arms_gr, type = "within")
    cnvkit <- data.frame(cnvkit[queryHits(cnvkit_ol),], arms[subjectHits(cnvkit_ol),])
    cnvkit <- cnvkit[-c(grep('.1', names(cnvkit)))]
    cnvkit <- mutate(cnvkit, cnv_size = end - start)
    cnvkit <- mutate(cnvkit, cnv_prop = cnv_size/arm_size)

    # Simplified CN state (gain or loss)
    cnvkit <- mutate(cnvkit, cn_simple = factor(ifelse(seg.mean > 0, 1, 2),
                                               labels = c('Gain', 'Del'),
                                               levels = c(1, 2)))
    cnvkit <- mutate(cnvkit, CNV = paste(cn_simple, '(', chrom, arm, ')', sep=''))
    return(cnvkit)
}
