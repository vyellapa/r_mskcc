library(GenomicRanges)
library('BSgenome.Hsapiens.UCSC.hg19')

args = commandArgs(TRUE)
bed.file = toString(args[1])
outfile = toString(args[2])

regions = read.table(file=bed.file, header=F)
regions.gr = with(regions, GRanges(V1, IRanges(V2, V3)))
seq = getSeq(Hsapiens, paste('chr',as.character(seqnames(regions.gr)),sep=""), start=start(regions.gr), end(regions.gr))
regions = cbind(regions,as.data.frame(seq))

write.table(file=outfile, regions, row.names=F, col.names=F, sep="\t", quote=F)
