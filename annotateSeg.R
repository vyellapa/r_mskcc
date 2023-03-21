library(GenomicRanges)

## Usage 
#Rscript --vanilla gtf2bed.R <input.seg> <seg_annotated.out>
args=commandArgs(TRUE)
s1=args[1] 
s2=args[2]
unlink(s2)

x1=read.table("~/local/downloads/Homo_sapiens.GRCh37.75.min_max.bed", sep="\t", header=F)
#x1=read.table(s1, sep="\t", header=F)

x2=read.table(s1, sep="\t", header=T)
x2[x2$chrom==23,]$chrom="X"


gr1=with(x1, GRanges(V1, IRanges(start=V2, end=V3, name=V6)))
cat("gr1 done")

gr2=with(x2, GRanges(chrom, IRanges(start=start, end=end, names=cnlr.median.clust)))


overLap=findOverlaps(query=gr1, subject=gr2)
overLap.df = data.frame(x1[queryHits(overLap),], x2[subjectHits(overLap),])


write.table(overLap.df,file=s2, append=TRUE, quote=F,sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
