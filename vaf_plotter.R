library(GenomicRanges)
library(stringr)
library(ggplot2)
args=commandArgs(TRUE)
vfile=args[1]
cavefile=args[2]
cnfile=args[3]


name=strsplit(vfile,'.caveman')[[1]][1]
vcf=read.table(gzfile(vfile), sep="\t", header=F)

caveman=read.table(gzfile("E-H-109096-T1-1-D1-1_vs_E-H-109096-N1-1-D1-1.caveman.tsv.gz"), sep="\t", header=T)
cn=read.table(gzfile("E-H-109096-T1-1-D1-1_subclones.txt"), sep="\t", header=T)
cn=cn[which(cn$LogR<0.1 & cn$LogR>-0.1), ]

cnRange=GRanges(seqnames= Rle( as.character(cn$chr)), ranges = IRanges(start = cn$startpos, end=cn$endpos))


vcf=vcf[vcf$V7=="PASS",]

ch=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
pattern="ASRD=[0-9]+.[0-9]+"
vcf$ASRD=gsub("ASRD=","",str_extract(vcf$V8, pattern))
vcf=vcf[vcf$ASRD>0.93,]
vcf=vcf[vcf$V1 %in% ch, ]
vcfRange=GRanges(seqnames= Rle( as.character(vcf$V1)), ranges = IRanges(start = (vcf$V2-1), end=vcf$V2), ID=as.character(vcf$V3))

hits=findOverlaps(vcfRange, cnRange)
over=vcf[queryHits(hits), c(1,2,3,4,5)]

plot=caveman[caveman$ID_VARIANT %in% over$V3,]$TARGET_VAF

plot=as.data.frame(caveman[caveman$ID_VARIANT %in% over$V3,]$TARGET_VAF)
ggplot(plot, aes(plot$VAF))+geom_density(fill="blue", alpha=0.6)+ggtitle(name)+annotate("text", label=paste("N =",nrow(plot)),x=0.1,y=3)
