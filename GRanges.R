library("BSgenome.Hsapiens.UCSC.hg19")
setwd("~/Desktop/SOHN_163/brass")

#args <- commandArgs(trailingOnly = TRUE)
#T1 <- args[1] #T1
#T2 <- args[2] #T2
#T3 <- args[3] #T3


getFasta <- function(svo) {
  cat(svo)
  svo
  s=getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, pos-15, pos+15 )
  len=length(s)
  header=paste("> ",chr,":",pos," length=",len,sep="")
  cat(header)
  cat(s)
}


T1="E-H-109099-T1-1-D1-1_vs_E-H-109099-N1-1-D1-1.annot.bedpe"
T2="E-H-109099-T2-1-D1-1_vs_E-H-109099-N1-1-D1-1.annot.bedpe"
  b=basename(T2)
  b=strsplit(b,"_")[[1]][1]
  p=paste("cat",T1,"|sed \"s/# chr/chr/g\" > temp",sep=" ")
  print(p)
  system(p)
  t1_sv=read.table("temp", sep="\t", header=T, stringsAsFactors = F,quote="~")
  
  p=paste("cat",T2,"|sed \"s/# chr/chr/g\" > temp",sep=" ")
  print(p)
  system(p)
  t2_sv=read.table("temp", sep="\t", header=T, stringsAsFactors = F,quote="~")
  
  t1_sv$merger=paste(t1_sv$chr1,t1_sv$start1,t1_sv$end1,t1_sv$chr2,t1_sv$start2,t1_sv$end2,sep=":")
  t2_sv$merger=paste(t2_sv$chr1,t2_sv$start1,t2_sv$end1,t2_sv$chr2,t2_sv$start2,t2_sv$end2,sep=":")
  
  g1=GRanges(seqnames= Rle( as.character(t1_sv$chr1)), ranges = IRanges(start = t1_sv$start1, end=t1_sv$end1), merger=t1_sv$merger)
  g2=GRanges(seqnames= Rle( as.character(t2_sv$chr1)), ranges = IRanges(start = t2_sv$start1, end=t2_sv$end1), merger=t2_sv$merger)
  
  
  g12=GRanges(seqnames= Rle( as.character(t1_sv$chr2)), ranges = IRanges(start = t1_sv$start2, end=t1_sv$end2), merger=t1_sv$merger)
  g22=GRanges(seqnames= Rle( as.character(t2_sv$chr2)), ranges = IRanges(start = t2_sv$start2, end=t2_sv$end2), merger=t2_sv$merger)
  
  g1bed=GRanges(seqnames= Rle( as.character(t1_sv$chr1)), ranges = IRanges(start = t1_sv$start1-100, end=t1_sv$end1100))
  g2bed=GRanges(seqnames= Rle( as.character(t2_sv$chr1)), ranges = IRanges(start = t2_sv$start1-100, end=t2_sv$end1+100))
  g22bed=GRanges(seqnames= Rle( as.character(t2_sv$chr2)), ranges = IRanges(start = t2_sv$start2-100, end=t2_sv$end2+100))
  
  hits=findOverlaps(g1,g2bed)
  new=as.data.frame(t1_sv[queryHits(hits),])
  new=new[!duplicated(new),]
  gnew=GRanges(seqnames= Rle( as.character(new$chr2)), ranges = IRanges(start = new$start2, end=new$end2))
  
  hits=findOverlaps(gnew,g22bed)
  new=as.data.frame(new[queryHits(hits),])
  new=new[!duplicated(new),]
  
  t1_sv$merger=paste(t1_sv$chr1,t1_sv$start1,t1_sv$end1,t1_sv$chr2,t1_sv$start2,t1_sv$end2,sep=":")
  t2_sv$merger=paste(t2_sv$chr1,t2_sv$start1,t2_sv$end1,t2_sv$chr2,t2_sv$start2,t2_sv$end2,sep=":")
  
  common=merged[!is.na(merged$chr1.x),]
  common=common[!is.na(common$chr1.y),]
  common=common[,-1]
  
  uniqT1=merged[is.na(merged$chr1.y),-1]
  uniqT2=merged[is.na(merged$chr1.x),-1]
  uniqT2=uniqT2[,c(47:92)]
  #sv$chr2=paste("chr",sv$chr2,sep="")
  colnames(uniqT1)=colnames(t1_sv)
  colnames(uniqT2)=colnames(t1_sv)
  
  write.table(common,file="E99_common.txt", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
  write.table(uniqT1,file="E99_uniqT1.txt", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)
  write.table(uniqT2,file="E99_uniqT2.txt", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)


  
  
  #library(RCircos);
  #library(scales)
  
  #library(BSgenome.Hsapiens.UCSC.hg19) # human genome
  #library(VariantAnnotation) # 
  
  
