library(ggplot2)
library(sciClone)
setwd("~/Desktop/scratch/")
for(i in grep("new",list.files(pattern = "T2.*_allDirichletProcessInfo.txt", recursive = T), value=T)) {
  T1=gsub("-T2-","-T1-",i)
  #j=list.files(pattern=T2,recursive = T)
  print(i)
  print(T1)
  
  name=strsplit(i,"/")[[1]][3]
  t2c=gsub("_allDirichletProcessInfo.txt",".*.caveman.tsv.gz",strsplit(i,"/")[[1]][4])
  print(t2c)
  t2c=list.files(pattern=t2c)
  t1c=gsub("-T2-", "-T1-",t2c)
  print(t2c)
  
  caveman_t1c=read.table(gzfile(t1c), sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(caveman_t1c)=paste("T1_CAVE",colnames(caveman_t1c),sep="_")
  caveman_t2c=read.table(gzfile(t2c), sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(caveman_t2c)=paste("T2_CAVE",colnames(caveman_t2c),sep="_")
  caveman_t1c$merger=paste(caveman_t1c$T1_CAVE_CHR, caveman_t1c$T1_CAVE_START,sep=":")
  caveman_t2c$merger=paste(caveman_t2c$T2_CAVE_CHR, caveman_t2c$T2_CAVE_START,sep=":")
  
  print(t2c)
  
  ii=read.table(sep="\t", header=T,T1)
  colnames(ii)=paste("T1",colnames(ii),sep="_")
  jj=read.table(sep="\t", header=T,i)
  colnames(jj)=paste("T2",colnames(jj),sep="_")
  ii$merger=paste(ii$T1_chr, ii$T1_end,sep=":")
  jj$merger=paste(jj$T2_chr, jj$T2_end,sep=":")
  merged=merge(ii, jj, by.x="merger", by.y="merger",all=T)
  merged=merge(merged, caveman_t1c,by.x="merger", by.y="merger", all.x = T )
  merged=merge(merged, caveman_t2c,by.x="merger", by.y="merger", all.x = T )
  
  merged[is.na(merged$T1_CAVE_TARGET_VAF),]$T1_CAVE_TARGET_VAF=0.00
  merged[is.na(merged$T2_CAVE_TARGET_VAF),]$T2_CAVE_TARGET_VAF=0.00
  
  v1=merged[,c("T1_chr", "T1_end", "T1_WT.count", "T1_mut.count", "T1_subclonal.fraction")]
  v2=merged[,c("T2_chr", "T2_end", "T2_WT.count", "T2_mut.count", "T2_subclonal.fraction")]
  colnames(v1)=c("chr","pos", "ref_reads", "var_reads", "vaf")
  colnames(v2)=c("chr","pos", "ref_reads", "var_reads", "vaf")
  v1$vaf=(v1$vaf/2)*100
  v2$vaf=(v2$vaf/2)*100
  
  
  cn1=merged[,c("T1_chr", "T1_start", "T1_end")]
  cn1$segMean=rep("2",nrow(cn1))
  colnames(cn1)=c("chr", "start", "stop", "segment_mean")
  cn1$start=cn1$stop
  cn1$stop=cn1$stop
  
  sc=sciClone(vafs=list(v1,v2), copyNumberCalls = list(cn1,cn1), sampleNames = c("T1","T2"))
  writeClusterTable(sc, paste(name,"clusters","txt",sep="."))
  sc.plot2d(sc, paste(name,"pdf",sep="."))
  
}

library(reshape)
library(ggplot2)
require(scales)
t=read.table("plot.txt", sep=" ", header=T)
m=melt(t)

m$value=m$value/1000000

ggplot(m, aes(NAME, value,fill=variable))+
  geom_bar(position="dodge",stat="identity",alpha=0.85)+theme_bw()+theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="top")+
  scale_fill_manual(values = c("#64572D", "#818854", "#DA8D3B","#E25800"))+
  ylab("Reads in Million")+xlab("")
  

