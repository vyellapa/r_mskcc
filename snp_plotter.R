library(reshape2)
library(ggplot2)
library(VennDiagram)
library("stringr")
setwd("/Users/yellapav/Desktop/SOHN_163/one_offs/SNP")
#args=commandArgs(TRUE)
#i=args[1]
#j=args[2]

i="E-H-109100-T2-1-D1-1_vs_E-H-109100-N1-1-D1-1.snps.ids.vcf.gz"
name=strsplit(i,"_",fixed=TRUE)[[1]][1]
snp=read.table(i,sep="\t",header=F)
snv=read.table("E-H-109100-T2-1-D1-1_vs_E-H-109100-N1-1-D1-1.caveman.muts.annot.vcf.gz",sep="\t",header=F)
snv=read.table(j,sep="\t",header=F)

pat<-'MP=[0-9,e,\\-,\\+,\\.]+'
snp$MP=gsub("MP=","",str_extract(snp$V8, pat))

foo <- data.frame(do.call('rbind', strsplit(as.character(snp$V10),':',fixed=TRUE)))
colnames(foo)=paste("NORMAL_",colnames(foo),sep="")
snp=cbind(snp,foo)
foo <- data.frame(do.call('rbind', strsplit(as.character(snp$V11),':',fixed=TRUE)))
colnames(foo)=paste("TUMOR_",colnames(foo),sep="")
snp=cbind(snp,foo)
snp$CALL="SNP"
#snp=snp[snp$NORMAL_X10<0.01]

pat<-'MP=[0-9,e,\\-,\\+,\\.]+'
snv$MP=gsub("MP=","",str_extract(snv$V8, pat))

foo <- data.frame(do.call('rbind', strsplit(as.character(snv$V10),':',fixed=TRUE)))
colnames(foo)=paste("NORMAL_",colnames(foo),sep="")
snv=cbind(snv,foo)
foo <- data.frame(do.call('rbind', strsplit(as.character(snv$V11),':',fixed=TRUE)))
colnames(foo)=paste("TUMOR_",colnames(foo),sep="")
snv=cbind(snv,foo)
snv$CALL="SNV"

snp=rbind(snp,snv)

snp$NORMAL_DP=as.integer(as.character(snp$NORMAL_X2))+as.integer(as.character(snp$NORMAL_X3))+as.integer(as.character(snp$NORMAL_X4))+as.integer(as.character(snp$NORMAL_X5))+as.integer(as.character(snp$NORMAL_X6))+as.integer(as.character(snp$NORMAL_X7))+as.integer(as.character(snp$NORMAL_X8))+as.integer(as.character(snp$NORMAL_X9))
snp$TUMOR_DP=as.integer(as.character(snp$TUMOR_X2))+as.integer(as.character(snp$TUMOR_X3))+as.integer(as.character(snp$TUMOR_X4))+as.integer(as.character(snp$TUMOR_X5))+as.integer(as.character(snp$TUMOR_X6))+as.integer(as.character(snp$TUMOR_X7))+as.integer(as.character(snp$TUMOR_X8))+as.integer(as.character(snp$TUMOR_X9))
snp$TUMOR_X10=as.numeric(as.character(snp$TUMOR_X10))
snp$NORMAL_X10=as.numeric(as.character(snp$NORMAL_X10))
snp$MP=as.numeric(as.character(snp$MP))




plot=snp[which(snp$CALL=="SNP"),]
plot=plot[which(plot$TUMOR_DP>9 & plot$NORMAL_DP>9),]
#plot=plot[which(plot$NORMAL_X1=="0|0" & (plot$NORMAL_X10<=0.001 | plot$NORMAL_X10 > 0.999)),]
plot=plot[which(plot$NORMAL_X10<=0.00),]
#plot=plot[which((plot$TUMOR_X10>0.05 & plot$TUMOR_X10<0.95)|((plot$TUMOR_X1=="1|1" | plot$TUMOR_X1=="0|1") & plot$TUMOR_X10>0.9)),]
plot=plot[which((plot$TUMOR_X10>0.05)),]

p=ggplot(plot,aes(TUMOR_X10,fill=TUMOR_X1))+geom_density(alpha=0.65,linetype="blank")+theme_bw()+
  geom_vline(xintercept = 0.05,color="blue")+geom_vline(xintercept = 0.95,color="blue")+scale_fill_brewer(palette = "Set1")+
  ggtitle(i)+xlab("Tumor VAF")+annotate("text",x=0.2,y=5,label=nrow(plot))

write.table(plot,file=paste(name,".txt",sep=""), append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=F)

plot=snp[which(snp$TUMOR_DP>9 & snp$NORMAL_DP>9),]
plot=plot[which(plot$NORMAL_X1=="0|0" & (plot$NORMAL_X10<0.001 | plot$NORMAL_X10 > 0.999)),]
plot=plot[which((plot$TUMOR_X10>0.05 & plot$TUMOR_X10<0.95)|((plot$TUMOR_X1=="1|1" | plot$TUMOR_X1=="0|1") & plot$TUMOR_X10>0.9)),]
q=ggplot(plot,aes(MP,fill=CALL))+geom_density(alpha=0.75,linetype="blank")+theme_bw()+
  scale_fill_brewer(palette = "Set1")+
  ggtitle(i)+coord_cartesian(ylim=c(0,1e+02))



ggsave(p,file=paste(name,"tumor.vaf.pdf",sep="."),width=8, height=8)
ggsave(q,file=paste(name,"MP.pdf",sep="."),width=8, height=8)
