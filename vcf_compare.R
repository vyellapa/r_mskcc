library(ggplot2)
library("stringr")


args=commandArgs(TRUE)
cn2vcf=args[1]
cn5vcf=args[2]
bam5=args[3]
bam2=bam5

#cn2vcf="~/Desktop/SOHN_163/one_offs/SNP/E-H-109099-T2-1-D1-1_vs_E-H-109099-N1-1-D1-1.caveman.muts.annot.vcf.gz"
#cn5vcf="~/Desktop/SOHN_163/one_offs/SNP/E-H-109099-T1-1-D1-1_vs_E-H-109099-N1-1-D1-1.caveman.muts.annot.vcf.gz"

CN2=read.table(cn2vcf,header=F,sep="\t",stringsAsFactors = F)
CN5=read.table(cn5vcf,header=F,sep="\t",stringsAsFactors = F)


pat<-"ASRD=[0-9,.]+"
CN2$m=paste(sep=":",CN2$V1,CN2$V2)
CN2$ASRD=gsub("ASRD=","",str_extract(CN2$V8, pat))
CN2$ASRD=as.numeric(as.character(CN2$ASRD))
cn2p=CN2[which(CN2$V7=="PASS" & CN2$ASRD>0.9),]

CN5$m=paste(sep=":",CN5$V1,CN5$V2)
CN5$ASRD=gsub("ASRD=","",str_extract(CN5$V8, pat))
CN5$ASRD=as.numeric(as.character(CN5$ASRD))
cn5p=CN2[which(CN5$V7=="PASS" & CN5$ASRD>0.9),]

common=cn5p[(cn5p$m %in% cn2p$m),]
common=common[,c(1,2,2,4,5)]
uniq_cn2=cn2p[!(cn2p$m %in% cn5p$m),]
uniq_cn2=uniq_cn2[,c(1,2,2,4,5)]
uniq_cn5=cn5p[!(cn5p$m %in% cn2p$m),]
uniq_cn5=uniq_cn5[,c(1,2,2,4,5)]



write.table(common,file="common.in", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=F)
write.table(uniq_cn2,file="uniq_cn2.in", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=F)
write.table(uniq_cn5,file="uniq_cn5.in", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=F)

command=paste(sep=" ","bam-readcount","-f","~/local/resources/hs37d5.fa","-l","common.in",bam5, "> common.out;paste common.in common.out |cut -f1-15 > common.out1")
system(command)
command=paste(sep=" ","bam-readcount","-f","~/local/resources/hs37d5.fa","-l","uniq_cn2.in",bam2, "> common.out;paste uniq_cn2.in common.out |cut -f1-15 > uniq_cn2.out1")
system(command)
command=paste(sep=" ","bam-readcount","-f","~/local/resources/hs37d5.fa","-l","uniq_cn5.in",bam5, "> common.out;paste uniq_cn5.in common.out |cut -f1-15 > uniq_cn5.out1")
system(command)

u.cn2=read.table("~/Desktop/uniq_cn2.out1",header=F,sep="\t",stringsAsFactors = F)
u.cn5=read.table("~/Desktop/uniq_cn5.out1",header=F,sep="\t",stringsAsFactors = F)
shared=read.table("~/Desktop/common.out1",header=F,sep="\t",stringsAsFactors = F)

u.cn2$STATUS="UNIQ_CN2"
u.cn5$STATUS="UNIQ_CN5"
shared$STATUS="SHARED"

all.rc=rbind(u.cn2,u.cn5,shared)

#Find out the REF base and readcount
all.rc$REF=all.rc$V11
all.rc[all.rc$V4=="C",]$REF=all.rc[all.rc$V4=="C",]$V12
all.rc[all.rc$V4=="G",]$REF=all.rc[all.rc$V4=="G",]$V13
all.rc[all.rc$V4=="T",]$REF=all.rc[all.rc$V4=="T",]$V14

#Find out the ALT base and readcount
all.rc$ALT=all.rc$V11
all.rc[all.rc$V5=="C",]$ALT=all.rc[all.rc$V5=="C",]$V12
all.rc[all.rc$V5=="G",]$ALT=all.rc[all.rc$V5=="G",]$V13
all.rc[all.rc$V5=="T",]$ALT=all.rc[all.rc$V5=="T",]$V14

all.rc=all.rc[,c(1,2,4,5,9,16,17,18)]
r=as.data.frame(do.call('rbind',strsplit(as.character(all.rc$REF), ':', fixed=TRUE)))
a=as.data.frame(do.call('rbind',strsplit(as.character(all.rc$ALT), ':', fixed=TRUE)))

all.rc=all.rc[,c(1,2,3,4,5,6)]
all.rc=cbind(all.rc,r,a)
colnames(all.rc)=c("CHR","POS","R1","A1","DP","CATEGORY","REF","COUNT","MAPQ","BASEQ","SE_MAPQ","POS_STRAND_READS","NEG_STRAND_READS","POS_AS_FRAC","MM_AS_FRAC","MMQS","NUM_Q2","AVG_DIST_Q2","CLIPPED_LEN","DIST_3P","ALT1","ALT_COUNT","ALT_MAPQ","ALT_BASEQ","ALT_SE_MAPQ","ALT_POS_STRAND_READS","ALT_NEG_STRAND_READS","ALT_POS_AS_FRAC","ALT_MM_AS_FRAC","ALT_MMQS","ALT_NUM_Q2","ALT_AVG_DIST_Q2","ALT_CLIPPED_LEN","ALT_DIST_3P")

all.rc$STRANDEDNESS=apply(all.rc,1,function(x) {p=as.numeric(as.character(x["ALT_POS_STRAND_READS"]));n=as.numeric(as.character(x["ALT_NEG_STRAND_READS"]));m=min(p,n)/(p+n); return(m)})
all.rc$VAF=apply(all.rc,1,function(x) {p=as.numeric(as.character(x["COUNT"]));n=as.numeric(as.character(x["ALT_COUNT"]));m=n/(p+n); return(m)})
all.rc$MMQS_DIFF=abs(as.numeric(as.character(all.rc$ALT_MMQS))-as.numeric(as.character(all.rc$MMQS)))

p=ggplot(a,aes(ALT_COUNT,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+coord_cartesian(xlim=c(0,150),ylim=c(0,0.15))+theme_bw()
ggsave(p, file="ALT_COUNT.pdf", width=6, height=4)

p=ggplot(a,aes(ALT_MAPQ,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(45,60),ylim=c(0,2))
ggsave(p, file="ALT_MAPQ.pdf", width=6, height=4)

p=ggplot(a,aes(ALT_BASEQ,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(ylim=c(0,0.4))
ggsave(p, file="ALT_BASEQ.pdf", width=6, height=4)

p=ggplot(a,aes(STRANDEDNESS,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggsave(p, file="STRANDEDNESS.pdf", width=6, height=4)

p=ggplot(a,aes(ALT_POS_AS_FRAC,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+geom_vline(xintercept = 0.01)
ggsave(p, file="ALT_POS_AS_FRAC.pdf", width=6, height=4)

p=ggplot(a,aes(ALT_MM_AS_FRAC,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(ylim=c(0,300),xlim=c(0.01,0.065))
ggsave(p, file="ALT_MM_AS_FRAC.pdf", width=6, height=4)

p=ggplot(a,aes(ALT_MMQS,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0,120),ylim=c(0,0.2))
ggsave(p, file="STRANDEDNESS.pdf", width=6, height=4)

p=ggplot(a,aes(ALT_AVG_DIST_Q2,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggsave(p, file="ALT_MMQS.pdf", width=6, height=4)

p=ggplot(a,aes(VAF,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggsave(p, file="VAF.pdf", width=6, height=4)

p=ggplot(a,aes(MMQS_DIFF,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0,100),ylim=c(0,0.15))
ggsave(p, file="MMQS_DIFF.pdf", width=6, height=4)

ggplot(a,aes(ALT_CLIPPED_LEN,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(125,150),ylim=c(0,0.75))
ggsave(p, file="ALT_CLIPPED_LEN.pdf", width=6, height=4)

ggplot(a,aes(ALT_DIST_3P,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()
ggsave(p, file="ALT_DIST_3P.pdf", width=6, height=4)

#ggplot(a,aes(MP,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0.80,1),ylim=c(0,25))
#ggsave(p, file="STRANDEDNESS.pdf", width=6, height=4)

ggplot(a,aes(DP,fill=CATEGORY))+geom_density(linetype="blank",alpha=0.6)+facet_wrap(~SAMPLE)+theme_bw()+coord_cartesian(xlim=c(0,750),ylim=c(0,0.015))
ggsave(p, file="DP.pdf", width=6, height=4)




