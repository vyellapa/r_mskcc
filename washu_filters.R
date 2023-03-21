library(reshape2)
library("stringr")
library(plyr)
setwd("~/Desktop")
z="E-H-109096-T1-1-D1-1_comsic.vcf"

name=strsplit(z,"_",fixed=TRUE)[[1]][1]
snp=read.table(z,sep="\t",header=F,stringsAsFactors = F)


pat<-'MP=[0-9,e,\\-,\\+,\\.]+'
snp$MP=gsub("MP=","",str_extract(snp$V8, pat))

foo <- data.frame(do.call('rbind', strsplit(as.character(snp$V10),':',fixed=TRUE)))
colnames(foo)=paste("NORMAL_",colnames(foo),sep="")
snp=cbind(snp,foo)
foo <- data.frame(do.call('rbind', strsplit(as.character(snp$V11),':',fixed=TRUE)))
colnames(foo)=paste("TUMOR_",colnames(foo),sep="")
snp=cbind(snp,foo)
snp$CALL="SNP"



snp$NORMAL_DP=as.integer(as.character(snp$NORMAL_X2))+as.integer(as.character(snp$NORMAL_X3))+as.integer(as.character(snp$NORMAL_X4))+as.integer(as.character(snp$NORMAL_X5))+as.integer(as.character(snp$NORMAL_X6))+as.integer(as.character(snp$NORMAL_X7))+as.integer(as.character(snp$NORMAL_X8))+as.integer(as.character(snp$NORMAL_X9))
snp$TUMOR_DP=as.integer(as.character(snp$TUMOR_X2))+as.integer(as.character(snp$TUMOR_X3))+as.integer(as.character(snp$TUMOR_X4))+as.integer(as.character(snp$TUMOR_X5))+as.integer(as.character(snp$TUMOR_X6))+as.integer(as.character(snp$TUMOR_X7))+as.integer(as.character(snp$TUMOR_X8))+as.integer(as.character(snp$TUMOR_X9))
snp$TUMOR_X10=as.numeric(as.character(snp$TUMOR_X10))
snp$NORMAL_X10=as.numeric(as.character(snp$NORMAL_X10))
snp$MP=as.numeric(as.character(snp$MP))

snp$TUMOR_X2=as.numeric(as.character(snp$TUMOR_X2))
snp$TUMOR_X3=as.numeric(as.character(snp$TUMOR_X3))
snp$TUMOR_X4=as.numeric(as.character(snp$TUMOR_X4))
snp$TUMOR_X5=as.numeric(as.character(snp$TUMOR_X5))
snp$TUMOR_X6=as.numeric(as.character(snp$TUMOR_X6))
snp$TUMOR_X7=as.numeric(as.character(snp$TUMOR_X7))
snp$TUMOR_X8=as.numeric(as.character(snp$TUMOR_X8))
snp$TUMOR_X9=as.numeric(as.character(snp$TUMOR_X9))
snp$TUMOR_X10=as.numeric(as.character(snp$TUMOR_X10))


snp$A_DP=snp$TUMOR_X2+snp$TUMOR_X6
snp$C_DP=snp$TUMOR_X3+snp$TUMOR_X7
snp$G_DP=snp$TUMOR_X4+snp$TUMOR_X8
snp$T_DP=snp$TUMOR_X5+snp$TUMOR_X9

snp$NORMAL_DP=as.numeric(as.character(snp$NORMAL_DP))
snp$TUMOR_DP=as.numeric(as.character(snp$TUMOR_DP))


snp_ids=snp[snp$CALL=="SNP",]

#snp_ids=snp_ids[snp_ids$NORMAL_DP>9 & snp_ids$TUMOR_DP>10 & snp_ids$NORMAL_X10<1 & snp_ids$TUMOR_X10>0.045,]


########### 1) SUBSET 1 #########################
subset1=snp_ids[snp_ids$NORMAL_DP>7 & snp_ids$TUMOR_DP>10 & snp_ids$NORMAL_X10<0.01 & snp_ids$TUMOR_X10>0.04,]


################## 4) REF-REF Calls Gunes ###############################
subset3=snp_ids[as.numeric(as.character(snp_ids$NORMAL_DP))>7,]
subset3=subset3[as.numeric(as.character(subset3$TUMOR_DP))>10,]
subset3=subset3[which(subset3$V4==subset3$V5 & subset3$NORMAL_X10>0.995  & subset3$TUMOR_X10<0.995),]


if(nrow(subset3)>0 & nrow(subset1)>0) {out=rbind(subset3,subset1);}
else if(nrow(subset3)>0) {out=subset3;}
else {out=subset1;}

write.table(out,file=paste(name,"_post.txt",sep=""), append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)





ggplot(a,aes(V2,fill=factor(V1)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 10)
ggplot(a,aes(V14,fill=factor(V1)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(V9,fill=factor(V1)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,0.1))
ggplot(a,aes(V2,fill=factor(V1)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 0)

a=read.table("~/Desktop/plot2",header=T,sep = " ",stringsAsFactors = F)
a$MMQS_DIFF=abs(a$AVG_MMQS-a$AVG_MMQS.1)
a$RLEN_DIFF=abs(a$AVG_CLIPPED_LENGTH.1-a$AVG_CLIPPED_LENGTH)
a$MAPQ_DIFF=abs(a$MAPQ.1-a$MAPQ)

ggplot(aa,aes(WASHU,MAPQ,fill=factor(CGP)))+geom_point(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(50,60),ylim = c(0,7))+geom_vline(xintercept = 32)+geom_vline(xintercept = 55)

ggplot(a,aes(MAPQ,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(50,60),ylim = c(0,7))+geom_vline(xintercept = 32)+geom_vline(xintercept = 55)
ggplot(a,aes(COUNT,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,250))
ggplot(a,aes(BASEQ,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 32)
ggplot(a,aes(SE_MQ,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(POS_STRAND,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(NEG_STRAND,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(POS_AS_FRAC,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 0.3)
ggplot(a,aes(NUM_MM_FRAC,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,0.1))+geom_vline(xintercept = 0.038)
ggplot(a,aes(AVG_MMQS,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,150))+geom_vline(xintercept = 65)
ggplot(a,aes(Q2_READS,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,15))
ggplot(a,aes(AVG_DIST_Q2,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 0.32)
ggplot(a,aes(AVG_CLIPPED_LENGTH,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(120,150))+geom_vline(xintercept = 145)
ggplot(a,aes(DIST_3P,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 0.3)

ggplot(a,aes(MAPQ.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(55,60))
ggplot(a,aes(COUNT.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(1,10))
ggplot(a,aes(BASEQ.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 32)
ggplot(a,aes(SE_MQ.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(POS_STRAND.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(NEG_STRAND.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(POS_AS_FRAC.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 0.3)
ggplot(a,aes(NUM_MM_FRAC.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,0.1))+geom_vline(xintercept = 0.038)
ggplot(a,aes(AVG_MMQS.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,150))+geom_vline(xintercept = 65)
ggplot(a,aes(Q2_READS.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,15))
ggplot(a,aes(AVG_DIST_Q2.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+geom_vline(xintercept = 0.32)
ggplot(a,aes(AVG_CLIPPED_LENGTH.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(120,150))
ggplot(a,aes(DIST_3P.1,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)
ggplot(a,aes(MMQS_DIFF,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,150))+geom_vline(xintercept = 60)
ggplot(a,aes(MAPQ_DIFF,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,10))+geom_vline(xintercept = 3)
ggplot(a,aes(RLEN_DIFF,fill=factor(CGP)))+geom_density(linetype="blank",alpha=0.5)+coord_cartesian(xlim=c(0,50))+geom_vline(xintercept = 5)




### RESCUE Elli NYGC COSMIC compare 
setwd("/Users/yellapav/Desktop/elli")
library(reshape2)

plot=read.table("plot",header=T,sep=" ",stringsAsFactors = F)
m=melt(plot)
ggplot(m,aes(Sample,value,fill=factor(variable)))+geom_bar(stat="identity",position = position_dodge())+
theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
      panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+ggtitle("Calls Added uniquly by Strelka and Overlap with NYGC")


cosmic=read.table("cosmic.plot", header=T, sep=" ",stringsAsFactors = F)
m=melt(cosmic)
ggplot(m,aes(Sample,value,fill=factor(variable)))+geom_bar(stat="identity",position = position_dodge())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+ggtitle("Calls Added uniquly by Strelka and Overlap with COSMIC")


rplot=read.table("rescue.plot",header=T,sep=" ",stringsAsFactors = F)
m=melt(rplot)
ggplot(m,aes(Sample,value,fill=factor(variable)))+geom_bar(stat="identity",position = position_dodge())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+ggtitle("Calls Rescued from SNPs and Overlap with NYGC 2Plus")

rcosmic=read.table("cosmic_rescue.plot",header=T,sep=" ",stringsAsFactors = F)
m=melt(rcosmic)
ggplot(m,aes(Sample,value,fill=factor(variable)))+geom_bar(stat="identity",position = position_dodge())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")+ggtitle("Calls Rescued from SNPs and Overlap with Cosmic")



#less ../SOHN_163/VCF_MAF/E-H-109*.filter.caveman.rc.gm.maf|awk -F'\t' '{OFS="\t";if($16~/-T1-/) {$7="CGP";print $0}}'|cut -f16-28,124 > allcgp 
#less ../SOHN_163/VCF_MAF/E-H-109*.filter.caveman.rc.gm.maf|awk -F'\t' '{OFS="\t";if($16~/-T2-/) {$7="CGP";print $0}}'|cut -f16-28,125 >> allcgp 
vaf=read.table("allvaf",header=F,sep="\t",stringsAsFactors = F)
rvaf=read.table("rescuevaf",header=F,sep="\t",stringsAsFactors = F)
allcgp=read.table("allcgp",header=F,sep="\t",stringsAsFactors = F)
head(vaf)
vaf$V1=gsub(".rc","",as.character(vaf$V1))
vaf$V2="ADDED"
plot=rbind(rvaf,vaf,allcgp)

ggplot(plot,aes(V14, fill=factor(V2)))+geom_density(linetype="blank", alpha=0.5)+facet_wrap(~V1,scales="free")+xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+coord_cartesian(xlim=c(0,1))


ggplot(plot,aes(V13,fill=factor(V2)))+geom_density(linetype="blank", alpha=0.5)+facet_wrap(~V1,scales="free")+xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+coord_cartesian(xlim=c(0,500))




setwd("~/Desktop/elli")
vcf_files=c("E-H-109096-T1-1-D1-1_sig.vcf","E-H-109096-T2-1-D1-1_sig.vcf","E-H-109097-T1-1-D1-1_sig.vcf","E-H-109097-T2-1-D1-1_sig.vcf","E-H-109098-T1-1-D1-1_sig.vcf","E-H-109098-T2-1-D1-1_sig.vcf","E-H-109099-T1-1-D1-1_sig.vcf","E-H-109099-T2-1-D1-1_sig.vcf","E-H-109100-T1-1-D1-1_sig.vcf","E-H-109100-T2-1-D1-1_sig.vcf","E-H-109101-T1-1-D1-1_sig.vcf","E-H-109101-T2-1-D1-1_sig.vcf")
sample_names=c("E-H-109096-T1","E-H-109096-T2","E-H-109097-T1","E-H-109097-T2","E-H-109098-T1","E-H-109098-T2","E-H-109099-T1","E-H-109099-T2","E-H-109100-T1","E-H-109100-T2","E-H-109101-T1","E-H-109101-T2")


ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = "hg19")

#vcf=read_vcfs_as_granges("c.vcf", "sample_names", genome = "hg19")

auto = extractSeqlevelsByGroup(species="Homo_sapiens",
                               style="UCSC",
                               group="auto")
vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
mutation_types(vcfs)
type_occurrences = mut_type_occurrences(vcfs, ref_genome)

test_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(test_matrix)