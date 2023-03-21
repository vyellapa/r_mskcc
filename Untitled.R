library(reshape2)
library(ggplot2)
library(viridis)
#myType ASH_2018
a=read.table("/Users/yellapav/Desktop/mytype_ash_2018/strelka/tgd_snv.maf.txt",header = T,sep = "\t", stringsAsFactors = FALSE, quote = "~")
a=read.table("/Users/yellapav/Desktop/mytype_ash_2018/strelka/wgs_snv.maf.txt",header = T,sep = "\t", stringsAsFactors = FALSE, quote = "~")
ut=read.table("/Users/yellapav/Desktop/mytype_ash_2018/strelka/unique_tgd",header = F,sep = "\t", stringsAsFactors = FALSE, quote = "~")
uw=read.table("/Users/yellapav/Desktop/mytype_ash_2018/strelka/unique_wgs",header = F,sep = "\t", stringsAsFactors = FALSE, quote = "~")
both=read.table("/Users/yellapav/Desktop/mytype_ash_2018/strelka/both",header = F,sep = "\t", stringsAsFactors = FALSE, quote = "~")


a$key=paste(a$Chromosome, a$Start_Position,a$End_Position,a$Reference_Allele,a$Tumor_Seq_Allele1,a$Tumor_Seq_Allele2,a$Tumor_Sample_Barcode,sep = "_")
a$key=gsub("-D1-1","",a$key)
a$key=gsub("-D1-2","",a$key)

m=merge(both,a, by.x="V1", by.y = "key")
m$ann=paste(m$Hugo_Symbol, m$HGVSp_Short,sep="_")
m=(m[,c(17:33)])
m$max.val=apply(m[,c(1:16)], 1, max)
m=m[m$max.val>=0.02,]
mm=melt(m)
mm$value=round((mm$value*100), digits = 1)

ggplot(mm, aes(variable, ann))+geom_tile(aes(fill=value))+xlab("")+ylab("")+geom_text(data=subset(mm, value > 0), aes(label=value))+
  theme_bw()+scale_fill_gradient2(low="#ffffff", high="#023858", limits=c(0.0, 50))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())


mm$type="both"
mm[grep("wgs",mm$variable),]$type="WGS"
mm[grep("tgd",mm$variable),]$type="TGD"
mm$variable=gsub(".D1.1_tgd","",mm$variable)
mm$variable=gsub(".D1.2_wgs","",mm$variable)
mm=mm[mm$variable!="max.val",]

######################################
######## VAF Correlations ############
see=cor(m[,1:15], use = "complete.obs")
see=melt(see)
see$g=paste(see$Var1,see$Var2,sep="__")
see=see[grep("T.*T",see$g, perl = TRUE),]
see$s1=substr(see$g,1,15)
see$a1=substr(see$g,22,24)
see$s2=substr(see$g,27,41)
see$a2=substr(see$g,48,50)
summary(see[which(see$s1==see$s2 & see$a1!=see$a2),]$value)

##########################################
tmaf=read.table("/Users/yellapav/Desktop/mytype_ash_2018/strelka/tgd.all.maf.txt",header = T,sep = "\t", stringsAsFactors = FALSE, quote = "~")
tmaf$max.val=as.numeric(apply(tmaf[,c(16:31)], 1, max))
tmaf=tmaf[tmaf$max.val>=0.02,]
tmaf=tmaf[!is.na(tmaf$Hugo_Symbol),]
summary(as.vector(table(tmaf$Tumor_Sample_Barcode)))


#######################################################################
######################  Purity           ##############################
#######################################################################

a=read.table("~/Desktop/mytype_ash_2018/purity.txt",sep = " ", header=T)

#######################################################################
#############  Coverage             ###################################
#######################################################################

hs=read.table("/Users/yellapav/Desktop/mytype_ash_2018/hs.plot", header=TRUE, sep=" ")
zzz=do.call(rbind.data.frame, (strsplit(as.character(hs$NAME), '[.]')))
colnames(zzz)=c("Sample","Type")
hs=cbind(hs,zzz)
summary(hs[hs$Type=="overall",c(2)])

### Mean Coverage
summary(hs[hs$Type=="cds",c(2)])
summary(hs[hs$Type=="finger",c(2)])
summary(hs[hs$Type=="focal",c(2)])
summary(hs[hs$Type=="igh",c(2)])

## 2X coverage
summary(hs[hs$Type=="cds",c(4)])
summary(hs[hs$Type=="finger",c(4)])
summary(hs[hs$Type=="focal",c(4)])
summary(hs[hs$Type=="igh",c(4)])

## 100X coverage
summary(hs[hs$Type=="cds",c(8)])
summary(hs[hs$Type=="finger",c(8)])
summary(hs[hs$Type=="focal",c(8)])
summary(hs[hs$Type=="igh",c(8)])

mm=read.table("~/Desktop/p220_wgs/misc/oncoplot/VCF_MAF/I-H-130719-T1-6-D1-2_vs_I-H-130719-N1-1-D1-2.caveman.maf", sep="\t", header = T, stringsAsFactors = F, quote='~')
mm=mm[mm$Chromosome==6,]
mm=mm[order(mm$Start_Position),]
head(mm$Start_Position)
mm$dist=0
for(i in 2:nrow(mm)) {
  mm[i,]$dist=mm[i,"Start_Position"]-mm[(i-1),"Start_Position"]
}
mm$dist=log(mm$dist+1)
ggplot(mm,aes(Start_Position,dist))+geom_point()

#TP, FP, FN, TN
df=as.data.frame(rbind(c(670,202,74,640),c(6700,2020,740,6400)))
bb=numeric(0)
for(i in 1:nrow(df)) {
  df1=df[i,]
dat <- as.table(matrix(c(df1$V1,df1$V2,df1$V3,df1$V4), nrow = 2, byrow = TRUE))
colnames(dat) <- c("Dis+","Dis-")
rownames(dat) <- c("Test+","Test-")
rval <- epi.tests(dat, conf.level = 0.95)
aa=cbind(df1,as.data.frame(cbind(round(rval$elements$sensitivity,digits=2),round(rval$elements$specificity,digits=2))))
if(exists("bb")) {bb=rbind(bb,aa)} else {bb=aa}
}
bb=as.data.frame(bb)
colnames(bb)=c(as.vector(colnames(df)),c("sens","sens.lower","sens.upper","spec","spec.lower","spec.upper"))
bb
638/703
