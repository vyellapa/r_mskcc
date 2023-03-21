library(ggplot2)
library(reshape)
library(gridExtra)

setwd("/Users/yellapav/Desktop/SOHN_163/VCF_MAF")

####################################
###      Coverage Plots       ######
####################################
med=read.table("/Users/yellapav/Desktop/SOHN_163/metrics/med", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
med$SAMPLE=gsub("-1-D1-1","", med$SAMPLE)
med$SAMPLE=gsub("-H-","", med$SAMPLE)

metrics=read.table("/Users/yellapav/Desktop/SOHN_163/metrics/metrics", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
metrics$SAMPLE=gsub("-1-D1-1","", metrics$SAMPLE)
metrics$SAMPLE=gsub("-H-","", metrics$SAMPLE)

q1=ggplot(med, aes(SAMPLE, MEDIAN_COVERAGE,fill="2"))+
  geom_bar(position = "stack", stat="identity", alpha=0.85)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(), legend.position="none", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  xlab("")+ylab("Median Coverage")+scale_y_continuous(breaks=seq(0, 275, 50))

m=melt(metrics)
mm=m[m$variable %in% c("PCT_DUPLICATES","PCT_BELOW_MAPQ20","PCT_BELOW_BASEQ20"),]
m=m[!(m$variable %in% c("PCT_DUPLICATES","PCT_BELOW_MAPQ20","PCT_BELOW_BASEQ20")),]
q2=ggplot(m, aes(SAMPLE, value,fill=variable))+
  geom_bar(position = "dodge", stat="identity", alpha=0.85, width=0.5)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(), legend.position="top", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("Percent")+scale_y_continuous(breaks=seq(0, 100, 20))

q3=ggplot(mm, aes(SAMPLE, value,fill=variable))+
  geom_bar(position = "dodge", stat="identity", alpha=0.85)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("Percent")

grid.arrange(q1,q2,q3)


####################################
###   Germline VAF Desnity    ######
####################################
vaf=read.table("plot_gm_099.tsv", sep="\t", header=F, stringsAsFactors = FALSE, quote="~")

gmD=ggplot(vaf,aes(V5,fill=V2))+geom_density(alpha=0.5, linetype="blank")+
  theme_bw()+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("VAF by Chromosome")+ylab("Density")




##################################################
######       FUNCTIONS                ############
##################################################

getDf <- function(maf) {
  muts_list=lapply(unique(maf$Variant_Classification),function(x) {
    n=nrow(maf[maf$Variant_Classification==x,]);
    s=maf[1,]$Tumor_Sample_Barcode; 
    print(c(s,x,n)) } )
  
  muts_num=as.data.frame(matrix(unlist(muts_list), nrow=length(unlist(muts_list))/3, byrow=T))
  
  colnames(muts_num)=c("Sample","Variant.Classification","Number")
  return(muts_num)
}



getDfType <- function(maf) {
  muts_list=lapply(unique(maf$Variant_Type),function(x) {
    n=nrow(maf[maf$Variant_Type==x,]);
    s=maf[1,]$Tumor_Sample_Barcode; 
    print(c(s,x,n)) } )
  
  muts_num=as.data.frame(matrix(unlist(muts_list), nrow=length(unlist(muts_list))/3, byrow=T))
  
  colnames(muts_num)=c("Sample","Variant.Classification","Number")
  return(muts_num)
}


indelInfo<- function(maf) {
  temp=maf
  
  ta1=temp[which(temp$Tumor_Seq_Allele2=="-" | temp$Tumor_Seq_Allele1=="-"),]
  ta2=temp[which(temp$Tumor_Seq_Allele2!="-" & temp$Tumor_Seq_Allele1!="-"),]
  ta1$LENGTH=(nchar(ta1$Tumor_Seq_Allele1)+nchar(ta1$Tumor_Seq_Allele2)-1)
  ta2$LENGTH=abs(nchar(ta2$Tumor_Seq_Allele1)-nchar(ta2$Tumor_Seq_Allele2))
  ta2$Variant_Type="COMPLEX"
  temp=rbind(ta1,ta2)
  return(temp)
  
}

plotIndel<-function(maf1,maf2) {
  temp1=indelInfo(maf1)
  colnames(temp1)
  temp2=indelInfo(maf2)
  temp=rbind(temp1,temp2)
  return(temp)
}



#Read MAFs and combine them
caveman1=read.table("E-H-109096-T1-1-D1-1.plus.strelka.rc.maf", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
pindel1=read.table("E-H-109096-T1-1-D1-1.v2.pindel.gm.maf", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
caveman2=read.table("E-H-109096-T2-1-D1-1.plus.strelka.rc.maf", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
pindel2=read.table("E-H-109096-T2-1-D1-1.v2.pindel.gm.maf", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")

setwd("/Users/yellapav/Desktop/SOHN_163/report_mar17/SOHN_ALL_VCFs_Mar282017/hiconf_sv/")
i="E-H-109096-T*-1-D1-1.filtered.cnv.annot.bedpe"

  p=paste(sep="","cut -f11-12 ",i,"|sed \'s/#\ chr/chr/g\'|grep -v \"^sample\"|sed 's/#\ chr/chr/g'|awk -F'\t' '{split($1,a,\",\"); print a[1],$2; }'|sort |uniq -c|awk '{print $1,$2,$3}' > temp")
  system(p)
sv_fil=read.table("temp", sep=" ", header=F, stringsAsFactors = FALSE, quote="~")
colnames(sv_fil)=c("count","Sample","SV_TYPE")

pindel1=pindel1[,c(1:123)]
pindel2=pindel2[,c(1:123)]

pindel1$VAF="NA"
#pindel1$E.H.109099.T2.1.D1.1="NA"

pindel2$VAF="NA"
#pindel2$E.H.109099.T1.1.D1.1="NA"

caveman1[,c("Variant_Type")]="SNV"
caveman2[,c("Variant_Type")]="SNV"
caveman1=caveman1[,c(1:124)]
caveman2=caveman2[,c(1:124)]

cave_pindel1=rbind(caveman1, pindel1)
cave_pindel2=rbind(caveman2, pindel2)
cave_pindel=rbind(cave_pindel1, cave_pindel2)



################################
# Plot 1 Number of Mutations 
################################
NS=c("Missense_Mutation","3'UTR", "5'UTR","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation")
muts_num1=getDf(cave_pindel1)
muts_num2=getDf(cave_pindel2)

muts_num=rbind(muts_num1, muts_num2)
ns_num=muts_num[muts_num$Variant.Classification %in% NS, ]

ns_num$Number=as.numeric(as.character(ns_num$Number))
muts_num$Number=as.numeric(as.character(muts_num$Number))
muts_num[muts_num$Variant.Classification=="RNA",]$Variant.Classification="Silent"
p1=ggplot(ns_num, aes(Sample, Number, label=Number,fill=Variant.Classification))+
  geom_bar(position = "stack", stat="identity", alpha=0.85)+
  geom_text(size = 4, position = position_stack(vjust = 0.5))+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")

pall1=ggplot(muts_num, aes(Sample, Number, fill=Variant.Classification))+
  geom_bar(position = "stack", stat="identity", alpha=0.95)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set3")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")


################################
## Plot 2 Indels
###############################
temp=plotIndel(pindel1,pindel2)

p2=ggplot(temp,aes(LENGTH, fill=Variant_Type))+geom_bar(stat="bin", binwidth=1, alpha=0.8)+
  coord_cartesian(xlim=c(0,50), ylim=c(0,200))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(),legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  scale_fill_brewer(palette = "Set1")+facet_wrap(~Tumor_Sample_Barcode)+ylab("")

#grid.arrange(p1,p2 ,ncol=2, nrow =1)

################################
## Plot 3 Raw Numbers
###############################
all=rbind(temp[,-c(125)], caveman1, caveman2)
samples=unique(temp$Tumor_Sample_Barcode)

s1=all[all$Tumor_Sample_Barcode==samples[1], ]
s2=all[all$Tumor_Sample_Barcode==samples[2], ]
ss1=as.data.frame(table(s1$Variant_Type))
ss1$Sample=rep(samples[1],nrow(ss1))
ss2=as.data.frame(table(s2$Variant_Type))
ss2$Sample=rep(samples[2],nrow(ss2))

ss=rbind(ss1,ss2)
colnames(ss)=c("Type", "Frequency","Sample")
p3=ggplot(ss,aes(Sample,Frequency,label=Frequency,fill=Type))+geom_bar(stat="identity", alpha=0.8)+
  geom_text(size = 4, position = position_stack(vjust = 0.5))+
theme_bw()+ 
scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")


################################
## Plot 4 VAF Density
###############################

caveman=rbind(caveman1, caveman2)
cavemanSub=caveman[, c("Tumor_Sample_Barcode","VAF")]
m=melt(cavemanSub)
p4=ggplot(m,aes(value, fill=variable))+geom_density(alpha=0.5, linetype="blank")+
  facet_wrap(~Tumor_Sample_Barcode,nrow=2)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("Variant Allele Fraction")+ylab("")


################################
## Plot 4 VAF Density
################################
sv=ggplot(sv_fil,aes(Sample,count,fill=SV_TYPE))+geom_bar(stat="identity")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")


grid.arrange(p1,p3,p4,p2, widths=c(1,2),layout_matrix=rbind(c(1,3), c(2,3),c(4,4)))
grid.arrange(p3,p1,p4,p2, widths=c(1,1,3),layout_matrix=rbind(c(1,2,3), c(1,2,3),c(4,4,4)))
grid.arrange(p3,pall1,p1,p4,p2, widths=c(0.75,0.75,0.75,2), heights=c(1,1,1.5),layout_matrix=rbind(c(1,2,3,4), c(1,2,3,4),c(5,5,5,5)))


sname=all$Tumor_Sample_Barcode[1]

bar=grid.arrange(p3,pall1,p1,sv, widths=c(1,1,1,1), layout_matrix=rbind(c(1,2,3,4)))
p4
p2
 
setwd("/Users/yellapav/Desktop/SOHN_163/report_mar17")
ggsave(bar, file=paste(sname,"_bar.pdf",sep=""), width=16, height=9)
ggsave(p4, file=paste(sname,"_vaf.pdf",sep=""), width=7, height=5)
ggsave(p2, file=paste(sname,"_indel_len.pdf",sep=""), width=7, height=5)






