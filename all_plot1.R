library(ggplot2)
library(reshape)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

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
  temp=temp[temp$Variant_Type!="SNV",]
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
combined=read.table("E-H-109099-T1-1-D1-1.plus.strelka.rc.maf", sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
colnames(combined)[124]="VAF"
combined=combined[,c(1:124)]
combined=caveman[which(caveman$Chromosome==27),]


for(i in list.files(pattern="E.*.plus.strelka.rc.maf")) {
  caveman=read.table(i, sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
  caveman[,c("Variant_Type")]="SNV"
  colnames(caveman)[124]="VAF"
  caveman=caveman[,c(1:124)]
  
  
  combined=rbind(combined,caveman)
  
}

for(i in list.files(pattern="E.*v2.pindel.gm.maf")) {
  pindel=read.table(i, sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
  pindel=pindel[,c(1:124)]
  colnames(pindel)[124]="VAF"
  
  combined=rbind(combined,pindel)
  
}


################################
# Plot 1 Number of Mutations 
################################
NS=c("Missense_Mutation","3'UTR", "5'UTR","Frame_Shift_Ins","Splice_Site","Frame_Shift_Del", "Nonsense_Mutation", "Nonstop_Mutation")
muts_num=table(combined[,c(16,9)])
muts_num=melt(muts_num)
ns_num=muts_num[muts_num$Variant_Classification %in% NS, ]

ns_num$value=as.numeric(as.character(ns_num$value))
muts_num$value=as.numeric(as.character(muts_num$value))

muts_num[muts_num$Variant_Classification=="RNA",]$Variant_Classification="Silent"

colnames(muts_num)=c("Sample","Variant.Classification","Number")
colnames(ns_num)=c("Sample","Variant.Classification","Number")


totals <- ns_num %>% group_by(Sample) %>% summarize(total = sum(Number))


p1=ggplot(ns_num, aes(Sample, Number, label=Number,fill=Variant.Classification))+
  geom_bar(position = "stack", stat="identity", alpha=0.85, vjust=-.1)+
  geom_text(aes(Sample, total, label = total, fill = NULL),vjust=-.05,color="#666666", data = totals)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")


totals <- muts_num %>% group_by(Sample) %>% summarize(total = sum(Number))
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
mycolors=colorRampPalette(brewer.pal(name="Accent", n = 9))(14)
pall1=ggplot(muts_num, aes(Sample, Number, fill=Variant.Classification))+
  geom_bar(position = "stack", stat="identity", alpha=0.95)+
  geom_text(aes(Sample, total, label = total, fill = NULL),vjust=-.05,color="#666666",data = totals)+
  theme_bw()+ 
  scale_fill_manual(values=mycolors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")


################################
## Plot 2 Indels
###############################

temp=indelInfo(combined)
#temp=plotIndel(pindel1,pindel2)

p2=ggplot(temp,aes(LENGTH, fill=Variant_Type))+geom_bar(stat="bin", binwidth=1, alpha=0.8)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(),legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  scale_fill_brewer(palette = "Set1")+facet_wrap(~Tumor_Sample_Barcode,scales="free")+ylab("")+coord_cartesian(xlim=c(0,50))

#grid.arrange(p1,p2 ,ncol=2, nrow =1)

################################
## Plot 3 Raw Numbers
###############################

counts=temp[,c(1:100)]
snv=combined[combined$Variant_Type=="SNV",c(1:100)]
counts=rbind(snv,counts)


p3=ggplot(counts,aes(Tumor_Sample_Barcode,fill=Variant_Type))+geom_bar(stat="count", alpha=0.8)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),axis.line.y = element_line(size = 0.001, colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")


################################
## Plot 4 VAF Density
###############################

caveman=combined[combined$Variant_Type=="SNV",]

p4=ggplot(caveman,aes(VAF, fill="1"))+geom_density(alpha=0.7, linetype="blank")+
  facet_wrap(~Tumor_Sample_Barcode,scales="free",nrow=3)+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+
  xlab("Variant Allele Fraction")+ylab("")



grid.arrange(p1,p3,p4,p2, widths=c(1,2),layout_matrix=rbind(c(1,3), c(2,3),c(4,4)))
grid.arrange(p3,p1,p4,p2, widths=c(1,1,3),layout_matrix=rbind(c(1,2,3), c(1,2,3),c(4,4,4)))
grid.arrange(p3,pall1,p1,p4,p2, widths=c(0.75,0.75,0.75,2), heights=c(1,1,1.5),layout_matrix=rbind(c(1,2,3,4), c(1,2,3,4),c(5,5,5,5)))









#####################################
## Venn Diagram 
#####################################

all_cave=rbind(caveman1,caveman2)
uniq=all_cave[!duplicated(all_cave[,c(5,6,7,11,13,16)]),]
uniq=uniq[,c(5,6,7,11,13)]

## Called in T1
venn=data.frame(sample="E099-T1")

venn$Called=(nrow(caveman1))-667

## Called in T2
venn$T2=nrow(caveman2)-667

## Called in Both
venn$Common=nrow(uniq[duplicated(uniq),])
venn$both2=nrow(uniq[duplicated(uniq),])

## Not Called but present in T1
venn$Present=nrow(caveman2[caveman2$E.H.109099.T1.1.D1.1>0.01,])-667

## Not Called but present in T2
venn$NC_T2=nrow(caveman1[caveman1$E.H.109099.T2.1.D1.1>0.01,])-667
m=melt(venn)
m$sample=as.character(m$sample)
m[m$variable=="T2",c("sample","variable")]=c("E099-T2","Called")
m[m$variable=="both2",c("sample", "variable")]=c("E099-T2","Common")
m[m$variable=="NC_T2",c("sample","variable")]=c("E099-T2","Present")

ggplot(m, aes(sample,value,fill=variable))+geom_bar(stat="identity")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")



############# SVs ###############

sv=read.table("/Users/yellapav/Desktop/SOHN_163/SV_CNV/bedpe.all",sep="\t",header=T,stringsAsFactors = F,quote = "~")
sv$sample=gsub(",E.*N1.*","",sv$sample)
ggplot(sv, aes(sample, fill=factor(svclass)))+geom_bar(stat = "count")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("")