#Usage:Rscript MM_VAF_plotter.R annotationsTejaMalin112116.txt 
#Load libraries
library(ggplot2)
library(gridExtra)

setwd("/Users/yellapav/Desktop/record/MM_Malin/")
args <- commandArgs(trailingOnly = TRUE)
#file <- args[1] #Input File

#Read the file
#anno=read.table(file, header=T, sep="\t",stringsAsFactors = FALSE)
anno=read.table("/Users/yellapav/Desktop/record/MM_Malin/annotationsTejaMalin112116.txt", header=T, sep="\t",stringsAsFactors = FALSE)
#Remove white spaces and remove "Ignore" mutations
anno$Classification=gsub(" ","",anno$Classification)
anno=anno[anno$Classification!="Ignore",]
anno[anno$Classification=="Possible",]$Classification="Unknown"

#Rename sample names
anno$Sample=rep("SAMPLE",nrow(anno))
anno[anno$Tumor_Sample_Barcode=="s_E_H_109074_T1_1_D1_1",]$Sample=rep("MGUS",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_109074_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112817_T1_1_D1_1",]$Sample=rep("Myeloma_S1",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112817_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112818_T1_1_D1_1",]$Sample=rep("Myeloma_S2",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112818_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112819_T1_1_D1_1",]$Sample=rep("Myeloma_S3",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112819_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112820_T1_1_D1_1",]$Sample=rep("Myeloma_S4",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112820_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112821_T1_1_D1_1",]$Sample=rep("Myeloma_S5",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112821_T1_1_D1_1",]))


#Plot the box plot
a=ggplot(anno, aes(Sample,Tumour.VAF))+geom_boxplot()+geom_point(alpha=0.6,size=5,aes(color=Classification))+
  #facet_wrap(~Sample,nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="left", strip.background  = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  xlab("")+ylim(c(0,100))+
  ylab("Variant Allele Fraction")+
  scale_color_brewer(palette = "Set1")+coord_flip()+xlim("Myeloma_S4", "Myeloma_S1", "Myeloma_S2", "Myeloma_S3", "Myeloma_S5","MGUS")

#Plot the density plot

anno$carb <- factor(anno$Sample, levels=c("MGUS","Myeloma_S5", "Myeloma_S3", "Myeloma_S2","Myeloma_S1","Myeloma_S4"))

b=ggplot(anno, aes(Tumour.VAF))+geom_density(fill="blue", alpha=0.5,linetype="blank")+
  facet_wrap(~carb,nrow=6)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="blank", strip.background  = element_blank(),
        legend.title=element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  xlab("Variant Allele Fraction")+
  ylab("")+xlim(c(0,100))+
  scale_fill_brewer(palette = "Set1")


#Arrange the plots
c=grid.arrange(a,b, widths=c(1,1),heights=c(1,1),layout_matrix=rbind(c(1,2), c(1,2)))

#print the plot to file
cat("Writing plot to MM_VAF.pdf in the current directory......\n")
ggsave(c, file="MM_VAF_Apr5.pdf", width=12, height=10)


plot=read.table("/Users/yellapav/Desktop/record/malin_interview/plot", header=T, sep="\t",stringsAsFactors = FALSE)
ggplot(plot,aes(x=reorder(Gene,-Frequency),y=Frequency,fill="blue"))+geom_bar(stat="identity")+
  theme_bw()+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("Frequency %")+scale_y_continuous( limits=c(0,30))

plot=read.table("/Users/yellapav/Desktop/record/malin_interview/plot_top100", header=T, sep="\t",stringsAsFactors = FALSE)
ggplot(plot,aes(x=reorder(Gene,-Frequency),y=Frequency,fill="blue"))+geom_bar(stat="identity")+
  theme_bw(base_size = 15)+ 
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", strip.background  = element_blank(),legend.title=element_blank(), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey80"),
        panel.grid.major.x = element_blank())+
  xlab("")+ylab("Frequency %")+scale_y_continuous( limits=c(0,30))

