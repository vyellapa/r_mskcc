library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)

anno=read.table("~/Desktop/record/MM_Malin/data_mutations_extended_AA_elli.tsv", header=T, sep="\t")

anno$Annotation="FILTER"
anno[anno$Normal.VAF>3,]$Annotation=paste(anno[anno$Normal.VAF>3,]$Annotation,"HIGH_N_VAF",sep=":")
anno[anno$Normal.VAF>3,]$Annotation

anno[, c(83:91)]=apply(anno[,c(83:91)],2,function(x) {gsub("[ACGT:]","",x)})
anno$max_1kg=apply(anno[,c(83:91)],1,function(x) {max(x, na.rm=T)})
anno$max_exac=apply(anno[,c(106:113)],1,function(x) {max(x, na.rm=T)})
anno$max_exac=gsub("-Inf","0",anno$max_exac)
anno[which(anno$max_1kg>0.001 | anno$max_exac >0.001),]$Annotation=paste(anno[which(anno$max_1kg>0.001 | anno$max_exac >0.001),]$Annotation,"GERMLINE_001",sep=":")
anno$Chromosome=paste("chr",anno$Chromosome, sep="")
anno$context=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, anno$Chromosome,anno$Start_Position-1,anno$Start_Position+1))




write.table(anno,file="~/Desktop/data_mutations_extended_AA_elli_anno_v2.txt", append=FALSE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)

library(ggplot2)
library(gridExtra)
anno=read.table("~/Desktop/MM_Malin/annotationsTejaMalin112116.txt", header=T, sep="\t",stringsAsFactors = FALSE)
anno$Classification=gsub(" ","",anno$Classification)
anno=anno[anno$Classification!="Ignore",]
anno[anno$Classification=="Possible",]$Classification="Unknown"
anno$Sample=rep("SAMPLE",nrow(anno))
anno[anno$Tumor_Sample_Barcode=="s_E_H_109074_T1_1_D1_1",]$Sample=rep("MGUS",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_109074_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112817_T1_1_D1_1",]$Sample=rep("Myeloma_S1",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112817_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112818_T1_1_D1_1",]$Sample=rep("Myeloma_S2",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112818_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112819_T1_1_D1_1",]$Sample=rep("Myeloma_S3",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112819_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112820_T1_1_D1_1",]$Sample=rep("Myeloma_S4",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112820_T1_1_D1_1",]))
anno[anno$Tumor_Sample_Barcode=="s_E_H_112821_T1_1_D1_1",]$Sample=rep("Myeloma_S5",nrow(anno[anno$Tumor_Sample_Barcode=="s_E_H_112821_T1_1_D1_1",]))



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
  xlab("")+
  ylab("Variant Allele Fraction")+ylim(c(0,100))+
  scale_color_brewer(palette = "Set1")+coord_flip()+xlim("Myeloma_S5", "Myeloma_S4", "Myeloma_S3", "Myeloma_S2", "Myeloma_S1","MGUS")

#anno$Samples= factor(anno$Samples, levels = c("MGUS", "Myeloma_S1", "Myeloma_S2", "Myeloma_S3", "Myeloma_S4", "Myeloma_S5"))
b=ggplot(anno, aes(Tumour.VAF))+geom_density(fill="blue", alpha=0.5,linetype="blank")+
  facet_wrap(~Sample,nrow=6)+
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
  ylab("")+
  scale_fill_brewer(palette = "Set1")+xlim(c(0,100))



grid.arrange(a,b, widths=c(1,1),heights=c(1,1),layout_matrix=rbind(c(1,2), c(1,2)))

