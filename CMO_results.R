setwd("~/Desktop/MM_Malin/")
library(ggplot2)
library(gridExtra)
library(reshape)
maf=read.table("data_mutations_extended.txt", sep="\t", header=T,quote="~")
maf=maf[maf$Validation_Status!="REDACTED",]
maf$Tumor_Sample_Barcode=gsub("s_","",maf$Tumor_Sample_Barcode)
maf$vaf=maf$t_alt_count/(maf$t_alt_count+maf$t_ref_count)
maf$tumor_depth=(maf$t_alt_count+maf$t_ref_count)
maf$normal_depth=(maf$n_alt_count+maf$n_ref_count)

ggplot(maf,aes(Tumor_Sample_Barcode,fill="1"))+geom_bar(stat="count",alpha=0.8)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="blank", strip.background  = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  xlab("")+
  ylab("# of Mutations")+
  scale_fill_brewer(palette = "Set1")

a=ggplot(maf, aes("1",vaf))+geom_boxplot()+geom_point(alpha=0.3,size=5,color="blue")+facet_wrap(~Tumor_Sample_Barcode,nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="blank", strip.background  = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  xlab("")+
  ylab("VAF")+
  scale_fill_brewer(palette = "Set1")

b=ggplot(maf, aes(vaf))+geom_density(fill="blue", alpha=0.5,line_type="blank")+facet_wrap(~Tumor_Sample_Barcode,nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="blank", strip.background  = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  xlab("")+
  ylab("VAF")+
  scale_fill_brewer(palette = "Set1")


grid.arrange(b,a, widths=c(1),heights=c(1,1),layout_matrix=rbind(c(1), c(2)))

m=melt(mafd)
ggplot(m, aes(value,fill=variable))+geom_density(alpha=0.5)+
  facet_wrap(~Tumor_Sample_Barcode,scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="blank", strip.background  = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(), 
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.y = element_blank())+
  xlab("Depth")+
  ylab("")+
  scale_fill_brewer(palette = "Set1")
