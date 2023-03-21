setwd("/Users/yellapav/Desktop/p212/uploads/oncoprint/")
library(ComplexHeatmap)


mat=read.table("p212_plot.txt", header = TRUE, stringsAsFactors = F, sep="\t")
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#4575b4", col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#d73027", col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#f46d43", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#B8DC3C", col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#1a9850", col = NA))
  }, 
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  Frame_Shift_Delet = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#EA3599", col = NA))
  }
  
)


anno=read.table("anno.txt", header=T)

ha_column = HeatmapAnnotation(df = anno, col = list(type1 = c("Black" =  "#33a02c", "White" = "#b2df8a")))

col = c("Frame_Shift_Ins" = "#f46d43", "Nonsense_Mutation" = "#d73027", "Del" = "#4575b4", "Splice_Site"="#B8DC3C", "Missense_Mutation" = "#1a9850", "Frame_Shift_Del"="#a6cee3" )

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, row_order=c("KRAS","BRAF","NRAS","DIS3","TP53","FAM46C","FGFR3","KDM6A","MAX","RASA2","SP140","KLHL6"),
          column_title = "", show_column_names = FALSE, bottom_annotation = ha_column, column_order=NULL,
          heatmap_legend_param = list(title = "Alternations", at = c("Nonsense_Mutation", "Del", "Frame_Shift_Ins", "Splice_Site", "Missense_Mutation", "Frame_Shift_Del"), 
                                      labels = c("Nonsense Mutation", "Deletion", "Frame Shift Ins", "Splice Site", "Missense", "Frame Shift Del")))


j=c("E-H-131112-T1-1-D1-1","E-H-131113-T1-1-D1-1","E-H-131114-T1-1-D1-1","E-H-131115-T1-1-D1-1","E-H-131121-T1-1-D1-1","E-H-131124-T1-1-D1-1","E-H-131126-T1-1-D1-1","E-H-131128-T1-1-D1-1","E-H-131140-T1-1-D1-1","E-H-131145-T1-1-D1-1","E-H-131146-T1-1-D1-1","E-H-131081-T1-1-D1-1","E-H-131087-T1-1-D1-1","E-H-131095-T1-1-D1-1","E-H-131098-T1-1-D1-1","E-H-131102-T1-1-D1-1","E-H-131106-T1-1-D1-1","E-H-131130-T1-1-D1-1","E-H-131132-T1-1-D1-1","E-H-131134-T1-1-D1-1","E-H-131135-T1-1-D1-1","E-H-131142-T1-1-D1-1","E-H-131143-T1-1-D1-1","E-H-131147-T1-1-D1-1","E-H-131149-T1-1-D1-1","E-H-131150-T1-1-D1-1","E-H-131152-T1-1-D1-1","E-H-131157-T1-1-D1-1")


setwd("/Users/yellapav/Desktop/p244/QC/picard_qc")
s=read.table("cov.plot",header=T)
p=read.table("pct.plot",header=T)

c=read.table("../concord.plot",header=T)
cc=read.table("../cont.plot",header=T)

m=melt(s)

mm=melt(p)
ggplot(m,aes(SAMPLE,value,fill=variable))+geom_bar(stat = "identity", position = position_dodge())+
  theme_bw(base_size = 10)+ 
  scale_fill_brewer(palette = "Blues")+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")

ggplot(mm,aes(SAMPLE,value,fill=variable))+geom_bar(stat = "identity", position = position_dodge())+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")


ggplot(c, aes(Sample, Concordance))+geom_bar(stat = "identity")+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")


ggplot(cc, aes(Sample, value,fill=Cont))+geom_bar(stat = "identity", position = position_dodge())+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")
