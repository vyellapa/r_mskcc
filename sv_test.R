setwd("/Users/yellapav/Desktop/p251/SV/")
library(ComplexHeatmap)
#mat=read.table("plot1", header = TRUE, stringsAsFactors = F, sep="\t")
mat=read.table("p251_SV_oncoplot.txt", header = TRUE, stringsAsFactors = F, sep="\t")
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EDEDED", col = NA))
  },
  TRA = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#377EB8", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E41A1C", col = NA))
  },
  Silent = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#4DAF4A", col = NA))
  },
  NS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#FF7F00", col = NA), name=c("hey"))
  },
  
  p.R683S = function(x, y, w, h) { grid.text("p.R683S",x, y, gp=gpar(fontsize=16))},
  p.A146T = function(x, y, w, h) { grid.text("p.A146T",x, y , gp=gpar(fontsize=16))},
  p.R683G = function(x, y, w, h) { grid.text("p.R683G",x, y, gp=gpar(fontsize=16)) },
  
  INV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#FFFF33", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#984EA3", col = NA))
  }
)

col = c("Silent" = "#4DAF4A", "AMP" = "#E41A1C", "TRA" = "#377EB8", "NS"="#FF7F00", "INV"="#FFFF33", "Fusion" = "#984EA3","NST"= "#984EA3", "p.R683S" = "#000000","p.R683G" = "#000000","p.A146T" = "#000000")



sampord=c("I-H-106917-T2-1","I-H-106917-T2-2","I-H-106917-T2-3","I-H-106917-T2-4","I-H-130718-T1-1","I-H-130718-T1-10","I-H-130718-T1-2","I-H-130718-T1-3","I-H-130718-T1-4","I-H-130718-T1-5","I-H-130718-T1-6","I-H-130718-T1-7","I-H-130718-T1-8","I-H-130718-T1-9","I-H-130719-T1-1","I-H-130719-T1-2","I-H-130719-T1-3","I-H-130719-T1-4","I-H-130719-T1-5","I-H-130719-T1-6","I-H-130720-T1-1","I-H-130720-T1-2","I-H-130720-T1-3","I-H-130720-T1-4","I-H-130720-T1-5","I-H-130720-T1-6","I-H-130720-T1-7","I-H-130720-T1-8")
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "IGH Rearrangements with a CGC partner", show_column_names = TRUE, column_order= NULL ,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "TRA", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Translocation", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))




oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL, show_pct = TRUE, pct_gp = gpar(fontsize=0), barplot_ignore = c("p.R683G","p.R683S","p.A146T"),
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))


mat=read.table("~/Desktop/SV_test/oncoplot_order_222.txt", header = TRUE, stringsAsFactors = F, sep="\t")
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)


alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EDEDED", col = NA))
  },
  TRA = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#377EB8", col = NA))
  }
)

col = c( "TRA" = "#377EB8" )



#sampord=c("I-H-106917-T2-1","I-H-106917-T2-2","I-H-106917-T2-3","I-H-106917-T2-4","I-H-130718-T1-1","I-H-130718-T1-10","I-H-130718-T1-2","I-H-130718-T1-3","I-H-130718-T1-4","I-H-130718-T1-5","I-H-130718-T1-6","I-H-130718-T1-7","I-H-130718-T1-8","I-H-130718-T1-9","I-H-130719-T1-1","I-H-130719-T1-2","I-H-130719-T1-3","I-H-130719-T1-4","I-H-130719-T1-5","I-H-130719-T1-6","I-H-130720-T1-1","I-H-130720-T1-2","I-H-130720-T1-3","I-H-130720-T1-4","I-H-130720-T1-5","I-H-130720-T1-6","I-H-130720-T1-7","I-H-130720-T1-8")
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "IGH Rearrangements identified using FISH and myType", show_column_names = FALSE ,row_order=NULL,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "TRA", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Translocation", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))




oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL, show_pct = TRUE, pct_gp = gpar(fontsize=0), barplot_ignore = c("p.R683G","p.R683S","p.A146T"),
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))




mat=read.table("~/Desktop/p222_fish_seq_transpose.txt", header = TRUE, stringsAsFactors = F, sep="\t")
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)
mat_nonas=gsub("NAS;","",mat)


alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EDEDED", col = NA))
  },
  Translocation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#377EB8", col = NA))
  },
  NAS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#000600", col = NA))
  }
)

col = c( "Translocation" = "#377EB8", "NAS" = "#000600" )

rows=c("FISHIGH_Rearr_ANY","FISHt(4.14)","t(14:4)_Delly","t(4:14)_BRASS","FISHt(6.14)","t(14:6)_Delly","t(6:14)_BRASS","FISHt(11.14)","t(14:11)_Delly","t(11:14)_BRASS","FISHt(14.16)","t(16:14)_Delly","t(14:16)_BRASS","FISHt(14.20) ","t(20:14)_Delly","t(14:20)_BRASS","t(14:8)_Delly","t(8:14)_BRASS")

#sampord=c("I-H-106917-T2-1","I-H-106917-T2-2","I-H-106917-T2-3","I-H-106917-T2-4","I-H-130718-T1-1","I-H-130718-T1-10","I-H-130718-T1-2","I-H-130718-T1-3","I-H-130718-T1-4","I-H-130718-T1-5","I-H-130718-T1-6","I-H-130718-T1-7","I-H-130718-T1-8","I-H-130718-T1-9","I-H-130719-T1-1","I-H-130719-T1-2","I-H-130719-T1-3","I-H-130719-T1-4","I-H-130719-T1-5","I-H-130719-T1-6","I-H-130720-T1-1","I-H-130720-T1-2","I-H-130720-T1-3","I-H-130720-T1-4","I-H-130720-T1-5","I-H-130720-T1-6","I-H-130720-T1-7","I-H-130720-T1-8")
oncoPrint(mat_nonas, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "IGH Rearrangements identified using FISH and myType", show_column_names = FALSE ,row_order=rows,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "Translocation", "NAS", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Translocation", "Not Available", "Non Synonymous Mutation",  "Inversion", "Fusion")))




oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL, show_pct = TRUE, pct_gp = gpar(fontsize=0), barplot_ignore = c("p.R683G","p.R683S","p.A146T"),
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))






t15=read.table("~/Desktop/SV_test/plot_1/by_project.tsv",header = T)
library(ggplot2)
library(reshape2)

m=melt(t15)

ggplot(m,aes(Project,value,fill=ANNOT))+geom_boxplot(width=0.5, alpha=0.01)+geom_jitter(size=6,position = position_dodge(width=0.5),aes(group=ANNOT, color=factor(ANNOT)), alpha=0.3)+facet_wrap( ~ variable, scales="free")+
  theme_bw(base_size = 15)+ 
  scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="right", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")+ylab("")
