setwd("/Users/yellapav/Desktop/p251/SV/")
library(ComplexHeatmap)
#mat=read.table("plot1", header = TRUE, stringsAsFactors = F, sep="\t")
mat=read.table("p251_SV_oncoplot.txt", header = TRUE, stringsAsFactors = F, sep="\t")
colnames(mat)=gsub(".T1.1.D1.1","",colnames(mat))
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EDEDED", col = NA))
  },
  Translocation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#ff7f00", col = NA))
  },
  tra_SimpleRep = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#a6cee3", col = NA))
  },
  tra_LINE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#b2df8a", col = NA))
  },
  Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#1f78b4", col = NA), name=c("hey"))
  },
  
  Hyperdiploid = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#33a02c", col = NA))
  },
  tra_SINE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#fb9a99", col = NA))
  },
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#e41a1c", col = NA))
  }
)
#D53E4F
col = c("Translocation" = "#ff7f00", "tra_SimpleRep" = "#a6cee3", "tra_LINE" = "#b2df8a", "Del"="#1f78b4", "Hyperdiploid"="#33a02c", "tra_SINE" = "#fb9a99","Amp"= "#e41a1c")



#sampord=c("I-H-106917-T2-1","I-H-106917-T2-2","I-H-106917-T2-3","I-H-106917-T2-4","I-H-130718-T1-1","I-H-130718-T1-10","I-H-130718-T1-2","I-H-130718-T1-3","I-H-130718-T1-4","I-H-130718-T1-5","I-H-130718-T1-6","I-H-130718-T1-7","I-H-130718-T1-8","I-H-130718-T1-9","I-H-130719-T1-1","I-H-130719-T1-2","I-H-130719-T1-3","I-H-130719-T1-4","I-H-130719-T1-5","I-H-130719-T1-6","I-H-130720-T1-1","I-H-130720-T1-2","I-H-130720-T1-3","I-H-130720-T1-4","I-H-130720-T1-5","I-H-130720-T1-6","I-H-130720-T1-7","I-H-130720-T1-8")
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Common MM Rearrangements", show_column_names = TRUE, column_order= NULL ,
          heatmap_legend_param = list(title = "Alternations", at = c("Translocation", "tra_SimpleRep", "tra_LINE", "Del", "Hyperdiploid", "tra_SINE","Amp"), 
                                      labels = c("Translocation", "tra_SimpleRep", "tra_LINE", "Del", "Hyperdiploid", "Trans_SINE","Amp")))




oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL, show_pct = TRUE, pct_gp = gpar(fontsize=0), barplot_ignore = c("p.R683G","p.R683S","p.A146T"),
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))


