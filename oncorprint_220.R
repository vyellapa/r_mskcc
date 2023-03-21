setwd("/Users/yellapav/Desktop/p220/post_process/mafs")
library(ComplexHeatmap)


mat=read.table("onco_all.plot", header = TRUE, stringsAsFactors = F, sep="\t")
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
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#d73027", col = NA))
  },
  Translocation = function(x, y, w, h) {
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

anno=read.table("annotation.tsv")
anno=read.table("anno.tsv", header=T)

ha_column = HeatmapAnnotation(df = anno, col = list(type1 = c("I-H-106917" =  "#33a02c", "I-H-130718" = "#b2df8a", "I-H-130719" = "#ff7f00", "I-H-130720" = "#fdbf6f")))

col = c("Translocation" = "#f46d43", "Amp" = "#d73027", "Del" = "#4575b4", "Splice_Site"="#B8DC3C", "Missense_Mutation" = "#1a9850", "Frame_Shift_Del"="#a6cee3" )

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in all Samples", show_column_names = TRUE, column_order=NULL,row_order=NULL,bottom_annotation = ha_column,
          heatmap_legend_param = list(title = "Alternations", at = c("Amp", "Del", "Translocation", "Splice_Site", "Missense_Mutation", "Frame_Shift_Del"), 
                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Frame Shift Del")))


