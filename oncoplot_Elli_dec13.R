setwd("/Users/yellapav/Desktop/SOHN_163")
library(ComplexHeatmap)
#mat=read.table("plot1", header = TRUE, stringsAsFactors = F, sep="\t")
mat=read.table("plot", header = TRUE, stringsAsFactors = F, sep="\t")
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EDEDED", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#377EB8", col = NA))
  },
  Focal.DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "#211BC7", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E41A1C", col = NA))
  },
  Focal.AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "#A31A48", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#F9C60B", col = NA))
  },
  
  Non.Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#FF7F00", col = NA), name=c("hey"))
  },
  Synonymous = function(x, y, w, h) {
    grid.rect(x, y, w*0.25, h*0.33, gp = gpar(fill = "#4DAF4A", col = NA))
  },
  INV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#FFFF33", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.75, gp = gpar(fill = "#984EA3", col = NA))
  }
)

col = c("Synonymous" = "#4DAF4A", "AMP" = "#E41A1C", "DEL" = "#377EB8", "Non.Synonymous"="#FF7F00", "INV"="#FFFF33", "Fusion" = "#984EA3", "Focal.AMP" = "#A31A48","Focal.DEL" = "#211BC7","Splice_Site" = "#F9C60B")

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL,show_pct = TRUE, pct_gp = gpar(fontsize=8),
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Focal.AMP", "Focal.DEL","Synonymous", "Non.Synonymous","INV", "Fusion","Splice_Site"), 
                                      labels = c("Amplification", "Deletion", "Focal-Amplification", "Focal-Deletion","Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion","Splice Site")))


#oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun, col = col, 
#          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL, show_pct = TRUE, pct_gp = gpar(fontsize=0), barplot_ignore = c("p.R683G","p.R683S","p.A146T"),
#          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
#                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))
