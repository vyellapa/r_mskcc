setwd("/Users/yellapav/Desktop/SOHN_163")
library(ComplexHeatmap)
#mat=read.table("plot1", header = TRUE, stringsAsFactors = F, sep="\t")
mat=read.table("plot_99", header = TRUE, stringsAsFactors = F, sep="\t")
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

col = c("Silent" = "#4DAF4A", "AMP" = "#E41A1C", "DEL" = "#377EB8", "NS"="#FF7F00", "INV"="#FFFF33", "Fusion" = "#984EA3","NST"= "#984EA3", "p.R683S" = "#000000","p.R683G" = "#000000","p.A146T" = "#000000")

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in all ALL Samples", show_column_names = TRUE, column_order=NULL,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))


oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Aberrations in Cancer Gene Census (CGC) Genes", show_column_names = TRUE, column_order=NULL, show_pct = TRUE, pct_gp = gpar(fontsize=0), barplot_ignore = c("p.R683G","p.R683S","p.A146T"),
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "NS","INV", "Fusion"), 
                                      labels = c("Amplification", "Deletion", "Synonymous Mutation", "Non Synonymous Mutation",  "Inversion", "Fusion")))
