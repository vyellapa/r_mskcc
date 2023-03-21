setwd("/Users/yellapav/Desktop/record/EHA_2018")
plot=read.table("eha_prcessed.txt",sep="\t",header=T)
plot.t=t(plot)


library(ComplexHeatmap)

#rownames(plot)=plot$sample
#plot=plot[,!(names(plot) %in% c("sample"))]
#plot=apply(plot, 2, function(x) gsub("0","",x))

#white=plot[,colnames(plot) %in% (race[race$Race=="White",]$LeukID)]
#black=plot[,colnames(plot) %in% (race[race$Race=="Black",]$LeukID)]

colnames(plot.t)=plot.t[1,]
plot.t=plot.t[-1,]
plot.t=plot.t[,c(1:45)]
mat=as.matrix(plot.t)

hey=rownames(mat)
hey=gsub("FISH.Del..1p..","FISH.Del(1p)",hey)
hey=gsub("SEQ.Del.1p.","SEQ.Del(1p)",hey)
hey=gsub("Gain_1q.FISH","FISH.Gain(1q)",hey)
hey=gsub("Gain_1q_seq","SEQ.Gain(1q)",hey)
hey=gsub("Gain_11q.FISH.MLL","FISH.Gain(11q)",hey)
hey=gsub("Amp_11q_seq","SEQ.Gain(11q)",hey)
hey=gsub("Del_12p.FISH","FISH.Del(12p)",hey)
hey=gsub("Del_12p_seq","SEQ.Del(12p)",hey)
hey=gsub("Del_13.13q.FISH","FISH.Del(13q)",hey)
hey=gsub("Del_13.13q_seq","SEQ.Del(13q)",hey)
hey=gsub("Del_16q.FISH","FISH.Del(16q)",hey)
hey=gsub("Del_16q_seq","SEQ.Del(16q)",hey)
hey=gsub("Del_17p.FISH","FISH.Del(17p)",hey)
hey=gsub("Del_17p_seq","SEQ.Del(17p)",hey)
hey=gsub("t6.14.FISH","FISH.t(6;14)",hey)
hey=gsub("SEQ.t.6.14.","SEQ.t(6;14)",hey)
hey=gsub("t11.14.FISH","FISH.t(11;14)",hey)
hey=gsub("t11.14.seq","SEQ.t(11;14)",hey)
hey=gsub("t14.16.FISH","FISH.t(14;16)",hey)
hey=gsub("t14.16.seq","SEQ.t(14;16)",hey)
hey=gsub("t14.20.FISH","FISH.t(14;20)",hey)
hey=gsub("t14.20.seq","SEQ.t(14;20)",hey)
hey=gsub("HRD.FISH","FISH.HRD",hey)
hey=gsub("HRD_seq","SEQ.HRD",hey)
hey=gsub("t4.14.seq","SEQ.t(4;14)",hey)
rownames(mat)=gsub("t4.14.FISH","FISH.t(4;14)",hey)


alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  YES = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#1f78b4", col = NA))
  },
  NOT.AVAILABLE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#a6cee3", col = NA))
  },
  translocation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#f46d43", col = NA))
  },
  inframe_codon_loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#EA3599", col = NA))
  }
  
)

col = c("translocation" = "#f46d43", "YES" = "#1f78b4", "NOT.AVAILABLE" = "#a6cee3", "inframe_codon_loss" = "#EA3599")

rorder=rownames(mat)
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, show_pct = FALSE,
          column_title = "Correlation of SVs and CNAs by FISH and myType", show_column_names = FALSE,row_order=rorder, remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alterations Present", at = c("YES", "NOT.AVAILABLE", "translocation","inframe_codon_loss"), 
                                      labels = c("YES", "NOT PROBED", "Translocation", "Inframe Codon Loss")))


