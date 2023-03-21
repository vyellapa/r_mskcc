## Bialleleic Part II after manual curation
library(reshape2)
plot=read.table("/Users/yellapav/Desktop/p292/scripts/data/fig4.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
plot=read.table("/Users/yellapav/Desktop/p292/scripts/data/plot.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
genes=c("TP53","FAM46C","RB1","TRAF3","CDKN2C","BIRC2","CYLD","PTEN","TNFAIP3","WWOX")
genes = c('TP53','RB1','DIS3','FAM46C','CYLD','TRAF3') 

plot = plot %>% dplyr::filter(gene %in% genes)


head(plot)

library(ComplexHeatmap)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#4575b4", col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#984ea3", col = NA))
  },
  BOTH = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#e41a1c", col = NA))
  }
  
)


col = c( "loss" = "#4575b4", "SNV"="#984ea3", "BOTH" = "#e41a1c" )

plot.c=dcast(plot,gene~ID,paste, collapse = ";")
plot.c= data.frame(lapply(plot.c, function(x) { gsub("Indel|SNV|SNV;Indel", "SNV", x) }))
plot.c= data.frame(lapply(plot.c, function(x) { gsub("loss;SNV", "BOTH", x) }))

#plot.c=plot.c[match(ord, plot.c$arm),]
mat=plot.c
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

#fi=paste("/ifs/res/leukgen/projects/292/RESULTS/misc/fish/figs/",as.character(i),"_d15.png",sep="")
png(file = "/Users/yellapav/Desktop/p292/scripts/fig4.png", res=300, width=5000, height=400)
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, show_pct = FALSE, right_annotation = NULL,top_annotation = NULL,
          column_title = "", show_column_names = FALSE, remove_empty_columns=FALSE, row_order=genes, column_names_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 14),
          heatmap_legend_param = list(title = "", at = c("BOTH","SNV", "loss" ), 
                                      labels = c("Deletion + Mutation", "Mutation","Deletion")))
dev.off()


