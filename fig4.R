suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(reshape2)))

df=read.table("/Users/yellapav/Desktop/p292/scripts/biallelic.txt",sep="\t",header=T)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  loss = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#377eb8", col = NA))
  },
  gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#d73027", col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.95, gp = gpar(fill = "#4daf4a", col = NA))
  },
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w*0.7, h*0.8, gp = gpar(fill = "#6a3d9a", col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#1a9850", col = NA))
  }, 
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.95, h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#2b8cbe", col = NA))
  }, 
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
    sequenced = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#A9A9A9", col = NA))
  }, 
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#EA3599", col = NA))
  }
  
)

col = c("SNV" = "#4daf4a", "gain" = "#d73027", "loss" = "#377eb8", "Indel"="#6a3d9a", "Missense_Mutation" = "#1a9850", "Frame_Shift_Del"="#a6cee3", "In_Frame_Del" = "#EA3599", "Nonsense_Mutation" = "#a6cee3", "Frame_Shift_Ins" = "#2b8cbe", "sequenced" = "#A9A9A9" )
row_order = c('13q','RB1','DIS3','17p','TP53','1p','FAM46C','6p','LTB','SP140')
row_order = c('13q','RB1','DIS3','17p','TP53','1p','FAM46C')
plot.c=dcast(df,arm~ID,paste, collapse = ";")
#plot.c=plot.c[match(ord, plot.c$arm),]

mat = plot.c
rownames(mat) = mat[,1]
mat = mat[,-1]
mat = as.matrix(mat)
mat = mat[row_order,]


#fi=paste("/ifs/res/leukgen/projects/292/RESULTS/misc/fish/figs/",as.character(i),"_d15.png",sep="")
png(file = "~/see.png", res=300, width=2800, height=600)
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "", show_column_names = FALSE, 
	remove_empty_columns=FALSE, 
	row_order=row_order, show_pct = FALSE, 
	column_names_gp = gpar(fontsize = 5),
	row_names_gp = gpar(fontsize = 8),
	#row_split = c(rep("13q",3),rep("17q",2), rep("1p",2),rep("6p",3)),
	row_split = c(rep("13q",3),rep("17q",2), rep("1p",2)), row_km=2,
          heatmap_legend_param = list(title = "Alternations", at = c("gain", "loss", "SNV", "Indel", "Missense_Mutation", "Frame_Shift_Del", "In_Frame_Del", "Nonsense_Mutation", "Frame_Shift_Ins", "sequenced"), 
                                      labels = c("Amplification", "Deletion", "SNV", "Indel", "Missense", "Frame Shift Del", "In-Frame Del", "Nonsense_Mutation", "Frame_Shift_Ins", "Sequenced")))

dev.off()



