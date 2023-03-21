#Maroon #752F3C
#steel blue #abb9c9
new=read.table("/Users/yellapav/Desktop/p292/scripts/new.txt", header=T,sep="\t")
library(ComplexHeatmap)

colors=c("#4CAE55","#4C4537","#E7AB1B")
colors=c("#06A77D","#035d46","#F1A208","#3E37A1","#005377","#D5C67A")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  BOTH = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#035d46", col = NA))
  },
  Pathology = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#752F3C", col = NA))
  },
  myTYPE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#F1A208", col = NA))
  },
   Sana = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#abb9c9", col = NA))
  },
  Sequenced = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#ffffff", col = NA))
  }
  
)


col = c("myTYPE" = "#F1A208", "Pathology" = "#752F3C", "BOTH" = "#035d46", "Sana" = "#abb9c9", "Sequenced" = "#ffffff" )


order=c('1p','1q','6q','8p','13q','14q','16q','17p','HRD','t(4;14)','t(6;14)','t(8;14)','t(11;14)','t(14;16)','t(14;20)')



n=new
colnames(n)=gsub(".T[1,2,3].1.D1.1","",colnames(n))
colnames(n)=gsub(".T1.2.D1.1",".2",colnames(n))
colnames(n)=gsub("I.H.","H",colnames(n))

mat=n
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

mat=gsub("myTYPE;BOTH","BOTH",mat)
mat=gsub("Pathology;BOTH","BOTH",mat)
mat=gsub("BOTH;Pathology","BOTH",mat)
mat=gsub("Pathology;myTYPE","BOTH",mat)
mat=gsub("myTYPE;Pathology","BOTH",mat)
mat=gsub("BOTH;Pathology","BOTH",mat)

fi=paste("/Users/yellapav/Desktop/p292/scripts/",as.character("fig2_"),"sample.png",sep="")
png(file = fi, res=300, width=2500, height=4000)
oncoPrint(t(mat), get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col, right_annotation = rowAnnotation(width = unit(5, "cm"),rbar = anno_oncoprint_barplot(bar_width=0.8)),
          column_title = "", show_row_names = FALSE,show_column_names = TRUE, remove_empty_columns=FALSE, column_names_gp = gpar(fontsize = 12),row_names_gp = gpar(fontsize = 7), column_order=order, show_pct=F, column_split = c(rep(c("CNA"), 9), rep(c("Translocations"), 6)), border = TRUE,
          heatmap_legend_param = list(title = "Assay", at = c("BOTH", "Pathology", "myTYPE", "Sana"), 
                                      labels = c("myTYPE & Array/FISH", "Array/FISH", "myTYPE", "NA")))

#draw(ht_list, col_split = c(rep(c("C"), 12), rep(c("D"), 3)))

dev.off()

