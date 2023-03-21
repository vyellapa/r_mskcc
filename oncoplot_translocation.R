setwd("/Users/yellapav/Desktop/p220/sciclone")
for(i in list.files(pattern="I-H-130720-T1-7-D1-1.scin")) {
  
  print(i)
  j=strsplit(i, "\\.")
  seg=paste(j[[1]][1],".seg", sep="")
  fig=paste(j[[1]][1],".png", sep="")
v1=read.table(i,header=F)
colnames(v1)=c("chr", "pos", "ref_reads", "var_reads", "vaf")

cn1 = read.table(seg)



sc = sciClone(vafs=v1, copyNumberCalls=cn1, sampleNames=j[[1]][1])

writeClusterTable(sc, paste(j[[1]][1], "out", sep="."))
#png(file = fig, res=300, width=1200, height=1200)
sc.plot1d(sc, fig)


}

dens=read.table("/Users/yellapav/Desktop/p220/sciclone/density.plot", sep = "\t")

library(ggplot2)

ggplot(dens, aes(V5))+geom_density()+facet_wrap(~V6)+coord_cartesian(ylim=c(0,0.2))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+xlab("VAFs of autosomes")

ggplot(dens, aes(V5))+geom_histogram()+facet_wrap(~V6)+coord_cartesian(ylim=c(0,40))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+xlab("VAFs of autosomes")

ggplot(dens, aes(V6,V5))+geom_boxplot(outlier.size=0)+geom_jitter(size=2,color="blue", alpha=0.3, position = position_jitter(width = .1))+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top", strip.background  = element_blank(), legend.title=element_blank(),
        panel.border = element_blank(), axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank())+ theme(legend.title=element_blank())+ylab("VAFs at Autosomes")+xlab("")





setwd("/Users/yellapav/Desktop/p220/")
library(ComplexHeatmap)
#mat=read.table("plot1", header = TRUE, stringsAsFactors = F, sep="\t")
mat=read.table("~/Desktop/p220/igh_rearr_oncoplot.txt", header = TRUE, stringsAsFactors = F, sep="\t")
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


