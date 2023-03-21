library("reshape2")
library(hash)
setwd("~/Desktop/p212/oncoprint/") 
cnv=read.table("/Users/yellapav/Desktop/p212/oncoprint/p212_cnv_all.tsv",sep="\t", header=TRUE, stringsAsFactors = F)

colnames(cnv) = c("sample","Del_8p","Del_6q","Del_1p","Del_17p","Del_16q","Del_13q","Del_12p","Amp_1q","Amp_11q")
cnv[cnv$Amp_1q==1,]$Amp_1q="AMP"
cnv[cnv$Amp_11q==1,]$Amp_11q="AMP"
temp=as.data.frame(apply(as.matrix(cnv), 2, function(x) (gsub("1", "DEL", x))))
temp$sample=cnv$sample
cnv=as.data.frame(t(temp))

cnv$sample=rownames(cnv)
rownames(cnv)=c()

colnames(cnv)=as.character(as.matrix(cnv[1,]))
cnv=cnv[cnv$sample!="Del_8p",]
cnv=cnv[-1,]

trans=read.table("/Users/yellapav/Desktop/p212/oncoprint/p212_trans_all.tsv",sep="\t", header=TRUE, stringsAsFactors = F)
trans$dict=paste(trans$t.CHR1.CHR2.,trans$SAMPLE,sep="__")
trans$value="translocation"

h=hash()
apply(trans, 1, function(x) {(h[[x[[3]]]]=x[[4]])})
df <- data.frame(sample=character(), 
                 translocation=character(),
                 value=character(), 
                 stringsAsFactors=FALSE) 

for(i in unique(trans$t.CHR1.CHR2.)) {
  for(j in grep("^E",colnames(cnv),value=T)) {
    #print(c(i,j))
    key=as.character(paste(i,j,sep="__"))
    #print(key)
    if(is.null(h[[key]])) {
      df=rbind(df,data.frame(sample=as.character(j),translocation=as.character(i),value=""))
    }
    else {
      df=rbind(df,data.frame(sample=as.character(j),translocation=as.character(i),value="translocation"))
    }
    
  }
}

aa=as.data.frame(t(as.matrix(dcast(df,formula = sample ~ translocation))))
aa$sample=rownames(aa)
rownames(aa)=c()
colnames(aa)=as.matrix(aa)[1,]
aa=aa[-1,]
#colnames(aa)=aa[1,]

snv <- data.frame(sample=character(), 
                 gene=character(),
                 ann=character(), 
                 stringsAsFactors=FALSE) 



caveman=read.table("/Users/yellapav/Desktop/p212/p212_caveman_consensus.txt",sep="\t", header=TRUE, stringsAsFactors = F)

#caveman=caveman[which(caveman$MANUAL_ANNOTATION %in% c("LIKELY","LIKELY ","ONCOGENIC","UNKNOWN")),]
caveman=caveman[which(caveman$MANUAL_ANNOTATION %in% c("LIKELY","LIKELY ","ONCOGENIC")),]
caveman=caveman[which(caveman$TARGET_NAME %in% colnames(cnv)),]
caveman=caveman[,c("GENE","EFFECT","TARGET_NAME")]

pindel=read.table("/Users/yellapav/Desktop/p212/p212_pindel_consensus.txt",sep="\t", header=TRUE, stringsAsFactors = F)
#pindel=pindel[which(pindel$MANUAL_ANNOTATION %in% c("LIKELY","ONCOGENIC","UNKNOWN")),]
pindel=pindel[which(pindel$MANUAL_ANNOTATION %in% c("LIKELY","ONCOGENIC")),]
pindel=pindel[which(pindel$TARGET_NAME %in% colnames(cnv)),]
pindel=pindel[,c("GENE","EFFECT","TARGET_NAME")]


caveman=rbind(caveman,pindel)


for(i in grep("^E",colnames(cnv),value=T)) {
  for(j in unique(caveman$GENE)) {
    
    annot=paste(unique(caveman[which(caveman$TARGET_NAME==i & caveman$GENE==j),]$EFFECT), collapse=";")
    snv=rbind(snv,data.frame(sample=as.character(i),gene=as.character(j),ann=as.character(annot)))
    
    
  }
}

snv$sample=as.character(snv$sample)
snv$gene=as.character(snv$gene)
snv$ann=as.character(snv$ann)

snv=dcast(snv,formula = gene ~ sample, value.var = "ann")
snv$sample=snv$gene
snv=snv[,-1]


plot=rbind(cnv,snv)
plot=rbind(plot,aa)

race=read.table("/Users/yellapav/Desktop/p212/p212_race.txt",sep="\t", header=TRUE, stringsAsFactors = F)


library(ComplexHeatmap)

rownames(plot)=plot$sample
plot=plot[,!(names(plot) %in% c("sample"))]
plot=apply(plot, 2, function(x) gsub("0","",x))

white=plot[,colnames(plot) %in% (race[race$Race=="White",]$LeukID)]
black=plot[,colnames(plot) %in% (race[race$Race=="Black",]$LeukID)]


mat=as.matrix(plot)
white=as.matrix(white)
black=as.matrix(black)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#e0e0e0", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#4575b4", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#d73027", col = NA))
  },
  translocation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#f46d43", col = NA))
  },
  splice_site_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#B8DC3C", col = NA))
  },
  non_synonymous_codon = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#1a9850", col = NA))
  }, 
  inframe_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  frameshift_variant = function(x, y, w, h) {
    grid.rect(x, y, w*0.6, h*0.6, gp = gpar(fill = "#F8CA00", col = NA))
  }, 
  stop_gained = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#fcc5c0", col = NA))
  }, 
  initiator_codon_change = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  inframe_codon_loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#490A3D", col = NA))
  }
  
)

#anno=read.table("annotation.tsv")
#anno=read.table("anno.tsv", header=T)

#ha_column = HeatmapAnnotation(df = anno, col = list(type1 = c("I-H-106917" =  "#33a02c", "I-H-130718" = "#b2df8a", "I-H-130719" = "#ff7f00", "I-H-130720" = "#fdbf6f")))

col = c("translocation" = "#f46d43", "AMP" = "#d73027", "DEL" = "#4575b4", "splice_site_variant"="#B8DC3C", "non_synonymous_codon" = "#1a9850", "inframe_variant"="#a6cee3", "frameshift_variant" = "#F8CA00", "stop_gained" = "#fcc5c0", "initiator_codon_change" = "#a6cee3","inframe_codon_loss" = "#490A3D")

#oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun, col = col, 
#          column_title = "Alterations in all Samples", show_column_names = TRUE, column_order=NULL,row_order=NULL,bottom_annotation = ha_column,
#          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
#                                      labels = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss")))



oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in all Samples", show_column_names = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss")))


rorder=c("Amp_1q", "Del_16q", "Del_13q", "Amp_11q", "Del_6q", "Del_1p", "Del_17p", "Del_12p", "t(11:14)", "t(14:16)", "t(4:14)", "t(8:14)","KRAS",  "BRAF",  "FAM46C",  "NRAS","TP53", "CCND1","EGR1","MAX", "CYLD", "DIS3", "FAT1", "FGFR3", "HIST1H1E","KDM6A", "ATM", "LTB", "MYC", "PTPN11","SP140" )

oncoPrint(black, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in African American Samples", show_column_names = FALSE,row_order = rorder,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Inframe Variant", "Frameshift Variant", "Stop Gain", "Initiator Codon Change","Inframe Codon Loss")))

oncoPrint(white, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in CEPH European Samples", show_column_names = FALSE,row_order = rorder,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Inframe Variant", "Frameshift Variant", "Stop Gain", "Initiator Codon Change","Inframe Codon Loss")))



oncoPrint(black, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in African American Samples", show_column_names = TRUE,row_order = rorder,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Inframe Variant", "Frameshift Variant", "Stop Gain", "Initiator Codon Change","Inframe Codon Loss")))

oncoPrint(white, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in CEPH European Samples", show_column_names = TRUE,row_order = rorder,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Inframe Variant", "Frameshift Variant", "Stop Gain", "Initiator Codon Change","Inframe Codon Loss")))



