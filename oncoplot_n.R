setwd("/Users/yellapav/Desktop/SOHN_163")

maf=read.table("~/Desktop/SOHN_163/VCF_MAF/E-H-109099-T1-1-D1-1.pindel.gm.maf",sep="\t", quote="~", header=TRUE,stringsAsFactors = FALSE)
cnv=read.table("~/Desktop/SOHN_163/SV_CNV/E-H-109099-T1-1-D1-1.cnv.tsv",sep="\t", stringsAsFactors = FALSE)
sv=read.table("~/Desktop/SOHN_163/SV_CNV/E-H-109099-T1-1-D1-1.sv.tsv", sep="\t",quote="~", stringsAsFactors = FALSE)
cen=read.table("~/local/resources/cancer_gene_census_nov7.3_51_26.2016.tsv", header=T,quote="~",sep="\t", stringsAsFactors = FALSE)

maf=maf[maf$Hugo_Symbol %in% cen$Gene.Symbol,]
sv=sv[sv$V48 %in% cen$Gene.Symbol,]
cnv=cnv[cnv$V9 %in% cen$Gene.Symbol,]

cnv$V10=abs(cnv$V6)
cnv_max=lapply(cnv$V9, function(x) {s=cnv[cnv$V9==x,]; m=max(s$V10); s=s[s$V10==m,]; return(s)})
cnv_max_df=do.call(rbind.data.frame, cnv_max)

colnames(cnv_max_df)=paste("CNV",colnames(cnv_max_df),sep="_")
merged=merge(maf,cnv_max_df,all=T,by.x="Hugo_Symbol",by.y="CNV_V9")

colnames(sv)=paste("SV",colnames(sv),sep="_")
merged=merge(merged,sv,all=T,by.x="Hugo_Symbol",by.y="SV_V48")
plot=merged[,c(1,9:10,16,37,117,134)]
plot=plot[(which((plot$CNV_V4>2.5 & plot$CNV_V4<1.5) | (!is.na(plot$SV_V12) | !is.na(plot$Variant_Type)))),]
plot$Tumor_Sample_Barcode=rep(plot$Tumor_Sample_Barcode[grep("E-",unique(plot$Tumor_Sample_Barcode))], nrow(plot))
m=melt(plot)



mat=read.table("plot", header = TRUE, stringsAsFactors = F, sep="\t")
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  Silent = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#B8DC3C", col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#A31A48", col = NA))
  },
  INV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#DB2464", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#EDDE45", col = NA))
  },
  Frame_Shift = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#EA3556", col = NA))
  },
  Flank5 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#2DD754", col = NA))
  },
  Flank3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.25, gp = gpar(fill = "#6C16B3", col = NA))
  }
)

col = c("Silent" = "#008000", "AMP" = "red", "DEL" = "blue", "Splice_Site"="#B8DC3C", "Missense" = "#A31A48", "INV"="#DB2464", "Fusion" = "#EDDE45", "Frame_Shift"="#EA3556", "Flank5"="#2DD754", "Flank3"="#6C16B3" )

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Sohn ALL", show_column_names = TRUE, column_order=NULL,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "Silent", "Splice_Site", "Missense", "INV", "Fusion", "Framse_Shift", "Flank5", "Flank3"), 
                                      labels = c("Amplification", "Deletion", "Sielnt", "Splice Site", "Missense", "Inversion", "Fusion", "Framse Shift", "5Flank", "3Flank")))


