setwd("~/Desktop/gunjan/")
require(maftools)

laml = read.maf(maf = "~/Desktop/gunjan/ia12_mmrf_all_maf001.maf")

#setwd("~/Desktop/gunjan/neg")
#laml = read.maf(maf = "~/Desktop/gunjan/neg/neg.maf")

#setwd("~/Desktop/gunjan/pos")
#laml = read.maf(maf = "~/Desktop/gunjan/pos/pos.maf")

write.mafSummary(maf = laml, basename = 'all')

pdf("all_summary.pdf",width=16,height=10,paper='special')
plotmafSummary(maf = laml, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

glist=c("NRAS","KRAS","DIS3","TP53", "BRAF","TRAF3","RYR2","ATM","MAX","EGR1","FAT3","FAT4","FAM46C","LTB","FGFR3","SP140","USH2A","PRKD2","LRP1B","FAT1","ABCA13","ZFHX4","GPR98","HIST1H1E","CSMD3","CSMD1","ZNF462")
png(file = "all_oncoplot.png", res=300, width=2200, height=1800)
oncoplot(maf = laml, top = 65, fontSize = 8, genes = glist)
dev.off()

pdf("all_oncoplot.pdf",width=14,height=10,paper='special')
oncoplot(maf = laml, top = 65, fontSize = 8, genes = glist)
dev.off()


png(file = "all_oncoplot_all.png", res=300, width=3200, height=2200)
oncoplot(maf = laml, top = 85, fontSize = 6)
dev.off()

pdf("all_oncoplot_all.pdf",width=16,height=12,paper='special')
oncoplot(maf = laml, top = 85, fontSize = 8)
dev.off()



laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
png(file = "all_titv.png", res=300, width=2200, height=1600)
plotTiTv(res = laml.titv)
dev.off()

pdf("all_titv.pdf",width=12,height=8,paper='special')
plotTiTv(res = laml.titv)
dev.off()




png(file = "all_somaInt.png", res=300, width=2200, height=1600)
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1), genes = glist)
dev.off()


pdf("all_somaInt.pdf",width=12,height=8,paper='special')
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1), genes = glist)
dev.off()



aa=read.table("/Users/yellapav/Desktop/gunjan/q1_neg", header=T, sep="\t")
bb=read.table("/Users/yellapav/Desktop/gunjan/q1_pos", header=T, sep="\t")

aa$Type="1q.WT"
bb$Type="1q+"
cc=rbind(aa,bb)
colnames(cc)=c("Tumor_Sample_Barcode","log2FC","Type")

ggplot(cc, aes(Type, log2FC))+geom_boxplot(width=0.75, outlier.size = 0)+geom_point(size=4,color="#0038b8",aes(alpha=0.05))+
  theme_bw(base_size = 15)+ scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", strip.background  = element_blank(), legend.title=element_blank(),
        axis.ticks = element_line(size = 0), panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor.y = element_line(colour = "grey90"),panel.grid.major.x = element_blank())+ theme(legend.title=element_blank())+xlab("")




cna=read.table("/Users/yellapav/Desktop/gunjan/data/MMRF_CoMMpass_IA12a_CNA_LongInsert_FISH_CN_All_Specimens.txt", sep = "\t", header = TRUE)
cnab=cna[,c(1,56:163,165)]

cnab.pos=cnab[which(cnab$Study_Visit_ID %in% bb$Tumor_Sample_Barcode),]
cnab.neg=cnab[which(cnab$Study_Visit_ID %in% aa$Tumor_Sample_Barcode),]

cn.pos=t(cnab.pos)
cn.neg=t(cnab.neg)
cn.all=t(cnab)


colnames(cn.pos)=cn.pos[1,]
colnames(cn.neg)=cn.neg[1,]
colnames(cn.all)=cn.all[1,]

cn.pos=cn.pos[-1,]
cn.neg=cn.neg[-1,]
cn.all=cn.all[-1,]


library(ComplexHeatmap)
plot=cn.pos
plot.all=cn.all

#rownames(plot)=plot$sample
#plot=plot[,!(names(plot) %in% c("sample"))]
plot=apply(plot, 2, function(x) gsub("0","",x))
plot.all=apply(plot, 2, function(x) gsub("1","DEL",x))


mat=as.matrix(plot)
mat.all=as.matrix(plot.all)


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
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, gp = gpar(fill = "#2b8cbe", col = NA))
  }, 
  stop_gained = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#fcc5c0", col = NA))
  }, 
  initiator_codon_change = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#a6cee3", col = NA))
  }, 
  inframe_codon_loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.95, gp = gpar(fill = "#EA3599", col = NA))
  }
  
)

#anno=read.table("annotation.tsv")
#anno=read.table("anno.tsv", header=T)

#ha_column = HeatmapAnnotation(df = anno, col = list(type1 = c("I-H-106917" =  "#33a02c", "I-H-130718" = "#b2df8a", "I-H-130719" = "#ff7f00", "I-H-130720" = "#fdbf6f")))

col = c("translocation" = "#f46d43", "AMP" = "#d73027", "DEL" = "#4575b4", "splice_site_variant"="#B8DC3C", "non_synonymous_codon" = "#1a9850", "inframe_variant"="#a6cee3", "frameshift_variant" = "#2b8cbe", "stop_gained" = "#fcc5c0", "initiator_codon_change" = "#a6cee3","inframe_codon_loss" = "#EA3599")



oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in all Samples", show_column_names = F,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("AMP", "CNA Positive", "Translocation", "Splice Site Variant", "Missense Variant", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss")))


oncoPrint(mat.all, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Alterations in all Samples", show_column_names = F,
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
                                      labels = c("AMP", "CNA Positive", "Translocation", "Splice Site Variant", "Missense Variant", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss")))


#rorder=c("amp_1q", "del_16q", "del_13q", "amp_11q", "del_6q", "del_1p", "del_17p", "del_12p", "t(14:11)", "t(14:16)", "t(4:14)", "t(8:14)","KRAS",  "BRAF",  "FAM46C",  "NRAS","TP53", "CCND1", ,"EGR1","MAX", "CYLD", "DIS3", "FAT1", "FGFR3", "HIST1H1E","KDM6A", "ATM", "LTB", "MYC", "PTPN11","SP140" )

#oncoPrint(black, get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun, col = col, 
#          column_title = "Alterations in African American Samples", show_column_names = FALSE,row_order = rorder,
#          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
#                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Inframe Variant", "Frameshift Variant", "Stop Gain", "Initiator Codon Change","Inframe Codon Loss")))

#oncoPrint(white, get_type = function(x) strsplit(x, ";")[[1]],
#          alter_fun = alter_fun, col = col, 
#          column_title = "Alterations in CEPH European Samples", show_column_names = FALSE,row_order = rorder,
#          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "DEL", "translocation", "splice_site_variant", "non_synonymous_codon", "inframe_variant", "frameshift_variant", "stop_gained", "initiator_codon_change","inframe_codon_loss"), 
#                                      labels = c("Amplification", "Deletion", "Translocation", "Splice Site", "Missense", "Inframe Variant", "Frameshift Variant", "Stop Gain", "Initiator Codon Change","Inframe Codon Loss")))


#cn.pos.num=as.data.frame(apply(as.matrix(cn.pos), 1, as.numeric))
#cn.pos.num$row_sum=as.data.frame(apply(cn.pos.num, 1, sum))
#pos.count=dim(cn.all.num)

#cn.all.num=as.data.frame(apply(as.matrix(cn.all), 1, as.numeric))
#cn.all.num$row_sum=as.data.frame(apply(cn.all.num, 1, sum))
 #class(apply(cn.all.num, 1, sum))
#total.count=dim(cn.all.num)[1]

#########################################################################
################  CNA 1q pos results       ##############################
#########################################################################

cn.results <- data.frame(cna=character(), 
                          pos.cases=numeric(), 
                          pos.total=numeric(),
                          neg.cases=numeric(),
                          neg.total=numeric(),
                          pval=numeric(),
                          odds.ratio=numeric(),
                          stringsAsFactors=FALSE)


cn.pos50=cn.pos[rownames(cn.pos) %in% grep("20percent",rownames(cn.pos),value=T,invert = T),]
cn.neg50=cn.neg[rownames(cn.neg) %in% grep("20percent",rownames(cn.neg),value=T,invert = T),]
cn.all50=cn.all[rownames(cn.all) %in% grep("20percent",rownames(cn.all),value=T,invert = T),]

for(i in 1:dim(cn.pos50)[1]) {
  #for(i in 1:2) {
  a=sum(as.numeric(cn.pos50[i,]))
  b=dim(cn.pos50)[2]
  c=sum(as.numeric(cn.neg50[i,]))
  d=dim(cn.neg50)[2]
  e=rownames(cn.pos50)[i]
  f.df=matrix(c(a,b,c,d),nrow = 2)
  ftest=fisher.test(f.df)
  j=paste(a,b,c,d,e,sep=",")
  print(j)
  cn.results <- rbind(cn.results,data.frame(cna=e,pos.cases=a,pos.total=b,neg.cases=c,neg.total=d,pval=ftest$p.value,odds.ratio=as.numeric(ftest$estimate)))
}
cn.results$pos.freq = (cn.results$pos.cases/cn.results$pos.total)
cn.results$neg.freq = (cn.results$neg.cases/cn.results$neg.total)
cn.results.sub=cn.results[cn.results$pos.freq>0.1,]
cn.results.sub$qval=p.adjust(cn.results.sub$pval, method = "fdr")
cn.signif=cn.results.sub[cn.results.sub$qval<0.1,]
cn.signif$cna=gsub("SeqWGS_Cp_","", cn.signif$cna)
cn.signif$cna=gsub("_50percent","", cn.signif$cna)

#write.table(cn.signif,file="/Users/yellapav/Desktop/gunjan/cna.signif.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)
write.table(cn.signif,file="/Users/yellapav/Desktop/gunjan/cna.pos.neg.signif.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)


#########################################################################
################   SV 1q pos results       ##############################
#########################################################################

sv=read.table("/Users/yellapav/Desktop/gunjan/data/MMRF_CoMMpass_IA12a_LongInsert_Canonical_Ig_Translocations.txt", sep = "\t", header = TRUE)
cal=as.vector(grep("CALL",colnames(sv)))
samp=c("1")
cal=c(samp,cal)
sv.call=sv[,as.numeric(cal)]
sv.call.pos=sv.call[sv.call$Study_Visit_iD %in% bb$Tumor_Sample_Barcode,]
sv.call.neg=sv.call[sv.call$Study_Visit_iD %in% aa$Tumor_Sample_Barcode,]
sv.call=t(sv.call)
sv.call.pos=t(sv.call.pos)
sv.call.neg=t(sv.call.neg)
colnames(sv.call)=sv.call[1,]
sv.call=sv.call[-1,]

colnames(sv.call.pos)=sv.call.pos[1,]
colnames(sv.call.neg)=sv.call.neg[1,]
sv.call.pos=sv.call.pos[-1,]
sv.call.neg=sv.call.neg[-1,]


sv.results <- data.frame(sv=character(), 
                         pos.cases=numeric(), 
                         pos.total=numeric(),
                         neg.cases=numeric(),
                         neg.total=numeric(),
                         pval=numeric(),
                         odds.ratio=numeric(),
                         stringsAsFactors=FALSE)


for(i in 1:dim(sv.call.pos)[1]) {
  #for(i in 1:2) {
  a=sum(as.numeric(sv.call.pos[i,]))
  b=dim(sv.call.pos)[2]
  c=sum(as.numeric(sv.call.neg[i,]))
  d=dim(sv.call.neg)[2]
  e=rownames(sv.call.pos)[i]
  f.df=matrix(c(a,b,c,d),nrow = 2)
  ftest=fisher.test(f.df)
  j=paste(a,b,c,d,e,sep=",")
  print(j)
  sv.results <- rbind(sv.results, data.frame(sv=e,pos.cases=a,pos.total=b,neg.cases=c,neg.total=d,pval=ftest$p.value,odds.ratio=as.numeric(ftest$estimate)))
}

sv.results$freq.pos = (sv.results$pos.cases/sv.results$pos.total)
sv.results$freq.neg = (sv.results$neg.cases/sv.results$neg.total)

sv.results$qval=p.adjust(sv.results$pval, method = "fdr")
#sv.signif=sv.results[sv.results$qval<0.1,]
sv.results$sv=gsub("SeqWGS_","", sv.results$sv)
#sv.results$cna=gsub("_50percent","", cn.signif$cna)

#write.table(sv.results,file="/Users/yellapav/Desktop/gunjan/sv.results.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)
write.table(sv.results,file="/Users/yellapav/Desktop/gunjan/sv.pos.neg.results.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE)

#########################################################################
################  CNA Somatic Interactions ##############################
#########################################################################

cn.pos.somint=melt(cn.pos)
cn.neg.somint=melt(cn.neg)
cn.all.somint=melt(cn.all)

write.table(cn.pos.somint,file="/Users/yellapav/Desktop/gunjan/data/cn.pos.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
write.table(cn.neg.somint,file="/Users/yellapav/Desktop/gunjan/data/cn.neg.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
write.table(cn.all.somint,file="/Users/yellapav/Desktop/gunjan/data/cn.all.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)

system("/Users/yellapav/Desktop/gunjan/data/cna2maf.sh",intern=TRUE)

cna.pos = read.maf(maf = "/Users/yellapav/Desktop/gunjan/data/cna.pos.maf")
cna.neg = read.maf(maf = "/Users/yellapav/Desktop/gunjan/data/cna.neg.maf")
cna.all = read.maf(maf = "/Users/yellapav/Desktop/gunjan/data/cna.all.maf")

png(file = "/Users/yellapav/Desktop/gunjan/pos.cna.somaInt.png", res=300, width=2200, height=1600)
somaticInteractions(maf = cna.pos, top = 35, pvalue = c(0.01, 0.5))
dev.off()
pdf("/Users/yellapav/Desktop/gunjan/pos.cna.somaInt.pdf",width=12,height=8,paper='special')
somaticInteractions(maf = cna.pos, top = 35, pvalue = c(0.01, 0.5))
dev.off()

png(file = "/Users/yellapav/Desktop/gunjan/neg.cna.somaInt.png", res=300, width=2200, height=1600)
somaticInteractions(maf = cna.neg, top = 35, pvalue = c(0.01, 0.5))
dev.off()
pdf("/Users/yellapav/Desktop/gunjan/neg.cna.somaInt.pdf",width=12,height=8,paper='special')
somaticInteractions(maf = cna.neg, top = 35, pvalue = c(0.01, 0.5))
dev.off()

png(file = "/Users/yellapav/Desktop/gunjan/all.cna.somaInt.png", res=300, width=2200, height=1600)
somaticInteractions(maf = cna.all, top = 35, pvalue = c(0.01, 0.05))
dev.off()
pdf("/Users/yellapav/Desktop/gunjan/all.cna.somaInt.pdf",width=12,height=8,paper='special')
somaticInteractions(maf = cna.all, top = 55, pvalue = c(0.01, 0.05))
dev.off()


#########################################################################
################  SV Somatic Interactions ##############################
#########################################################################


#cna=e,pos.cases=a,pos.total=b,all.cases=c,all.total=d,pval=ftest$p.value,odds.ratio=ftest$estimate

#########################################################################
################       RNA-Seq DESeq2           #########################
#########################################################################


#For Gene annotation
library("biomaRt")
library("RColorBrewer")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ens_v75=getBM(attributes=c('ensembl_gene_id','external_gene_id','chromosome_name'), mart = ensembl_75)

#Read counts and remove ENSGs if max value is < 10
stewd("/Users/yellapav/Desktop/gunjan/data")
getCts <- function(filepath) {
  counts=read.table(filepath, header=T, sep="\t",stringsAsFactors=FALSE)
  counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
  counts=counts[counts$maxVal>25, ]
  counts=counts[,-dim(counts)[2]]
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  return(counts)
}


#DESEQ2
getResults <- function(ycts, conds) {
  ycoldata=as.data.frame(cbind(colnames(ycts),conds,rep("paired",length(conds))))
  rownames(ycoldata)=ycoldata[,1]
  ycoldata=ycoldata[,-1]
  colnames(ycoldata)=c("condition", "type")
  
  library("DESeq2")
  
  ydds <- DESeqDataSetFromMatrix(countData = ycts, colData = ycoldata, design = ~ condition)
  ydds <- ydds[ rowSums(counts(ydds)) > 20, ]
  
  ydds <- DESeq(ydds)
  yres <- results(ydds)
  l1=list(dds=ydds, res=yres)
  return(l1)
}

annotate <- function(yres, ens_v75, outFile) {
  y=as.data.frame(yres$res)
  y=y[which(y$padj<0.1),]
  y$ensg=rownames(y)
  y=as.data.frame(y)
  m=merge(y, ens_v75, by.x="ensg", by.y="ensembl_gene_id", all.x= TRUE)
  write.table(m,file=outFile, append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)
  return(m)
} 


# read counts and seperate 1q+ and WT
ycts=getCts("/Users/yellapav/Desktop/gunjan/data/MMRF_CoMMpass_IA12a_E74GTF_HtSeq_Gene_Counts.txt.gz")
ycts.neg=ycts[,(colnames(ycts) %in% aa$Tumor_Sample_Barcode)]
ycts.pos=ycts[,(colnames(ycts) %in% bb$Tumor_Sample_Barcode)]
sv.call.melt = melt(sv.call)

# get results for whsc1+ cohort of 1q+ and WT
for(i in grep("CCND2", unique(sv.call.melt$Var1), value = T, invert = T))
{
  print(i)
  sv.call.melt.sub = sv.call.melt[which(sv.call.melt$Var1==i & sv.call.melt$value=="1"),]
  ycts.neg.sub=ycts.neg[,(colnames(ycts.neg) %in% sv.call.melt.sub$Var2)]
  ycts.pos.sub=ycts.pos[,(colnames(ycts.pos) %in% sv.call.melt.sub$Var2)]
  
  stat=c("NEG/POS",dim(ycts.neg.sub)[2],"/",dim(ycts.pos.sub)[2])
  
  ycts.sub=cbind(ycts.neg.sub, ycts.pos.sub)
  conds=c(rep("NEG",dim(ycts.neg.sub)[2]),rep("POS",dim(ycts.pos.sub)[2]))
  #conds=c("WT","Hotspot","WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","Hotspot","Hotspot","WT","WT","Hotspot","WT","WT","WT","WT","WT","WT")
  yres=getResults(ycts.sub, conds)
  outf=paste("/Users/yellapav/Desktop/gunjan/results/", i,"_deseq2.tsv",sep = "")
  statf=paste("/Users/yellapav/Desktop/gunjan/results/", i,"_stats.tsv",sep = "")
  all.res=annotate(yres, ens_v75, outf)
  write.table(stat,file=statf, append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)
  
  
  #rld <- rlog(yres$dds, blind=FALSE)
  #ntd <- normTransform(yres$dds)
  
  #####rm later #######
  #see=as.data.frame(yres$res)
  #write.table(see[see$padj<0.1,],file=outf, append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)
  ######
}

i="SeqWGS_MYC_CALL"
print(i)
sv.call.melt.sub = sv.call.melt[which(sv.call.melt$Var1==i & sv.call.melt$value=="1"),]
ycts.neg.sub=ycts.neg[,(colnames(ycts.neg) %in% sv.call.melt.sub$Var2)]
ycts.pos.sub=ycts.pos[,(colnames(ycts.pos) %in% sv.call.melt.sub$Var2)]
ycts.sub=cbind(ycts.neg.sub, ycts.pos.sub)
print(c(dim(ycts.neg.sub)[2],dim(ycts.pos.sub)[2]))
conds=c(rep("NEG",dim(ycts.neg.sub)[2]),rep("POS",dim(ycts.pos.sub)[2]))
#conds=c("WT","Hotspot","WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","Hotspot","Hotspot","WT","WT","Hotspot","WT","WT","WT","WT","WT","WT")
yres=getResults(ycts.sub, conds)
outf=paste("/Users/yellapav/Desktop/gunjan/results/", i,"_deseq2.tsv",sep = "")
##all.res=annotate(yres$res, ens_v75, outf)
#####rm later #######
see=as.data.frame(yres$res)
write.table(see[see$padj<0.1,],file=outf, append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)


rld <- rlog(yres$dds, blind=FALSE)
library("pheatmap")
nmlzd = counts(yres$dds,normalized=TRUE)
cc=assay(rld)
zz=transform(cc,SD=apply(cc,1,sd,na.rm=TRUE))
zz_ordered=zz[order((zz$SD),decreasing=TRUE),]
zz_ordered=zz_ordered[1:250,]
zz_ordered=subset(zz_ordered,select=-SD)


df <- as.data.frame(colData(yres$dds)[,c("condition","type")])
df$annotation1=conds
# df$annotation2=annotation2
df <- df[,c("condition","annotation1")]

#zz_ordered=t(scale(t(zz_ordered), center = TRUE, scale = TRUE))
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
ann_colors = colorRampPalette( rev(brewer.pal(9, "Set1")) )
pheatmap((zz_ordered), cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, col=colors,scale="row",clustering_distance_rows="manhattan",
         clustering_distance_cols="manhattan", clustering_method="ward")


sampleDists=dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
df$Sample=paste(sub(".1.R1.1","",rownames(df)), df$condition, df$annotation1,sep="-")
colnames(sampleDistMatrix) <- df$Sample

colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=rev(colors),clustering_distance_cols=sampleDists, main="L2 distance between samples", border_color=c("#FFFFFF"))


#Calcualte pairwise distances and plot
sampleDists=dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
df$Sample=paste(sub(".1.R1.1","",rownames(df)), df$condition, df$annotation1,sep="-")
colnames(sampleDistMatrix) <- df$Sample

colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=rev(colors),clustering_distance_cols=sampleDists, main="L2 distance between samples", border_color=c("#FFFFFF"))


rld$type=annotation1
plotPCA(rld, intgroup=c("condition", "type"))

pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape=condition)) +
  geom_point(size=5, alpha=0.65) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + ggtitle("PCA on Sample Distances") +
  theme_bw(base_size=15) + scale_colour_brewer(palette='Set1')





################################################################
######################## Gowers Clustering #####################
###############################################################

maf.all<-read.table("~/Desktop/gunjan/ia12_mmrf_all_maf001.maf", sep="\t", header=T, quote = "~")
#maf.all = maf.all[,c(1:32)]
NS=c("Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","Splice_Site","Frame_Shift_Ins","In_Frame_Del","Translation_Start_Site","Nonstop_Mutation","In_Frame_Ins")
glist=c("DIS3","LRP1B","ATM","TRAF3","USH2A","NRAS")
cnlist=c("13q14","13q34","16q22","CDKN2C","FAM46C","RB1","11p15","11q23","17q23","18p11","18q21","1q21","2p23","6p22")
#svlist=c("A")

rna.list=c()
maf.all = maf.all[which(maf.all$Variant_Classification %in% NS),]


see = maf.all[,c("Hugo_Symbol","Variant_Type","Tumor_Sample_Barcode")]
s = dcast(see,Hugo_Symbol~Tumor_Sample_Barcode )
mut.matrix = s[which(s$Hugo_Symbol %in% glist),]


#Getting normalised counts ###########################################
ycts=getCts("/Users/yellapav/Desktop/gunjan/data/MMRF_CoMMpass_IA12a_E74GTF_HtSeq_Gene_Counts.txt.gz")
conds=c(rep("NEG",440),rep("POS",441))
yres=getResults(ycts, conds)
nmlzd = counts(yres$dds,normalized=TRUE)
write.table(nmlzd,file="/Users/yellapav/Desktop/gunjan/results/normalized_counts.txt", append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)
rna.sig=read.table("/Users/yellapav/Desktop/gunjan/data/gower_rna.in", sep = "\t", header = FALSE)
rna.matrix=nmlzd[which(rownames(nmlzd) %in% rna.sig$V2),]


# CNA

cn.all.20p=cn.all[grep("20percent",rownames(cn.all)),]
rownames(cn.all.20p)=gsub("SeqWGS_Cp_","",rownames(cn.all.20p))
rownames(cn.all.20p)=gsub("_20percent","",rownames(cn.all.20p))
#maf.all = maf.all[which(maf.all$Hugo_Symbol %in% glist),]

############################################################################




######################### Gowers Method ###################################
library(gplots)
library(cluster)
library(ggplot2)
library(reshape2)
library(ggdendro)

setwd("~/Desktop/icluster_v2")
system("./s1")
sample1<-read.table("processed.tsv", sep="\t", header=T,row.names=1)
#tempData=log(sample1+0.1)
tempData=sample1
pdf("gower.pdf", width=20,height=15)
#panhm=heatmap.2(as.matrix(tempData), dendrogram="both",Rowv=TRUE,Colv=TRUE, distfun=function(tempData) daisy(tempData,metric = "gower"), hclustfun= function(tempData) hclust(tempData,method = 'ward'),colsep=c(1), ColSideColors= cd, col = bluered, trace="none",scale=c("row"))
panhm=heatmap.2(as.matrix(tempData), dendrogram="both",Rowv=TRUE,Colv=TRUE, distfun=function(tempData) daisy(tempData,metric = "gower"), hclustfun= function(tempData) hclust(tempData,method = 'ward'),colsep=c(1), col = bluered, trace="none",scale=c("row"))
row=panhm$rowInd
col=panhm$colInd

dev.off()

m=as.matrix(tempData)
m=m[row,col]
write.table(m,file="m", append=FALSE, quote=F,sep="\t", eol="\n", row.names=TRUE, col.names=TRUE)


