library("biomaRt")
library("RColorBrewer")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ens_v75=getBM(attributes=c('ensembl_gene_id','external_gene_id','chromosome_name'), mart = ensembl_75)

stewd("~/Desktop/TP53_take2/")
getCts <- function(filepath) {
  counts=read.table(filepath, header=T, sep="\t",stringsAsFactors=FALSE)
  counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
  counts=counts[counts$maxVal>10, ]
  counts=counts[,-dim(counts)[2]]
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  return(counts)
}

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
  y=y[which(y$padj<0.05),]
  y$ensg=rownames(y)
  y=as.data.frame(y)
  m=merge(y, ens_v75, by.x="ensg", by.y="ensembl_gene_id", all.x= TRUE)
  write.table(m,file=outFile, append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)
  return(m)
} 




### Hotspot #####
ycts=getCts("/Users/yellapav/Desktop/TP53_take2/hotspot_nonComplex.mat")
conds=c("WT","Hotspot","WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","Hotspot","Hotspot","WT","WT","Hotspot","WT","WT","WT","WT","WT","WT")
#c("Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y")
yres=getResults(ycts, conds)
see=annotate(yres$res, ens_v75, "~/Desktop/TP53_take2/WT_vs_Hotspot_nonComplex.tsv")


#### MIssense ####
ycts=getCts("/Users/yellapav/Desktop/TP53_take2/missense_nonComplex.mat")
conds=c("WT","missense","WT","WT","missense","missense","WT","WT","WT","missense","missense","WT","WT","WT","WT","WT","missense")

yres=getResults(ycts, conds)
see=annotate(yres, ens_v75, "~/Desktop/TP53_take2/WT_vs_Missense_nonComplex.tsv")


##### FrameShift and Hostspot ######
ycts=getCts("/Users/yellapav/Desktop/record/RNA_TP53/fs.mat")
conds=c("WT","frameshift","WT","frameshift","WT","WT","WT","WT","WT","WT","frameshift","WT","WT","WT")

yres=getResults(ycts, conds)
see=annotate(yres, ens_v75, "~/Desktop/TP53_take2/WT_vs_Truncating_nonComplex.tsv")


########## All ####
ycts=getCts("/Users/yellapav/Desktop/TP53_take2/all.mat")
conds=c("Complex","Complex","Complex","Non-Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Non-Complex","Non-Complex","Non-Complex","Complex","Complex","Complex","Non-Complex","Complex","Non-Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Complex","Non-Complex","Non-Complex","Complex","Complex","Complex","Complex","Non-Complex")
kar=c("Y","Y","Y","N","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","N","N","N","Y","Y","Y","N","Y","N","Y","Y","Y","Y","Y","Y","Y","Y","N","N","Y","Y","Y","Y","N")
annotation1=c("WT","Hotspot","frameshift","Hotspot","Missense","WT","Hotspot","Hotspot","frameshift","WT","frameshift and Missense","Hotspot","Hotspot","Hotspot","Hotspot","Hotspot","Missense","Missense","Missense","Missense","WT","Missense","WT","Missense","Hotspot and Missense","WT","Missense","Missense","Missense and Nonesense","WT","WT","frameshift","Missense","Missense","WT","WT","WT","Missense","Missense")
annotation2=c("WT","Hotspot","frameshift","Hotspot","Missense","WT","Hotspot","Hotspot","frameshift","WT","Missense","Hotspot","Hotspot","Hotspot","Hotspot","Hotspot","Missense","Missense","Missense","Missense","WT","Missense","WT","Missense","Hotspot","WT","Missense","Missense","Missense","WT","WT","frameshift","Missense","Missense","WT","WT","WT","Missense","Missense")
yres=getResults(ycts, conds)

see=annotate(yres$res, ens_v75, "~/Desktop/TP53_take2/all_res.tsv")
rld <- rlog(yres$dds, blind=FALSE)
ntd <- normTransform(yres$dds)


library("pheatmap")
#cc=counts(yres$dds,normalized=TRUE)
cc=assay(rld)
zz=transform(cc,SD=apply(cc,1,sd,na.rm=TRUE))
zz_ordered=zz[order((zz$SD),decreasing=TRUE),]
zz_ordered=zz_ordered[1:250,]
zz_ordered=subset(zz_ordered,select=-SD)


df <- as.data.frame(colData(yres$dds)[,c("condition","type")])
df$annotation1=annotation1
df$annotation2=annotation2
df <- df[,c("condition","annotation1","annotation2")]

#zz_ordered=t(scale(t(zz_ordered), center = TRUE, scale = TRUE))
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
ann_colors = colorRampPalette( rev(brewer.pal(9, "Set1")) )
pheatmap((zz_ordered), cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, col=colors,scale="row",clustering_distance_rows="manhattan", clustering_distance_cols="manhattan", clustering_method="ward")


sampleDists=dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
df$Sample=paste(sub(".1.R1.1","",rownames(df)), df$condition, df$annotation1,sep="-")
colnames(sampleDistMatrix) <- df$Sample

colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=rev(colors),clustering_distance_cols=sampleDists, main="L2 distance between samples", border_color=c("#FFFFFF"))


rld$type=annotation2
plotPCA(rld, intgroup=c("condition", "type"))

pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape=condition)) +
  geom_point(size=5, alpha=0.65) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + ggtitle("PCA on Sample Distances") +
  theme_bw(base_size=15) + scale_colour_brewer(palette='Set1')


######################################

plotMA(yres, ylim=c(-3,3))
plotMA(yresLFC, ylim=c(-2,2))
plotCounts(ydds, gene=which.min(yres$padj), intgroup="condition")



rld <- rlog(yres$dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(yres$dds, blind=FALSE)
ntd <- normTransform(yres$dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))

select <- order(rowMeans(counts(yres$dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(yres$dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col = df)



sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition,colnames(rld), sep="-")
colnames(sampleDistMatrix) <- paste(rld$condition, colnames(rld), sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, col=colors, cluster_rows = TRUE, cluster_cols = TRUE)
plotPCA(rld, intgroup=c("condition", "type"))




