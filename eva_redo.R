library("biomaRt")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="mmusculus_gene_ensembl")
ens_v75=getBM(attributes=c('ensembl_gene_id','external_gene_id','chromosome_name'), mart = ensembl_75)

setwd("/Users/yellapav/Desktop/record/RNA_TP53/matahi")
load("GE.RData")
counts=GE
counts$maxVal=apply(counts[,1:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]

counts_WT=counts[,9:17]
counts_NULL=counts[,c(1:8,12:17)]

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
  y=yres[which(yres$padj<0.5),]
  y$ensg=rownames(y)
  y=as.data.frame(y)
  m=merge(y, ens_v75, by.x="ensg", by.y="ensembl_gene_id", all.x= TRUE)
  write.table(m,file=outFile, append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)
  return(m)
} 


conds=c("WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","MT","MT","MT","MT","MT","MT")

conds_NULL=c("WT","WT","WT","WT","WT","WT","WT","WT","MT","MT","MT","MT","MT","MT")
conds_WT=c("WT","WT","WT","MT","MT","MT","MT","MT","MT")

yres=getResults(counts, conds)
see=annotate(yres$res, ens_v75, "mouse_deseq2.tsv")

yres=getResults(counts_NULL, conds_NULL)
see=annotate(yres$res, ens_v75, "mouse_deseq2_NULL.tsv")

yres=getResults(counts_WT, conds_WT)
see=annotate(yres$res, ens_v75, "mouse_deseq2_WT.tsv")

library("pheatmap")
rld <- rlog(yres$dds, blind=FALSE)
ntd <- normTransform(yres$dds)
cc=assay(rld)
zz=transform(cc,SD=apply(cc,1,sd,na.rm=TRUE))
zz_ordered=zz[order((zz$SD),decreasing=TRUE),]
zz_ordered=zz_ordered[1:250,]
zz_ordered=subset(zz_ordered,select=-SD)

df <- as.data.frame(colData(yres$dds)[,c("condition","type")])
annotation1=c("NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","WT","WT","WT","MT","MT","MT","MT","MT","MT")
df$annotation1=annotation1
df <- df[,c("condition","annotation1")]
pheatmap(zz_ordered, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, scale="row",clustering_distance_rows="manhattan", clustering_distance_cols="manhattan", clustering_method="ward")

sampleDists=dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
df$Sample=paste(sub(".1.R1.1","",rownames(df)), df$condition, df$annotation1,sep="-")
colnames(sampleDistMatrix) <- df$Sample

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
pheatmap((sampleDistMatrix), clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, scale="column",col=colors,main="L2 distance between samples", clustering_method="ward")

sampleDists <- dist(t(assay(rld)))
#sampleDists <- log2(sampleDists+0.1)
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition,colnames(rld), sep="-")
colnames(sampleDistMatrix) <- paste(rld$condition, colnames(rld), sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, col=colors, cluster_rows = TRUE, cluster_cols = TRUE)
plotPCA(rld, intgroup=c("condition", "type"))
