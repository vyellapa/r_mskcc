library("biomaRt")
library("RColorBrewer")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ens_v75=getBM(attributes=c('ensembl_gene_id','external_gene_id','chromosome_name'), mart = ensembl_75)

setwd("~/Desktop/TP53_take2/")
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




##### FrameShift and Hostspot ######
ycts=getCts("fshs.mat")
conds=c("H","H","H","H","H","H","H","F","F","F")

yres=getResults(ycts, conds)
see=annotate(yres, ens_v75, "~/Desktop/TP53_take2/HS_vs_FS_nonComplex.tsv")



