library(reshape2)
library(ggplot2)
library(DESeq)
library("biomaRt")
setwd("~/Desktop")
trans=read.table("~/Desktop/MMRF/MMRF_CoMMpass_IA9b_LongInsert_Canonical_Ig_Translocations.txt", header=T, sep="\t")
cyto=read.table("~/Desktop/MMRF/MMRF_CoMMpass_IA9b_CNA_LongInsert_FISH_CN_Unique_Visits.txt", header=T, sep="\t")
#maf=read.table("~/Desktop/MMRF/MMRF_Canonical_GLFilter_uniq.maf", sep = "\t", stringsAsFactors=FALSE, encoding= "utf-8")



fpkm=read.table("~/Desktop/MMRF/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt", header=T, sep="\t",stringsAsFactors=FALSE)
counts=read.table("~/Desktop/MMRF/MMRF_CoMMpass_IA9_E74GTF_HtSeq_Gene_Counts.txt", header=T, sep="\t",stringsAsFactors=FALSE)
#counts$rowsums=rowSums(counts[,2:dim(counts)[2]])
#counts=counts[counts$rowsums>49, ]
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]

rownames(counts)=counts[,1]
counts=counts[,-1]


newFpkm=as.data.frame(t(fpkm[fpkm$GENE_ID %in% rownames(counts), ]))
#geneId=c("ENSG00000110092","ENSG00000002586","ENSG00000110651","ENSG00000057657")
#newFpkm=as.data.frame(t(fpkm[fpkm$GENE_ID %in% geneId, ]))
#newFpkm$GENE_ID=row.names(newFpkm)
#colnames(newFpkm)=c("ENSG00000002586","ENSG00000057657","ENSG00000110092","ENSG00000110651","GENE_ID")
#colnames(newFpkm)=c("ENSG00000002586","ENSG00000110092","GENE_ID")
colnames(newFpkm)=apply(newFpkm[1,],1,as.character)
newFpkm=(newFpkm[c(-1,-2),])



transPos=trans[trans[,"SeqWGS_CCND1_CALL"]>0, c("Study_Visit_ID","SeqWGS_CCND1_CALL")]
transNeg=trans[trans[,"SeqWGS_CCND1_CALL"]==0, c("Study_Visit_ID","SeqWGS_CCND1_CALL")]
tPosId=transPos[,"Study_Visit_ID"]
tNegId=transNeg[,"Study_Visit_ID"]

## DESEQ
dcountNeg=counts[,colnames(counts) %in% tNegId ]
dcountNeg=dcountNeg[,1:10]
dcountPos=counts[,colnames(counts) %in% tPosId ]
dcountPos=dcountPos[,1:10]
n=rep("NEG",ncol(dcountNeg))
p=rep("POS",ncol(dcountPos))
conds=c(n,p)
dcountAll=cbind(dcountNeg,dcountPos)
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
library( "DESeq" )
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "NEG", "POS")

ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
all75=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id', values = ensg, mart = ensembl_75)
newRes=merge(res,all75, by.x=c("id"), by.y=c("ensembl_gene_id"))

newFpkmPos=newFpkm[rownames(newFpkm) %in% tPosId, ]
newFpkmNeg=newFpkm[rownames(newFpkm) %in% tNegId, ]

temp=sapply(newFpkmPos, function(x) {log2(as.numeric(as.character(x))+1) })
row.names(temp)=row.names(newFpkmPos)
newFpkmPos=temp
CASE=rep("CCND1-POS",nrow(newFpkmPos))
newFpkmPos=cbind.data.frame(newFpkmPos,CASE)


temp=sapply(newFpkmNeg, function(x) {log2(as.numeric(as.character(x))+1) })
row.names(temp)=row.names(newFpkmNeg)
newFpkmNeg=temp
CASE=rep("CCND1-NEG",nrow(newFpkmNeg))
newFpkmNeg=cbind.data.frame(newFpkmNeg,CASE)

plot=rbind.data.frame(newFpkmNeg,newFpkmPos)
melted_plot=melt(plot)
melted_plot$variable=as.character(melted_plot$variable)
melted_plot$variable[melted_plot$variable=="ENSG00000110092"]="CCND1"
melted_plot$variable[melted_plot$variable=="ENSG00000002586"]="CD99"
melted_plot$variable[melted_plot$variable=="ENSG00000057657"]="PRDM1"
melted_plot$variable[melted_plot$variable=="ENSG00000110651"]="CD81"

ggplot(melted_plot, aes(variable, value,fill=CASE))+geom_boxplot(alpha=0.6)+geom_point(alpha=0.3, position = position_dodge(width=0.75))+scale_fill_brewer(palette = "Set1")+theme_bw()+ylab("log 2 FPKM")+xlab("")


## MAF
ns=c("Missense_Mutation","Splice_Site","Frame_Shift_Del","Nonsense_Mutation","Frame_Shift_Ins","Translation_Start_Site","Nonstop_Mutation")
nsmaf=maf[maf$Variant_Classification %in% ns,]
primSamps=unique(nsmaf[grep("MMRF_...._1_BM",nsmaf$Tumor_Sample_Barcode,perl=TRUE),c(16)])
nsmafPrim=nsmaf[nsmaf$Tumor_Sample_Barcode %in% primSamps,]
#Check which BCL genes are present
unique(nsmaf[grep("BCL",nsmaf$Hugo_Symbol,perl=TRUE),c(1)])



nfpkm=fpkm[,c(-2)]
nfpkm=nfpkm[!duplicated(nfpkm$GENE_ID), ]
rownames(nfpkm)=nfpkm[,1]
nfpkm=nfpkm[,c(-1)]
l=log2(nfpkm+1)


#Explore
fpkmUniq=fpkm[!duplicated(fpkm$GENE_ID), ]
rownames(fpkmUniq)=fpkmUniq$GENE_ID
fpkmUniq=fpkmUniq[,c(-1,-2)]
fpkmUniq=as.data.frame(t(log2(fpkmUniq+1)))
fpkmUniq$SAMPLE=rownames(fpkmUniq)

fpkmUniqPos=fpkmUniq[fpkmUniq$SAMPLE %in% tPosId, ]
fpkmUniqNeg=fpkmUniq[fpkmUniq$SAMPLE %in% tNegId, ]

tttest <- function(gene){
  library("biomaRt")
  p=fpkmUniqPos
  n=fpkmUniqNeg
  CASE=rep("CCND1-NEG",nrow(n))
  n=cbind.data.frame(n,CASE)
  CASE=rep("CCND1-POS",nrow(p))
  p=cbind.data.frame(p,CASE)
  
  
  plot=rbind.data.frame(n,p)
  
  #plot=plot[,which(colnames(plot)==gene)]
  #cat(head(plot))
  #print(gene)
  x=t.test(plot[,gene]~plot$CASE)
  plot=cbind.data.frame(plot[,gene],plot$CASE)
  m=melt(plot)
  m$variable=gene
  setwd("/Users/yellapav/Desktop/MM_CD99/figs")
  if((x$p.value<1) & (abs(x$estimate[1]-x$estimate[2])>0)) {
  #if((x$p.value<0.0001) & (abs(x$estimate[1]-x$estimate[2])>1)) {
    print(x$p.value)
    ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    all75=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id', values = gene, mart = ensembl_75)
    hgnc=all75$hgnc_symbol
  colnames(m)=c("CASE","variable","value")
  val=paste("t-test p-value",x$p.value,sep="=")
  p=ggplot(m, aes(variable, value,fill=CASE))+geom_boxplot(alpha=0.6)+geom_point(alpha=0.3, position = position_dodge(width=0.75))+scale_fill_brewer(palette = "Set1")+theme_bw()+ylab("log 2 FPKM")+xlab(hgnc)+ylim(c(0,14))+annotate("text", x = 1, y = 13, label = val)
  filePdf=paste(hgnc,"pdf",sep=".")
  ggsave(p, file=filePdf)
  }
  return(x)
}







