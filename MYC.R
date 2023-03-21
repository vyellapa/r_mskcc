library(reshape2)
library(ggplot2)
library(DESeq)
setwd("~/Desktop/MM_CD99")
trans=read.table("data/MMRF_CoMMpass_IA9b_LongInsert_Canonical_Ig_Translocations.txt", header=T, sep="\t")
cyto=read.table("data/MMRF_CoMMpass_IA9b_CNA_LongInsert_FISH_CN_Unique_Visits.txt", header=T, sep="\t")
#fpkm=read.table("data/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt", header=T, sep="\t",stringsAsFactors=FALSE)
counts=read.table("data/MMRF_CoMMpass_IA9_E74GTF_HtSeq_Gene_Counts.txt", header=T, sep="\t",stringsAsFactors=FALSE)

counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


#newFpkm=as.data.frame(t(fpkm[fpkm$GENE_ID %in% rownames(counts), ]))
#colnames(newFpkm)=apply(newFpkm[1,],1,as.character)
#newFpkm=(newFpkm[c(-1,-2),])



transPos=trans[trans[,"SeqWGS_MYC_CALL"]>0, c("Study_Visit_ID","SeqWGS_MYC_CALL")]
transNeg=trans[trans[,"SeqWGS_MYC_CALL"]==0, c("Study_Visit_ID","SeqWGS_MYC_CALL")]
tPosId=transPos[,"Study_Visit_ID"]
tNegId=transNeg[,"Study_Visit_ID"]

## DESEQ
dcountNeg=counts[,colnames(counts) %in% tNegId ]
dcountNeg=dcountNeg[,1:130]
dcountPos=counts[,colnames(counts) %in% tPosId ]

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
write.table(res,file="SeqWGS_MYC_CALL.tsv", append=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=TRUE)
