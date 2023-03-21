setwd("~/Desktop/RNA_TP53")

library( "DESeq" )

counts=read.table("/Users/yellapav/Desktop/RNA_TP53/hotspot.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","Hotspot","Hotspot","Hotspot","WT","WT","Hotspot","WT","WT","WT","WT","WT","WT")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "Hotspot")
write.table(res,file="WT_vs_Hotspot.tsv", append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)






### Missense #####
counts=read.table("/Users/yellapav/Desktop/RNA_TP53/missense.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","missense","WT","WT","missense","missense","missense","missense","WT","missense","WT","missense","WT","missense","missense","WT","WT","missense","missense","WT","WT","WT","missense","missense")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "missense")
write.table(res,file="WT_vs_Missense.tsv", append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)





### Frame Shift and Non-sense #####
counts=read.table("/Users/yellapav/Desktop/RNA_TP53/fs.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","frameshift","WT","frameshift","WT","WT","WT","WT","WT","WT","frameshift","WT","WT","WT")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "frameshift")
write.table(res,file="WT_vs_Truncating.tsv", append=FALSE, sep="\t", eol="\n", row.names=TRUE, quote=F, col.names=TRUE)



### Frame Shift and Hotspot #####
counts=read.table("/Users/yellapav/Desktop/RNA_TP53/HSFS.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("Hotspot","frameshift","Hotspot","Hotspot","Hotspot","frameshift","frameshift","Hotspot","Hotspot","Hotspot","Hotspot","Hotspot","Hotspot","frameshift")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "frameshift", "Hotspot")
write.table(res,file="frameshift_vs_Hotspot.tsv", append=FALSE, sep="\t", eol="\n", row.names=TRUE, quote=F, col.names=TRUE)


### Missense+ Hotspot vs WT #####
counts=read.table("/Users/yellapav/Desktop/RNA_TP53/mut_wt_fr.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>100, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","missense","missense","missense","missense","missense","WT","missense","missense","missense","missense","missense","missense","missense","missense","missense","WT","missense","WT","missense","missense","WT","missense","missense","WT","WT","missense","missense","WT","WT","WT","missense","missense")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "missense")
write.table(res,file="WT_vs_mutated_fr.tsv", append=FALSE, sep="\t", eol="\n", row.names=TRUE, quote=F, col.names=TRUE)



### Antisense ####


counts=read.table("/Users/yellapav/Desktop/RNA_TP53/hotspot_fr.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","WT","Hotspot","Hotspot","Hotspot","Hotspot","Hotspot","WT","WT","Hotspot","WT","WT","WT","WT","WT","WT")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "Hotspot")
write.table(res,file="WT_vs_Hotspot_fr.tsv", append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)



### Missense #####
counts=read.table("/Users/yellapav/Desktop/RNA_TP53/missense_fr.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","missense","WT","WT","missense","missense","missense","missense","WT","missense","WT","missense","WT","missense","missense","WT","WT","missense","missense","WT","WT","WT","missense","missense")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "missense")
write.table(res,file="WT_vs_Missense_fr.tsv", append=FALSE, sep="\t", eol="\n", quote=F, row.names=TRUE, col.names=TRUE)




### Frame Shift and Non-sense #####
counts=read.table("/Users/yellapav/Desktop/RNA_TP53/fs_fr.mat", header=T, sep="\t",stringsAsFactors=FALSE)
counts$maxVal=apply(counts[,2:dim(counts)[2]], 1, max)
counts=counts[counts$maxVal>10, ]
counts=counts[,-dim(counts)[2]]
rownames(counts)=counts[,1]
counts=counts[,-1]


## DESEQ
conds=c("WT","frameshift","WT","frameshift","WT","WT","WT","WT","WT","WT","frameshift","WT","WT","WT")
dcountAll=counts
dcountDesign=data.frame(row.names = colnames(dcountAll), condition=conds, libType = c(rep("paired",ncol(dcountAll))))
countTable=dcountAll
condition=dcountDesign$condition
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
#cds=counts(cds,normalized=TRUE) 
cds=estimateDispersions(cds)
res=nbinomTest(cds, "WT", "frameshift")
write.table(res,file="WT_vs_Truncating_fr.tsv", append=FALSE, sep="\t", eol="\n", row.names=TRUE, quote=F, col.names=TRUE)





