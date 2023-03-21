setwd("~/Documents/Work/Postdoc/projects/Eva")

sample.list <- c( "W0017438F", # TP53 NULL 
                  "W0017439F", # TP53 NULL
                  "W0017440F", # TP53 NULL
                  "W0017441F", # TP53 NULL
                  "W0017442F", # TP53 NULL
                  "W0017443F", # TP53 NULL
                  "W0017444F", # TP53 NULL
                  "W0017445F", # TP53 NULL
                  "W0017446F", # TP53 WT 
                  "W0017447F", # TP53 WT
                  "W0017448F", # TP53 WT
                  "W0017449F", # TP53 MT  
                  "W0017450F", # TP53 MT
                  "W0017451F", # TP53 MT
                  "W0017452F", # TP53 MT
                  "W0017453F", # TP53 MT
                  "W0017454F") # TP53 MT


# CountData
# GE <- Reduce("cbind", lapply(1:length(sample.list), function(n)
#                                   {
#                                           GE <- read.table(paste0("data/output/", sample.list[n], "_ReadsPerGene.out.tab"), row.names=1, header=F, sep="\t")
# 
#                                           return(GE[-(1:4),1,drop=F])
#                                   }))
# colnames(GE) <- sample.list
# 
# save(GE, file="data/processed/GE.RData")

GE <- get(load("./data/processed/GE.RData"))

# Color scales and other
library(gplots)
library(RColorBrewer)
library(scales)
colors.condition <- c("NULL"="#7fc97f",
                  "WT"= "#beaed4",
                  "MT"= "#fdc086")
colmap <- greenred

# Infos about the samples
assoc <- read.csv("data/raw/info/assoc.csv", header=T)
Batch <- as.character(assoc$Batch)
Condition <- rep(c("NULL","WT","MT"), c(8,3,6))
SampleNames <- as.character(assoc[,5])

# Verify that TP53 is not expressed for the first samples
TP53 <- "ENSMUSG00000059552"
TP53 %in% rownames(GE)
GE[TP53,]

TP53.df <- data.frame(TP53= t(as.vector(GE[TP53,])), Condition=Condition, Batch=Batch, Sample=SampleNames)
TP53.df$Sample <- factor(TP53.df$Sample, levels=SampleNames)
TP53.df$Condition <- factor(TP53.df$Condition, levels=c("NULL","WT","MT"))
colnames(TP53.df)[1] <- "TP53"

library(ggplot2)
pdf("results/TP53.pdf")
print(ggplot(TP53.df) + geom_point(aes(x=Sample,y=TP53, colour=Condition,shape=Batch), size=3) + scale_colour_manual(values=colors.condition) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(0,5000))
dev.off()

pdf("results/TP53_boxplot.pdf")
print(ggplot(TP53.df) + geom_boxplot(aes(x=Condition,y=TP53)) + geom_jitter(aes(x=Condition,y=TP53,colour=Condition)) + scale_colour_manual(values=colors.condition)  + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(0,5000))
dev.off()

# FOXH1
FOXH1 <- "ENSMUSG00000033837"
FOXH1 %in% rownames(GE)
GE[FOXH1,]

FOXH1.df <- data.frame(FOXH1= t(as.vector(GE[FOXH1,])), Condition=Condition, Batch=Batch, Sample=SampleNames)
FOXH1.df$Sample <- factor(FOXH1.df$Sample, levels=SampleNames)
FOXH1.df$Condition <- factor(FOXH1.df$Condition, levels=c("NULL","WT","MT"))

colnames(FOXH1.df)[1] <- "FOXH1"

library(ggplot2)
pdf("results/FOXH1.pdf")
print(ggplot(FOXH1.df) + geom_point(aes(x=Sample,y=FOXH1, colour=Condition,shape=Batch), size=3) + scale_colour_manual(values=colors.condition) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) )
dev.off()

pdf("results/FOXH1_boxplot.pdf")
print(ggplot(FOXH1.df) + geom_boxplot(aes(x=Condition,y=FOXH1)) + geom_jitter(aes(x=Condition,y=FOXH1,colour=Condition)) + scale_colour_manual(values=colors.condition)  + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) )
dev.off()

# Do a regular PCA to see if data cluster
dat <- as.matrix(GE)
filter.sd <- which(apply(dat,1,sd)==0)
dat.filter <- dat[-filter.sd,]
# pc.cr <-  princomp(log(dat.filter + 0.01), center=T, scale.=T)
pc.cr <-  princomp(log(dat.filter + 0.01), cor=T)

pc.df <- data.frame(PC1=pc.cr$loadings[,1], PC2=pc.cr$loadings[,2], Condition=Condition, Batch=Batch, Sample=SampleNames)
pc.df$Sample <- factor(pc.df$Sample, levels=SampleNames)
pc.df$Condition <- factor(pc.df$Condition, levels=c("NULL","WT","MT"))

# Add PC values
# Eigenvalues
eig <- (pc.cr$sdev)^2
# Variances in percentage
variance <-eig*100/sum(eig)
#

pdf("results/PCA_full.pdf")
# print(ggplot(pc.df) + geom_point(aes(x=PC1,y=PC2, colour=Condition), size=3) + scale_colour_manual(values=colors.condition) +  theme_bw() + xlab(paste0("PC1 (Variance Explained=", floor(variance[1]) , " %)" )) + ylab(paste0("PC2 (Variance Explained=", floor(variance[2]) , " %)" )))
print(ggplot(pc.df, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Condition), size=3) + geom_text(nudge_x=0.0003, nudge_y=0.001, aes(label=Sample)) + scale_colour_manual(values=colors.condition) +  theme_bw() + xlab(paste0("PC1 (Variance Explained=", floor(variance[1]) , " %)" )) + ylab(paste0("PC2 (Variance Explained=", floor(variance[2]) , " %)" )))
dev.off()

# Do a regular PCA to see if data cluster
filter.TP53 <- Condition %in% c("NULL","MT")
GE.final <- GE[,filter.TP53]

info.df <- data.frame(Batch= Batch[filter.TP53], Condition=Condition[filter.TP53], SampleNames=SampleNames[filter.TP53])
rownames(info.df) <- sample.list[filter.TP53]

dat <- as.matrix(GE.final)
filter.sd <- which(apply(dat,1,sd)==0)
dat.filter <- dat[-filter.sd,]
pc.cr <-  princomp(log(dat.filter + 0.01), center=T, scale.=T)
# pc.cr <-  princomp(log(dat.filter + 0.01), cor=T)

pc.df <- data.frame(PC1=pc.cr$loadings[,1], PC2=pc.cr$loadings[,2], Condition=info.df$Condition, Batch=info.df$Batch, Sample=info.df$SampleNames)
pc.df$Sample <- factor(pc.df$Sample, levels=SampleNames)
pc.df$Condition <- factor(pc.df$Condition, levels=c("NULL","MT"))

eig <- (pc.cr$sdev)^2
# Variances in percentage
variance <-eig*100/sum(eig)
#

pdf("results/PCA_NULL_vs_MT.pdf")
print(ggplot(pc.df) + geom_point(aes(x=PC1,y=PC2, colour=Condition), size=3) + scale_colour_manual(values=colors.condition) +  theme_bw() + xlab(paste0("PC1 (Variance Explained=", round(variance[1],1) , " %)" )) + ylab(paste0("PC2 (Variance Explained=", round(variance[2],1) , " %)" )))
dev.off()

# Do a regular PCA to see if data cluster
filter.TP53 <- Condition %in% c("WT","MT")
GE.final <- GE[,filter.TP53]

info.df <- data.frame(Batch= Batch[filter.TP53], Condition=Condition[filter.TP53], SampleNames=SampleNames[filter.TP53])
rownames(info.df) <- sample.list[filter.TP53]

dat <- as.matrix(GE.final)
filter.sd <- which(apply(dat,1,sd)==0)
dat.filter <- dat[-filter.sd,]
pc.cr <-  princomp(log(dat.filter + 0.01), center=T, scale.=T)
# pc.cr <-  princomp(log(dat.filter + 0.01), cor=T)

eig <- (pc.cr$sdev)^2
# Variances in percentage
variance <-eig*100/sum(eig)
#

pc.df <- data.frame(PC1=pc.cr$loadings[,1], PC2=pc.cr$loadings[,2], Condition=info.df$Condition, Batch=info.df$Batch, Sample=info.df$SampleNames)
pc.df$Sample <- factor(pc.df$Sample, levels=SampleNames)
pc.df$Condition <- factor(pc.df$Condition, levels=c("WT","MT"))

pdf("results/PCA_WT_vs_MT.pdf")
print(ggplot(pc.df) + geom_point(aes(x=PC1,y=PC2, colour=Condition), size=3) + scale_colour_manual(values=colors.condition) +  theme_bw() + xlab(paste0("PC1 (Variance Explained=", floor(variance[1]) , " %)" )) + ylab(paste0("PC2 (Variance Explained=", floor(variance[2]) , " %)" )))
dev.off()

# Differential GE on TP53 WT and MT
######################################
source("src/lib/run_DESeq.R")

# filter list
filter.list <- list(TP53=c("WT","MT"),
                    non_mut=c("WT","NULL"),
                   no_WT= c("NULL", "MT"))
names(filter.list) <- c("WT_vs_MT", "WT_vs_NULL", "NULL_vs_MT")

library(DESeq2)

res <- lapply(1:length(filter.list), function(n)
              {
                      out <- run_DESeq(GE, filter.list[[n]], Batch, Condition, out.name=names(filter.list)[n])
              })

# comparison MT vs NULL
gene.list <- lapply(1:length(res), function(n)
                    {
                            genes <- rownames(res[[n]]$DESeq)[which(res[[n]]$DESeq$padj < 0.05)]
                    })
names(gene.list) <- names(filter.list)

source("./src/convert_gene_ids.R")





library(VennDiagram)
venn.out <- venn.diagram(gene.list, filename=NULL)

pdf("results/venn.pdf")
grid.draw(venn.out)
dev.off()

### differentially expressed in WT vs MT but also WT vs NULL:
gene.interest <- intersect(gene.list[[1]], gene.list[[3]])

FOXH1 <- "ENSMUSG00000033837"







#
FOXH1 %in% gene.interest
FOXH1 %in% gene.list[[3]]

# WT_vs_MT
GE.norm <- counts(res[[3]]$dds)
colors.groups <- c("NULL"="#7fc97f",
                   #"WT"= "#beaed4",
                   "MT"= "#fdc086")
colmap <- greenred

# Full genes
genes.notzero <- which(res[[3]]$DESeq$baseMean > quantile(res[[3]]$DESeq$baseMean,0.9))

GE.full <- GE.norm[genes.notzero, ]
sd.infos <- apply(GE.full,1,sd)
GE.full <- GE.full[which(sd.infos !=0),]
GE.scaled <- scale(t(GE.full))

#####

pdf("results/heatmap_gene_interest_WT_NULL_full.pdf")
#p <- heatmap.2(as.matrix(GE.filtered), scale="row", dendrogram="row",hclustfun=function(x){hclust(x,method="complete")},
p <- heatmap.2(t(GE.scaled), dendrogram="col",hclustfun=function(x){hclust(x,method="complete")},
          trace="none", labCol=NA, key=T, labRow=NA,
          ColSideColors=colors.groups[Condition[Condition %in% filter.list[["NULL_vs_MT"]]]],
          col=colmap)
legend("bottomleft", names(colors.groups), fill=colors.groups)

dev.off()


# DE Genes
GE.full <- GE.norm[gene.interest,]
# GE.filtered <- GE[gene.interest,Condition %in% filter.list[["NULL_vs_MT"]]]

GE.bis <- scale(t(GE.full))

pdf("results/heatmap_gene_interest_WT_NULL.pdf")
#p <- heatmap.2(as.matrix(GE.filtered), scale="row", dendrogram="row",hclustfun=function(x){hclust(x,method="complete")},
p <- heatmap.2(t(GE.bis), dendrogram="col",hclustfun=function(x){hclust(x,method="complete")},
          #trace="none", labCol=NA, key=T, labRow=NA,
          trace="none", labCol=SampleNames[Condition %in% filter.list[["NULL_vs_MT"]]], key=T, labRow=NA,
          ColSideColors=colors.groups[Condition[Condition %in% filter.list[["NULL_vs_MT"]]]],
          col=colmap)
legend("bottomleft", names(colors.groups), fill=colors.groups)
dev.off()

# VOLCANO PLOT
gene.interest <- intersect(gene.list[[1]], gene.list[[3]])

dat.df <- data.frame(log2FC=res[[3]]$DESeq$log2FoldChange, log10p=log10(res[[3]]$DESeq$padj), genes= rownames(res[[3]]$DESeq), gene_symbol= res.symbol.full[[3]])

GO.up <- read.table("./GO/GO_up.txt", header=T, sep="\t")
GO.k <- 1

GO.genes <- strsplit(as.character(GO.up[GO.k, "geneID"]), "/")[[1]]

cols <- rep("grey", nrow(dat.df))
cols[dat.df$gene_symbol %in% GO.genes] <- "green"

dat.df$cols <- cols

# cols <- sapply(1:nrow(dat.df), function(n)
#                {
#                        if (is.na(dat.df$log2FC[n])|is.na(dat.df$log10p[n]))
#                        {
#                                return("grey")
#                        } else if ((abs(dat.df$log2FC[n])>=2)&(dat.df$log10p[n]<log10(0.05)))
#                        {
#                                return("green")
#                        } else if ((abs(dat.df$log2FC[n])<2)&(dat.df$log10p[n]<log10(0.05)))
#                        {
#                                return("red")
#                        } else if ((abs(dat.df$log2FC[n])>=2)&(dat.df$log10p[n]>log10(0.05)))
#                        {
#                                return("orange")
#                        } else
#                        {
#                                return("black")
#                        }
#                })
# dat.df$cols <- cols

dat.df$alpha <- 1
dat.df$alpha[dat.df$cols=="grey"] <- 0.6

cols.manual <- c("grey"= "grey",
                 "green" = "#91cf60",
                 "red" = "#fc8d59",
                 "orange" = "#ffffbf",
                 "black" = "#999999")

### Remove pseudo-genes
ps.index <- grep("ENS", dat.df$gene_symbol)

dat.df.bis <- dat.df[ -ps.index,]




library(ggplot2)
pdf("results/volcano_gene_interest_WT_NULL.pdf")
#print(ggplot(dat.df) + geom_point(aes(x=log2FC,y=-log10p,color=cols)) + scale_color_manual(values=cols.manual) + theme_bw() + theme(legend.position="none"))
print(ggplot(dat.df) + geom_point(aes(x=log2FC,y=-log10p,color=cols)) + geom_text(data=dat.df[which(dat.df$log10p< -20),   ], aes(x=log2FC + 0.05, y=-log10p + 1, label= gene_symbol, color=cols)) + scale_color_manual(values=cols.manual) + theme_bw() + theme(legend.position="none"))
dev.off()

library(ggplot2)
pdf("results/volcano_gene_interest_WT_NULL_GO.pdf")
#print(ggplot(dat.df.bis) + geom_point(aes(x=log2FC,y=-log10p,color=cols)) + geom_text(data=dat.df.bis[which(dat.df.bis$cols=="green"),   ], aes(x=log2FC + 0.05, y=-log10p + 1, label= gene_symbol, color=cols)) + scale_color_manual(values=cols.manual) + theme_bw() + theme(legend.position="none"))
print(ggplot(dat.df.bis) + geom_point(aes(x=log2FC,y=-log10p,color=cols, alpha=alpha))  + scale_color_manual(values=cols.manual) + theme_bw() + theme(legend.position="none"))
dev.off()


# pp <- which.min(dat.df$log10p)

# rownames(res[[3]]) <- res.symbol.dup[[3]]
res[[3]]$DESeq$Symbol <- res.symbol.full[[3]]

res.bis <- res[[3]]$DESeq[ -ps.index, ]

down.genes <- res[[3]]$DESeq[which(res[[3]]$DESeq$log2FoldChange <0 & res[[3]]$DESeq$padj < 0.05),]
up.genes <- res[[3]]$DESeq[which(res[[3]]$DESeq$log2FoldChange >0 & res[[3]]$DESeq$padj < 0.05),]
all.genes <- res[[3]]$DESeq[which(res[[3]]$DESeq$padj < 0.05),]

down.genes <- down.genes[order(abs(down.genes$log2FoldChange),decreasing=T),]
up.genes <- up.genes[order(abs(up.genes$log2FoldChange),decreasing=T),]
all.genes <- all.genes[order(abs(all.genes$log2FoldChange),decreasing=T),]

write.table(down.genes, file=paste0("results/", names(filter.list)[3] , "_down.txt"), col.names=T, row.names=T, quote=F)
write.table(up.genes, file=paste0("results/", names(filter.list)[3] , "_up.txt"), col.names=T, row.names=T, quote=F)
write.table(all.genes, file=paste0("results/", names(filter.list)[3] , "_all.txt"), col.names=T, row.names=T, quote=F)

### No Pseudo genes
down.genes <- res.bis[which(res.bis$log2FoldChange <0 & res.bis$padj < 0.05),]
up.genes <- res.bis[which(res.bis$log2FoldChange >0 & res.bis$padj < 0.05),]
all.genes <- res.bis[which(res.bis$padj < 0.05),]

down.genes <- down.genes[order(abs(down.genes$log2FoldChange),decreasing=T),]
up.genes <- up.genes[order(abs(up.genes$log2FoldChange),decreasing=T),]
all.genes <- all.genes[order(abs(all.genes$log2FoldChange),decreasing=T),]

write.table(down.genes, file=paste0("results/", names(filter.list)[3] , "_down_no_pseudo_genes.txt"), col.names=T, row.names=T, quote=F)
write.table(up.genes, file=paste0("results/", names(filter.list)[3] , "_up_no_pseudo_genes.txt"), col.names=T, row.names=T, quote=F)
write.table(all.genes, file=paste0("results/", names(filter.list)[3] , "_all_no_pseudo_genes.txt"), col.names=T, row.names=T, quote=F)


################
# WT vs NULL
n <- 2
GE.norm <- counts(res[[n]]$dds)
colors.groups <- c("NULL"="#7fc97f",
                   "WT"= "#beaed4")
                   #"MT"= "#fdc086")
colmap <- greenred

# Full genes
genes.filter <- gene.list[[n]]
GE.full <- GE.norm[genes.filter, ]
sd.infos <- apply(GE.full,1,sd)
GE.full <- GE.full[which(sd.infos !=0),]
GE.scaled <- scale(t(GE.full))

# DE Genes

pdf("results/heatmap_gene_interest_WT_NULL.pdf")
#p <- heatmap.2(as.matrix(GE.filtered), scale="row", dendrogram="row",hclustfun=function(x){hclust(x,method="complete")},
p <- heatmap.2(t(GE.scaled), dendrogram="col",hclustfun=function(x){hclust(x,method="complete")},
          trace="none", labCol=NA, key=T, labRow=NA,
          ColSideColors=colors.groups[Condition[Condition %in% filter.list[[n]]]],
          col=colmap)
legend("bottomleft", names(colors.groups), fill=colors.groups)

dev.off()

# List genes 
res[[n]]$DESeq$Symbol <- res.symbol.full[[n]]

down.genes <- res[[n]]$DESeq[which(res[[n]]$DESeq$log2FoldChange <0 & res[[n]]$DESeq$padj < 0.05),]
up.genes <- res[[n]]$DESeq[which(res[[n]]$DESeq$log2FoldChange >0 & res[[n]]$DESeq$padj < 0.05),]
all.genes <- res[[n]]$DESeq[which(res[[n]]$DESeq$padj < 0.05),]

down.genes <- down.genes[order(abs(down.genes$log2FoldChange),decreasing=T),]
up.genes <- up.genes[order(abs(up.genes$log2FoldChange),decreasing=T),]
all.genes <- all.genes[order(abs(all.genes$log2FoldChange),decreasing=T),]


write.table(down.genes, file=paste0("results/", names(filter.list)[n] , "_down.txt"), col.names=T, row.names=T, quote=F)
write.table(up.genes, file=paste0("results/", names(filter.list)[n] , "_up.txt"), col.names=T, row.names=T, quote=F)
write.table(all.genes, file=paste0("results/", names(filter.list)[n] , "_all.txt"), col.names=T, row.names=T, quote=F)

###



# # SampleNames[p$colInd]
# SampleNames[p$colInd]
# 
# # Order over expressed and under-expressed
# res.filtered <- res[[1]][[1]][gene.interest,]
# neg.Gene <- which(res.filtered$log2FoldChange > 0)
# pos.Gene <- which(res.filtered$log2FoldChange < 0)
# 
# order.genes <- c(pos.Gene, neg.Gene)
# 
# pdf("results/heatmap_gene_interest_order.pdf")
# p <- heatmap.2(as.matrix(GE.filtered[order.genes,]), hclustfun=function(x){hclust(x,method="complete")},
#                scale="row", Rowv=NULL, dendrogram="none",
#           trace="none", labCol=NA, key=T, labRow=NA,
#           ColSideColors=colors.groups[Condition],
#           col=colmap)
# legend("topright", names(colors.groups), fill=colors.groups)
# dev.off()
# 
# genelist.neg <- res.filtered[which(res.filtered$log2FoldChange > 0),]
# genelist.pos <- res.filtered[which(res.filtered$log2FoldChange < 0),]
# 
# FOXH1 %in% rownames(genelist.pos)
# 
# nrow(genelist.pos)
# nrow(genelist.neg)
# 
# write.table(rownames(genelist.pos), file="results/genepos.txt", col.names=F, row.names=F, quote=F, sep="\t")

## ### differentially expressed in WT vs MT but not in WT vs NULL:
## gene.interest <- intersect(gene.list[[1]], gene.list[[3]])
## 
## FOXH1 <- "ENSMUSG00000033837"
## 
## #
## FOXH1 %in% gene.interest
## 
## #
## GE.norm <- counts(dds)
## 
## library(gplots)
## library(RColorBrewer)
## colors.groups <- c("NULL"="#7fc97f",
##                   "WT"= "#beaed4",
##                   "MT"= "#fdc086")
## colmap <- greenred
## 
## GE.filtered <- GE[gene.interest,]
## 
## pdf("results/heatmap_gene_interest.pdf")
## p <- heatmap.2(as.matrix(GE.filtered), scale="row", dendrogram="row",hclustfun=function(x){hclust(x,method="complete")},
##           trace="none", labCol=NA, key=T, labRow=NA,
##           ColSideColors=colors.groups[Condition],
##           col=colmap)
## legend("topright", names(colors.groups), fill=colors.groups)
## dev.off()
## 
## # SampleNames[p$colInd]
## SampleNames[p$colInd]
## 
## # Order over expressed and under-expressed
## res.filtered <- res[[1]][gene.interest,]
## neg.Gene <- which(res.filtered$log2FoldChange > 0)
## pos.Gene <- which(res.filtered$log2FoldChange < 0)
## 
## order.genes <- c(pos.Gene, neg.Gene)
## 
## pdf("results/heatmap_gene_interest_order.pdf")
## p <- heatmap.2(as.matrix(GE.filtered[order.genes,]), hclustfun=function(x){hclust(x,method="complete")},
##                scale="row", Rowv=NULL, dendrogram="none",
##           trace="none", labCol=NA, key=T, labRow=NA,
##           ColSideColors=colors.groups[Condition],
##           col=colmap)
## legend("topright", names(colors.groups), fill=colors.groups)
## dev.off()
## 
## genelist.neg <- res.filtered[which(res.filtered$log2FoldChange > 0),]
## genelist.pos <- res.filtered[which(res.filtered$log2FoldChange < 0),]
## 
## FOXH1 %in% rownames(genelist.pos)
## 
## nrow(genelist.pos)
## nrow(genelist.neg)
## 
## write.table(rownames(genelist.pos), file="results/genepos.txt", col.names=F, row.names=F, quote=F, sep="\t")




########################################################################
# FOXH1
FOXH1 <- "ENSMUSG00000033837"

# Normalized counts
GE.norm <- counts(dds)

info.df <- data.frame(Batch= Batch[filter.TP53], Condition=Condition[filter.TP53], SampleNames=SampleNames[filter.TP53])
rownames(info.df) <- sample.list[filter.TP53]

dat <- as.matrix(GE.norm)
filter.sd <- which(apply(dat,1,sd)==0)
dat.filter <- dat[-filter.sd,]
pc.cr <-  princomp(log(dat.filter + 0.01), center=T, scale.=T)

pc.df <- data.frame(PC1=pc.cr$loadings[,1], PC2=pc.cr$loadings[,2], Condition=info.df$Condition, Batch=info.df$Batch, Sample=info.df$SampleNames)
pc.df$Sample <- factor(pc.df$Sample, levels=SampleNames)
pc.df$Condition <- factor(pc.df$Condition, levels=c("WT","MT"))

pdf("results/PCA_GE_normalized.pdf")
print(ggplot(pc.df) + geom_point(aes(x=PC1,y=PC2, colour=Condition,shape=Batch)) + theme_bw() )
dev.off()



