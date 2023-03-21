library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(gplots)
mycolors1=colorRampPalette(brewer.pal(name="Spectral", n = 11))(14)
mycolors2=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000') 
mycolors=colorRampPalette(brewer.pal(name="Set1", n = 9))(9)
mycolors=viridis_pal(option = "D")(12)
mycolors=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")



setwd("/ifs/res/leukgen/projects/220/RESULTS/misc/dirichlet_step1/I-H-130720-T1-2-D1-2_DPoutput_10000iters_1000burnin/")
#setwd("/ifs/res/leukgen/projects/220/RESULTS/misc/dirichlet_sep/I-H-130719-T1-2-D1-2_DPoutput_10000iters_1000burnin/")
for(i in list.files(pattern=".*_10000iters_1000burnin_bestConsensusAssignments.bed",recursive=T)) {
  name=gsub("-D1-2_10000iters_1000burnin_bestConsensusAssignments.bed", "",i)
  name=gsub('-','.',name)

  CLUSTER=read.table(i, sep="\t", header=T, stringsAsFactors = FALSE)
  CLUSTER=unique(CLUSTER)
  CLUSTER$merger=paste(CLUSTER$chr, CLUSTER$end, sep=":")
  clusters=as.data.frame(table(CLUSTER$cluster))
  CLUSTER=CLUSTER[!is.na(CLUSTER$cluster),]
  aa=clusters[clusters$Freq>100,"Var1"]
  CLUSTER=CLUSTER[which(CLUSTER$cluster %in% aa),]
  
  colnames(CLUSTER)=paste(name, colnames(CLUSTER), sep="_")
  merged=CLUSTER
    }


for(j in list.files(pattern="*_allDirichletProcessInfo.txt",recursive=T)) {
      
        print(j)
      
      T1=read.table(j, sep="\t", header=T, stringsAsFactors = FALSE)
      T1$merger=paste(T1$chr, T1$end,sep=":")
      jname=gsub("_allDirichletProcessInfo.txt","",j)
      jname=gsub('-','.',jname)
      colnames(T1)=paste(jname,colnames(T1),sep="_") 
      merged=merge(merged, T1,by.x=paste(name, "merger", sep="_"), by.y=paste(jname, "merger", sep="_"), all.x = T )
      
    
  }

m=merged[, colnames(merged) %in% grep("cluster|subclonal.fraction",colnames(merged), value=T)]
mm=m
mm$I.H.130720.T1.2_cluster=factor(mm$I.H.130720.T1.2_cluster)
mm=melt(mm)
mm$variable=gsub("D1.2_subclonal.fraction","",mm$variable)

clust=merged[, colnames(merged) %in% grep("_merger|_cluster|subclonal.fraction",colnames(merged), value=T)]
clust=unique(clust)
clust$I.H.130720.T1.2_merger=make.unique((clust$I.H.130720.T1.2_merger))
cols=clust[,c(1,2)]
rownames(clust)=clust[,1]
#clust=clust[,c(-1,-2)]
sidecols=as.data.frame(cbind(as.data.frame(seq(1:14)), as.data.frame(mycolors1)))
colnames(sidecols)=c("cluster","hex")
#heatmap.2(as.matrix(clust), dendrogram="both",Rowv=TRUE,Colv=TRUE, distfun=function(clust) dist(clust,method = 'manhattan'), hclustfun= function(clust) hclust(clust,method = 'ward'),colsep=c(1), ColSideColors= cd, col = bluered, trace="none")
co=merge(clust,sidecols, by.x="I.H.130720.T1.2_cluster",by.y="cluster", sort = F)
clust=(co[,colnames(co) %in% grep("_merger|subclonal.fraction",colnames(co),value=T)])
sidecols=as.character((co[,colnames(co) %in% grep("hex",colnames(co),value=T)]))
rownames(clust)=clust[,1]
clust=clust[,c(-1)]
#head(clust[10100:10200,])
#head(sidecols[10100:10200,])

#rownames(clust)


colnames(clust)=c('I720.2.L.Arm.Subcut','I720.3.Rib6','I720.4.Sternum','I720.5.Liver','I720.8.R.Lung.Lobe')

png(file = "I130720_d15.png", res=300, width=2000, height=3000)
heatmap.2(as.matrix(clust), dendrogram="row",Rowv=TRUE,Colv=FALSE, scale="row",distfun=function(clust) dist(clust,method = 'euclidean'), hclustfun= function(clust) hclust(clust,method = 'ward'),colsep=c(1), RowSideColors= sidecols, col = bluered, trace="none", margins = c(12,14), labRow = FALSE)
dev.off()

#png(file = "d15_log.png", res=300, width=2000, height=2000)
#heatmap.2(as.matrix(log(clust+0.01)), dendrogram="row",Rowv=TRUE,Colv=FALSE, scale="row",distfun=function(clust) dist(log(clust+0.01),method = 'euclidean'), hclustfun= function(clust) hclust(log(clust+0.01),method = 'ward'),colsep=c(1), RowSideColors= sidecols, col = bluered, trace="none")
#dev.off()


