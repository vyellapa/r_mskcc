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



setwd("/ifs/res/leukgen/projects/220/RESULTS/misc/dirichlet_sep/I-H-130718/")
j="I-H-130718-T1-1-D1-2_allDirichletProcessInfo.txt"
T1=read.table(j, sep="\t", header=T, stringsAsFactors = FALSE)
T1$merger=paste(T1$chr, T1$start,T1$end,sep=":")
merged=as.data.frame(T1$merger)
colnames(merged)=c("merger")
for(j in list.files(pattern="*_allDirichletProcessInfo.txt",recursive=T)) {
      
        print(j)
      T1=read.table(j, sep="\t", header=T, stringsAsFactors = FALSE)
      T1$merger=paste(T1$chr,T1$start, T1$end,sep=":")
      jname=gsub("_allDirichletProcessInfo.txt","",j)
      jname=gsub('-','.',jname)
      colnames(T1)=paste(jname,colnames(T1),sep="_") 
      merged=merge(merged, T1,by.x="merger", by.y=paste(jname, "merger", sep="_"), all.x = T )
    }

m=merged[, colnames(merged) %in% grep("merger|subclonal.fraction",colnames(merged), value=T)]
mm=m
rownames(mm)=make.unique(as.character(mm[,1]))
mm=mm[,c(-1)]
colnames(mm)=c("I718.1.L.Calf","I718.10.Gallbladder","I718.2.L.Inguinal","I718.4.R.Pleural","I718.6.R.Lung","I718.9.Liver")
zz=mm
zz[is.na(zz)] <- 0
zz$maxVal=apply(zz[,1:dim(zz)[2]], 1, max)
zz=zz[zz$maxVal>0,]
zz=zz[,c(-7)]

clust=zz

png(file = "d15.png", res=300, width=2000, height=3000)
heatmap.2(as.matrix(clust), dendrogram="row",Rowv=TRUE,Colv=FALSE, scale="row",distfun=function(clust) dist(clust,method = 'euclidean'), hclustfun= function(clust) hclust(clust,method = 'ward'),colsep=c(1), col = bluered, trace="none", margins = c(12,16), labRow = FALSE)
dev.off()

#png(file = "d15_log.png", res=300, width=2000, height=2000)
#heatmap.2(as.matrix(log(clust+0.01)), dendrogram="row",Rowv=TRUE,Colv=FALSE, scale="row",distfun=function(clust) dist(log(clust+0.01),method = 'euclidean'), hclustfun= function(clust) hclust(log(clust+0.01),method = 'ward'),colsep=c(1), RowSideColors= sidecols, col = bluered, trace="none")
#dev.off()


