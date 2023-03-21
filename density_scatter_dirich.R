library(ggplot2)
library(gridExtra)
setwd("/Users/yellapav/Desktop/SOHN_163")

cgc=read.table("~/local/resources/cancer_gene_census_nov7.3_51_26.2016.tsv",sep="\t",header=T)

for(i in list.files(pattern=".*99.*_bestConsensusAssignments.bed",recursive=T)) {
  name=gsub("_bestConsensusAssignments.bed", "",strsplit(i,"/")[[1]][3])
  T2=gsub("-T1-", "-T2-", name)
  
  T1_CCF=list.files(pattern=paste(name,".*_allDirichletProcessInfo.txt",sep=""), recursive=T)
  T2_CCF=list.files(pattern=paste(T2,".*_allDirichletProcessInfo.txt",sep=""), recursive=T)
  
  
  CLUSTER=read.table(i, sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(CLUSTER)=paste("cluster",colnames(CLUSTER),sep="_")
  T1_CCF=read.table(T1_CCF, sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(T1_CCF)=paste("T1",colnames(T1_CCF),sep="_")
  T2_CCF=read.table(T2_CCF, sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(T2_CCF)=paste("T2",colnames(T2_CCF),sep="_")
  
  
  
  T1_CAVE=list.files(pattern=paste(name,".filter.caveman.rc.gm.maf",sep=".*"), recursive=T)
  T2_CAVE=list.files(pattern=paste(T2,".filter.caveman.rc.gm.maf",sep=".*"), recursive=T)
  T1_CAVE=read.table(T1_CAVE, sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
  colnames(T1_CAVE)=paste("T1_CAVE",colnames(T1_CAVE),sep="_")
  T2_CAVE=read.table(T2_CAVE, sep="\t", header=T, stringsAsFactors = FALSE, quote="~")
  colnames(T2_CAVE)=paste("T2_CAVE",colnames(T2_CAVE),sep="_")
  
  
  T1_CCF$merger=paste(T1_CCF$T1_chr, T1_CCF$T1_end,sep=":")
  T2_CCF$merger=paste(T2_CCF$T2_chr, T2_CCF$T2_end,sep=":")
  CLUSTER$merger=paste(CLUSTER$cluster_chr, CLUSTER$cluster_end, sep=":")
  T1_CAVE$merger=paste(T1_CAVE$T1_CAVE_Chromosome, T1_CAVE$T1_CAVE_Start_Position,sep=":")
  T2_CAVE$merger=paste(T2_CAVE$T2_CAVE_Chromosome, T2_CAVE$T2_CAVE_Start_Position,sep=":")
  
  
  merged=merge(CLUSTER, T1_CCF,by.x="merger", by.y="merger", all.x = T )
  merged=merge(merged, T2_CCF,by.x="merger", by.y="merger", all.x = T )
  
  merged=merged[merged$cluster_cluster %in% 1:7,]
  merged=merge(merged, T1_CAVE,by.x="merger", by.y="merger", all.x = T )
  merged=merge(merged, T2_CAVE,by.x="merger", by.y="merger", all.x = T )
  #print(c(i,T2, T1_CCF, T2_CCF))
  
  merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T1_CAVE_Hugo_Symbol=merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T2_CAVE_Hugo_Symbol
  merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T1_CAVE_Hugo_Symbol=merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T2_CAVE_Hugo_Symbol
  merged=merged[,c(1,5,21,37,39,47,152,153,267,268)]
  
  if(nrow(merged[is.na(merged$T1_subclonal.fraction),])>0) {merged[is.na(merged$T1_subclonal.fraction),]$T1_subclonal.fraction=0.00}
  if(nrow(merged[is.na(merged$T2_subclonal.fraction),])>0) {merged[is.na(merged$T2_subclonal.fraction),]$T2_subclonal.fraction=0.00}
  
  
  cgc_merged=merged[merged$T1_CAVE_Hugo_Symbol %in% cgc$Gene.Symbol, ]
  
  
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
       theme(plot.background = element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
 panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(),
 axis.text.y = element_blank(),axis.ticks = element_blank() )
  
  dens_top <- ggplot(merged, aes(T1_subclonal.fraction,fill=factor(cluster_cluster)))+geom_density(alpha=0.4, linetype="blank")+
    theme_bw()+xlab("")+ylab("")+theme(legend.position="none")+scale_fill_brewer(palette = "Set1")+coord_cartesian(ylim=c(0,15), xlim=c(0,1.5))
  dens_right <- ggplot(merged, aes(T2_subclonal.fraction,fill=factor(cluster_cluster)))+geom_density(alpha=0.4, linetype="blank")+
    theme_bw()+xlab("")+ylab("")+theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+coord_flip()+scale_fill_brewer(palette = "Set1")+coord_flip(ylim=c(0,10),xlim=c(0,1.5))
  mid=ggplot(merged, aes(T1_subclonal.fraction, T2_subclonal.fraction, color=factor(cluster_cluster)))+geom_point(alpha=0.2,size=4)+coord_cartesian(xlim=c(0,1.5), ylim=c(0,1.5))+
    theme_bw()+ 
    scale_color_brewer(palette = "Set1")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(), strip.background  = element_blank(), 
     axis.ticks = element_line(size = 0.5), panel.grid.major = element_line(colour = "grey80"))+
    theme(legend.position=c(1,1),legend.justification=c(1,1)) +
    xlab(name)+ylab(T2)+annotate("text", x = 0, y = 0.76, label = "",fontface="italic")+annotate("text", x = 0.00, y = 1.112, label = "JAK2\np.R683S",fontface="italic")+
    annotate("text", x = 1.16, y = 0, label = "JAK2:p.R683G",fontface="italic")+annotate("text", x = 0.0, y = 0.626, label = "ZFHX3\n5'UTR",fontface="italic")+
    annotate("text", x = 0, y = 0.96, label = "KRAS\np.A146T",fontface="italic")+annotate("text", x = 0.00, y = 0.78, label = "ACVR1\n5'UTR",fontface="italic")
  
  grid.arrange(dens_top, empty, mid, dens_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  
  
  hdng=merged[,c(2,3,4)]
  colnames(hdng)=c("cluster","prim.vaf","relapse.vaf")
}



x <- infer.clonal.models(variants=hdng,
                         cluster.col.name="cluster",
                         vaf.col.names=c("prim.vaf","relapse.vaf"),
                         subclonal.test="bootstrap",
                         subclonal.test.model="non-parametric",
                         cluster.center="mean",
                         num.boots=1000,
                         founding.cluster=1,
                         min.cluster.vaf=0.01,
                         p.value.cutoff=0.01,
                         alpha=0.1,
                         random.seed=63108)

f = generateFishplotInputs(results=x)
samples=c("T1","T2")
fishes = createFishPlotObjects(f)
#plot with fishplot
pdf('fish.pdf', width=8, height=5)
for (i in 1:length(fishes)){
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
           vlines=seq(1, length(samples)), vlab=samples, pad.left=0.5)
}
dev <- dev.off()