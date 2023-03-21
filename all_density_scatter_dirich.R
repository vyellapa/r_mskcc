library(ggplot2)
library(gridExtra)
library(RColorBrewer)
mycolors=colorRampPalette(brewer.pal(name="Paired", n = 12))(20)
mycolors=colorRampPalette(brewer.pal(name="Set1", n = 9))(9)
setwd("/Users/yellapav/Desktop/SOHN_163")

cgc=read.table("~/local/resources/cancer_gene_census_nov7.3_51_26.2016.tsv",sep="\t",header=T)

for(i in list.files(pattern=".*096.*_bestConsensusAssignments.bed",recursive=T)) {
#for(i in list.files(pattern=".*_bestConsensusAssignments.bed",recursive=T)) {
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
  
  
  
  T1_CAVE=grep("bkup",list.files(pattern=paste(name,".plus.strelka.rc.maf",sep=".*"), recursive=T),value=TRUE,invert = T)
  T2_CAVE=grep("bkup",list.files(pattern=paste(T2,".plus.strelka.rc.maf",sep=".*"), recursive=T),value=TRUE,invert = T)
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
  

  #merged=merged[merged$cluster_cluster %in% 1:20,]
  #096
  merged=merged[merged$cluster_cluster %in% c(1,2,3,5),]
  #097
  #merged=merged[merged$cluster_cluster %in% c(1,2,3,4,5,6,10),]
  #098
  #merged=merged[merged$cluster_cluster %in% c(1,2,3,4,5,6),]
  #099
  #merged=merged[merged$cluster_cluster %in% c(1,2,3,4,5,6),]
  #100
  #merged=merged[merged$cluster_cluster %in% c(1,2,3,4,5),]
  #101
  #merged=merged[merged$cluster_cluster %in% c(1,2,3,4),]
  merged=merge(merged, T1_CAVE,by.x="merger", by.y="merger", all.x = T )
  merged=merge(merged, T2_CAVE,by.x="merger", by.y="merger", all.x = T )
  #print(c(i,T2, T1_CCF, T2_CCF))
  see=merged
  see[is.na(see$T1_CAVE_Hugo_Symbol),"T1_CAVE_Hugo_Symbol"]=see[is.na(see$T1_CAVE_Hugo_Symbol),"T2_CAVE_Hugo_Symbol"]
  
  
  merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T1_CAVE_Hugo_Symbol=merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T2_CAVE_Hugo_Symbol
  merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T1_CAVE_Hugo_Symbol=merged[is.na(merged$T1_CAVE_Hugo_Symbol),]$T2_CAVE_Hugo_Symbol
  #merged=merged[,c(1,5,21,37,39,47,152,153,267,268)]
  merged=merged[,c(1,5,21,37,39,47,75,163,171,199,268)]
  
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
    theme(legend.position=c(1,1),legend.justification=c(1,1)) # +
    #xlab(name)+ylab(T2)+annotate("text", x = 0, y = 0.76, label = "",fontface="italic")+annotate("text", x = 0.00, y = 1.112, label = "JAK2\np.R683S",fontface="italic")+
    #annotate("text", x = 1.16, y = 0, label = "JAK2:p.R683G",fontface="italic")+annotate("text", x = 0.0, y = 0.626, label = "ZFHX3\n5'UTR",fontface="italic")+
    #annotate("text", x = 0, y = 0.96, label = "KRAS\np.A146T",fontface="italic")+annotate("text", x = 0.00, y = 0.78, label = "ACVR1\n5'UTR",fontface="italic")
  
  
  ################# Colors Not Enough #########################
  
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(plot.background = element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
          panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(),axis.ticks = element_blank() )
  
  dens_top <- ggplot(merged, aes(T1_subclonal.fraction,fill=factor(cluster_cluster)))+geom_density(alpha=0.4, linetype="blank")+
    theme_bw()+xlab("")+ylab("")+theme(legend.position="none")+scale_fill_manual(values=mycolors)+coord_cartesian(ylim=c(0,15), xlim=c(0,1.5))
  dens_right <- ggplot(merged, aes(T2_subclonal.fraction,fill=factor(cluster_cluster)))+geom_density(alpha=0.4, linetype="blank")+
    theme_bw()+xlab("")+ylab("")+theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+coord_flip()+scale_fill_manual(values=mycolors)+coord_flip(ylim=c(0,10),xlim=c(0,1.5))
  mid=ggplot(merged, aes(T1_subclonal.fraction, T2_subclonal.fraction, color=factor(cluster_cluster)))+geom_point(alpha=0.3,size=4)+coord_cartesian(xlim=c(0,1.5), ylim=c(0,1.5))+
    theme_bw()+ 
    scale_color_manual(values=mycolors)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(), strip.background  = element_blank(), 
          axis.ticks = element_line(size = 0.5), panel.grid.major = element_line(colour = "grey80"))+
    theme(legend.position=c(1,1),legend.justification=c(1,1))+
  
  #096
  xlab(name)+ylab(T2)+annotate("text", x = 0.36, y = 0, label = "JAK2\np.L611S",fontface="italic",size=4)+annotate("text", x = 0, y = 0.47, label = "JAK2\np.R683T",fontface="italic")
  #097
  #xlab(name)+ylab(T2)+annotate("text", x = 0.42, y = 0.87, label = "JAK2:T875N",fontface="italic",size=4)+annotate("text", x = 0.62, y = 1.02, label = "IKZF1:S41*",fontface="italic")+annotate("text", x = 0, y = 1.01, label = "KRAS\np.G13D",fontface="italic")+annotate("text", x = 0.0, y = 0.60, label = "TP53\np.R248Q",fontface="italic")
  #098
  #xlab(name)+ylab(T2)+annotate("text", x = 0.56, y = 0.0, label = "KRAS:Q61H",fontface="italic",size=4)
  #099
  #xlab(name)+ylab(T2)+annotate("text", x = 0, y = 0.76, label = "",fontface="italic")+annotate("text", x = 0.00, y = 1.11, label = "JAK2\np.R683S",fontface="italic")+annotate("text", x = 1.16, y = 0, label = "JAK2:p.R683G",fontface="italic")+annotate("text", x = 0.05, y = 0.96, label = "KRAS:p.A146T",fontface="italic",size=3)+annotate("text", x = 0.03, y = 0.90, label = "FOXO1\np.S22W",fontface="italic",size=3)
  #100
  #xlab(name)+ylab(T2)+annotate("text", x = 0.88, y = 0.95, label = "KMT2D\np.X5446_splice",fontface="italic",size=4)+annotate("text", x = 0.45, y = 0.81, label = "LZTR1\np.G405D",fontface="italic") 
  #101
  #xlab(name)+ylab(T2)+annotate("text", x = 0.41, y = 0.15, label = "FLT3\np.D835Y",fontface="italic",size=4)
  ###########################################################
  
  see=merged
  see[is.na(see$T1_CAVE_Hugo_Symbol),"T1_CAVE_Hugo_Symbol"]=see[is.na(see$T1_CAVE_Hugo_Symbol),"T2_CAVE_Hugo_Symbol"]
  see=see[see$T1_CAVE_Hugo_Symbol %in% cgc$Gene.Symbol,]
  
  
  
  grid.arrange(dens_top, empty, mid, dens_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  ggsave(s,file=paste(name,".pdf",sep = ""), width = 12, height = 12)
  
 # hdng=merged[,c(2,3,4)]
  #colnames(hdng)=c("cluster","prim.vaf","relapse.vaf")
  
  
}


merged=merged[merged$T1_CAVE_Hugo_Symbol %in% cgc$Gene.Symbol,]
merged=merged[merged$T1_CAVE_Variant_Classification!="Intron",]
mmm=merge(merged,T1_CAVE,by.x="merger", by.y="merger",all.x=T)
mmm=mmm[,c(1:6,47)]
#View(mmm)
