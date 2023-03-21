library(cowplot)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(gplots)
library(gtools)
library(dplyr)
mycolors1=colorRampPalette(brewer.pal(name="Spectral", n = 11))(14)
mycolors2=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000') 
mycolors3=colorRampPalette(brewer.pal(name="Set1", n = 9))(23)
mycol4=colorRampPalette(c("#a800ff","#743699","#1d0255","#0079ff","#2c61e3","#022394","#00f11d","#00a228","#ffef00","#fdd60f","#FABE28","#ff7f00","#EC8711","#ff0900","#8B4513"))(23)
mycol5=colorRampPalette(c("#BF2222","#EC8711","#F9FF00","#35A845","#3361E1"))(24)
mycol6=colorRampPalette(c("#BF2222","#EC8711","#F9FF00","#35A845","#3361E1"))(8)
mycolors=viridis_pal(option = "D")(12)
mycolors=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#8B4513","#a800ff","#19DD89","#1d0255","#0079ff","#2c61e3","#022394","#00f11d")

#cgc=read.table("~/local/resources/cancer_gene_census_nov7.3_51_26.2016.tsv",sep="\t",header=T)


setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/diagnostics/I-H-130718_step2/bkup/")
#setwd("/ifs/res/leukgen/projects/220/RESULTS/misc/dirichlet_triple/I-H-130720_step2/DP_output_merged/")

for(i in list.files(pattern=".*_allClusterassignmentsFromParallelRuns.txt",recursive=T)) {
  name=gsub("_allClusterassignmentsFromParallelRuns.txt", "",i)
  print(c(i,name))
  CLUSTER=read.table(i, sep="\t", header=T, stringsAsFactors = FALSE)
  CLUSTER=unique(CLUSTER)
  CLUSTER_NUM = (CLUSTER) %>% group_by(cluster) %>% summarize(number=n()) %>% right_join(CLUSTER, by="cluster") 
  CLUSTER_NUM$key = paste(CLUSTER_NUM$chr,CLUSTER_NUM$end,sep=":")
  setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/diagnostics")
  
  files = grep(name,list.files(pattern=".*_allDirichletProcessInfo.txt",recursive=T), value = T)
  
  #For boxplot 
  all = CLUSTER_NUM
  #all=read.table(files[1], sep="\t", header=T, stringsAsFactors = FALSE)
  #all$key=paste(all$chr, all$end,sep=":")
  
  for(j in grep(name,list.files(pattern=".*_allDirichletProcessInfo.txt",recursive=T), value = T)) {
    
    T1=read.table(j, sep="\t", header=T, stringsAsFactors = FALSE)
    T1$merger=paste(T1$chr, T1$end,sep=":")
    jname=gsub("_allDirichletProcessInfo.txt","",basename(j))
    jname=gsub('-','.',jname)
    colnames(T1)=paste(jname,colnames(T1),sep="_")
    all=merge(all, T1,by.x="key", by.y=paste(jname, "merger", sep="_"), all.x = T )
    all = all[,colnames(all) %in% c("key","cluster",grep("_subclonal.fraction$",colnames(all),value=T))]
    
    for(k in grep(j,files, value = T, invert = T)) {
      
      #setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/diagnostics")
      print("Reading ....")
      print(c(j,k))
      
      
      
      
      T2=read.table(k, sep="\t", header=T, stringsAsFactors = FALSE)
      T2$merger=paste(T2$chr, T2$end,sep=":")
      kname=gsub("_allDirichletProcessInfo.txt","",basename(k))
      kname=gsub('-','.',kname)
      colnames(T2)=paste(kname,colnames(T2),sep="_")
      
      print("here..")
      merged=merge(CLUSTER_NUM, T1,by.x="key", by.y=paste(jname, "merger", sep="_"), all.x = T )
      merged=merge(merged, T2,by.x="key", by.y=paste(kname, "merger", sep="_"), all.x = T )
      
      merged=merged %>% dplyr::filter(number>50)
      
      
      rm_file=paste(kname,"_vs_",jname,"pdf",sep = '.')
      
      if(file.exists(rm_file)) {
        print(c("Skipping ",rm_file))
        break
        
      }
      
      #Plotting
        filler=grep("_cluster", colnames(merged),value=T)
      
        empty <- ggplot()+geom_point(aes(1,1), colour="white") +
        theme(plot.background = element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
            panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(),axis.ticks = element_blank() )
      
     
      
        dens_top <- ggplot(merged, aes_string(paste(jname,"subclonal.fraction", sep="_")))+geom_density(alpha=0.4, linetype="blank", aes(fill=factor(cluster)))+
          theme_bw()+xlab("")+ylab("")+theme(legend.position="none", axis.text.x = element_blank() )+scale_fill_manual(values = mycolors)+coord_cartesian(ylim=c(0,10), xlim=c(0,2))
      
        dens_right <- ggplot(merged, aes_string(paste(kname,"subclonal.fraction", sep="_")))+geom_density(alpha=0.4, linetype="blank", aes(fill=factor(cluster)))+
        theme_bw()+xlab("")+ylab("")+theme(legend.position="none", axis.text.y = element_blank())+coord_flip()+scale_fill_manual(values = mycolors)+coord_flip(ylim=c(0,10),xlim=c(0,2))
      
      
      #mid <- ggplot(merged, aes_string(paste(jname,"subclonal.fraction", sep="_"), paste(kname,"subclonal.fraction", sep="_"), fill=factor(filler)))+geom_point(alpha=0.35,size=3, aes(color=factor(filler)))+coord_cartesian(xlim=c(0,2), ylim=c(0,2))+
       mid <- ggplot(merged, aes_string(paste(jname,"subclonal.fraction", sep="_"), paste(kname,"subclonal.fraction", sep="_")))+geom_point(alpha=0.5,size=1.25, aes(color=factor(cluster)))+coord_cartesian(xlim=c(0,2), ylim=c(0,2))+
         theme_bw()+xlab(jname)+ylab(kname)+
         scale_color_manual(values = mycolors)+
         theme(axis.text.x = element_text(angle = 0, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 1) , legend.position="bottom",legend.direction="horizontal",legend.title=element_blank(), strip.background  = element_blank(), 
                      axis.ticks = element_line(size = 0.5), axis.title=element_text(size=12,face="bold"), axis.text=element_text(size=12), panel.grid.major = element_line(colour = "grey80"))+
          theme(legend.position=c(1,1),legend.justification=c(1,1))
      
      
        s=grid.arrange(dens_top, empty, mid, dens_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
        s=plot_grid(dens_top, NULL, mid,dens_right, align = "hv", nrow = 2, rel_heights = c(0.2, 0.8, 0.8), rel_widths = c(0.8,0.2))
      #setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/diagnostics/plots")
        ggsave(s,file=paste(jname,"_vs_",kname,"pdf",sep = '.'), width = 12, height = 12)
      
        rm_file=paste(kname,"_vs_",jname,"pdf",sep = '.')
        unlink(rm_file)
      
      
      pp=ggplot(merged, aes(as.character(cluster),fill=factor(chr)))+geom_bar(stat="count")
      qq=ggplot(merged, aes(as.character(cluster),fill=factor(chr)))+geom_bar(stat="count", position = "fill")
      
      png(file = paste(jname,"bar.png",sep="_"), res=300, width=3200, height=1800)
        gridExtra::grid.arrange(pp,qq,ncol=2)
      dev.off()
      
    }
  }
  
}


alll=all[,colnames(all) %in% c("key","cluster",grep("_subclonal.fraction$",colnames(all),value=T))]
alll = alll %>% dplyr::select(-key) %>% melt(id.vars="cluster") 

alll$variable=gsub(".D1.2|_subclonal.fraction","",alll$variable)
setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/diagnostics/plots")
s=ggplot(alll, aes(variable,value,fill=variable))+geom_boxplot()+facet_wrap(~cluster,nrow=1)+theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")+scale_fill_manual(values = mycolors)+xlab("")+ylab("CCF")
ggsave(s,file=paste(jname,"boxplot","pdf",sep = '.'), width = 18, height = 8)

s=ggplot(alll, aes(variable,value,fill=variable))+geom_violin(scale="width")+facet_wrap(~cluster,nrow=1)+theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")+scale_fill_manual(values = mycolors)+xlab("")+ylab("CCF")
ggsave(s,file=paste(jname,"violin","pdf",sep = '.'), width = 18, height = 8)

all2 = all %>% dplyr::filter(cluster==2) %>% unique() %>% distinct(key, .keep_all = TRUE)
colnames(all2)=gsub(".D1.2|_subclonal.fraction|.H.|.D1.1|T1.","",colnames(all2))


row.names(all2)=all2$key
all2=all2 %>% dplyr::select(-cluster,-key)

tempData=all2
heatmap.2(as.matrix(tempData), dendrogram="both",scale="row",Rowv=TRUE,Colv=TRUE, distfun=function(tempData) dist(tempData,method = 'euclidean'), hclustfun= function(tempData) hclust(tempData,method = 'ward'),colsep=c(1), col = bluered, trace="none")

tempData = tempData %>% dplyr::filter(I130718.1<0.4 & I130718.11<0.4 & I130718.12<0.4 & I130718.2<0.4 & I130718.4<0.4 & I130718.6<0.4 & I130718.9<0.4)
par(mar=c(7,4,4,2)+0.1) 
heatmap.2(as.matrix(tempData), dendrogram="both",Rowv=TRUE,Colv=TRUE, distfun=function(tempData) dist(tempData,method = 'euclidean'), hclustfun= function(tempData) hclust(tempData,method = 'ward'),colsep=c(1), col = bluered, trace="none",margins=c(6.6,1.5))



