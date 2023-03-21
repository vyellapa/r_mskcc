library(ape)
library(phylobase)
library(reshape2)
library(ggplot2)



a = "((PD26403.3:282, (PD26403.2:1248,PD26403.10:1) PD26403.4:173, PD26403.1:3105)PD26403.7:5678) I:3;"

####################################
#I-H-106917
stri="I-H-106917"
a="((((I6917.8:186,I6917.1001:3)I6917.10:345,(I6917.5:283,I6917.100:11)I6917.6:152)I6917.11:128,(I6917.3:361,(I6917.1:221,I6917.102:1)I6917.2:275)I6917.4:243)I6917.13:9085)I:1;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/p220_2019/dirichlet/c2/%s_tree.pdf",stri), useDingbats=FALSE,width=14,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()


####################################
#I-H-130719
stri="I-H-130719"
a="(((((I719.102:1,I719.7:94)I719.10:269,I719.101:1)I719.11:185,((I719.2:931,I719.100:1)I719.3:393,I719.1:665)I719.4:63)I719.12:666,(I719.5:224,I719.103:1)I719.6:138)I719.14:8134)I:1;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/p220_2019/dirichlet/c2/%s_tree.pdf",stri), useDingbats=FALSE,width=14,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()


####################################
#I-H-130720
stri="I-H-130720"
a="(((I720.10:284,I720.5:740,(I720.1:2967,I720.101:1)I720.2:611,I720.3:1013)I720.11:143,(I720.7:697,I720.100:1)I720.8:554)I720.13:3739)I1:2000;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/p220_2019/dirichlet/c2/%s_tree.pdf",stri), useDingbats=FALSE,width=14,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()


####################################
#I-H-130718
stri="I-H-130718"
a="(((I718.4:1779,I718.7:1993)I718.11:474,((I718.13:1622,I718.6:2760)I718.10:1564,(I718.1:3156,(I718.5:2792,I718.8:2484,I718.9:2149)I718.12:284)I718.16:125)I718.15:65)I718.3:5891)I1:1;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/p220_2019/dirichlet/c2/%s_tree.pdf",stri), useDingbats=FALSE,width=14,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()





a="(((I718.4:1779,I718.7:1993)I718.11:474,((I718.13:1622,I718.6:2760)I718.10:1564,(I718.1:3156,I718.5:2792,I718.8:2484,I718.9:2149)I718.16:125)I718.15:65)I718.3:5891)I1:1;"


tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/p220_2019/dirichlet/c2/%s_tree_rm12.pdf",stri), useDingbats=FALSE,width=14,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()



##### Signatures ############################

sig=read.table("~/Downloads/DP_clusters_mmsig_manual_for_teja.txt",header=T,sep="\t")
sig=sig[sig$X!="remove",]
col=c(RColorBrewer::brewer.pal(8, "Dark2")[c(8)], RColorBrewer::brewer.pal(8, "Paired")) 
col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Paired"))
melter <- function(stri) {
  
  s=sig[grep(stri,sig$sample_Id),]
  s$SBS1=s$SBS1*s$mutations
  s$SBS2=s$SBS2*s$mutations
  s$SBS5=s$SBS5*s$mutations
  s$SBS8=s$SBS8*s$mutations
  s$SBS9=s$SBS9*s$mutations
  s$SBS13=s$SBS13*s$mutations
  s$SBS18=s$SBS18*s$mutations
  s$SBS.MM1=s$SBS_MM1*s$mutations
  s$SBS35=s$SBS35*s$mutations
  s$sample_Id=gsub(" ","_",s$sample_Id)
  s$sample_Id=gsub("-H-|106917|130917|130720|130718","",s$sample_Id)
  s=s[,c("sample_Id","SBS1","SBS2","SBS5","SBS8","SBS9","SBS13","SBS18","SBS.MM1","SBS35")]
  
  m=melt(s)
  head(m)
  pdf(sprintf("~/Desktop/%s_bargg.pdf",stri), useDingbats=FALSE)
  p = ggplot(m,aes(sample_Id,value,fill=variable))+geom_bar(stat="identity")+theme_minimal()+scale_fill_manual(values = col)
  print(p)
  #ggsave(p, file=sprintf("~/Desktop/%s_bargg.pdf",stri))
  dev.off()
  
  return(m)
}

m = melter("I-H-106917")
m = melter("I-H-130718")
m = melter("I-H-130719")
m = melter("I-H-130720")





##### CNV Changes ##############
library(GenomicRanges)
library(dplyr)
seg = read.table("~/Desktop/p220_2019/genome_plots/input/ranges.txt",header=F,sep="\t")
cyto = read.table("~/Desktop/p220_2019/genome_plots/input/Myeloma_SeqFISH_Probe_Locations.txt", header = F, sep = "\t")
head(seg)
head(cyto)


gr.cyto = GRanges(seqnames=Rle(cyto$V1), IRanges(cyto$V2, cyto$V3), arm=cyto$V4)
gr.seg = GRanges(seqnames=Rle(seg$V2), IRanges(seg$V3, seg$V4))

overlapGenes <- findOverlaps(gr.cyto, gr.seg, minoverlap=100)
df = data.frame(seg[subjectHits(overlapGenes),], cyto[queryHits(overlapGenes),])

colnames(df) = c("sample","chr","start","stop","ntot","cell.frac","maj","min","clonal","chrom","st","sp","aber","known")
df$patient=substr((df$sample),1,10)

df = df %>% dplyr::filter(!is.na(maj)) %>%
  dplyr::filter(patient!="I-H-106917AAAA") %>%
  dplyr::mutate(key=paste(patient,clonal,aber,maj,min,sep="_")) %>%
  unique() %>% dplyr::mutate(mb=(stop-start)/1000000)
rec = df %>% 
  dplyr::group_by(key) %>%
  summarise(recurrance = n()) 

df = df %>% left_join(rec,by="key") %>% dplyr::filter(!(maj==1 & min==1))

head(df)
head(rec)

write.table(df,file="~/Desktop/p220_2019/genome_plots/input/recurrance.txt", append=FALSE, sep="\t", eol="\n", row.names=F, col.names=TRUE, quote=FALSE)





##################################################
######        Create bins              ###########
##################################################
cyto.bed=read.table("/Users/yellapav/SV_project/data/cytoband_b37.bed",sep="\t")
colnames(cyto.bed)=c("chr","start","stop","band","direction")
cyto <- data.frame(chrom=character(), 
                   start=numeric(), 
                   stop=numeric(),
                   arm=character(), 
                   stringsAsFactors=FALSE)


#Create cytoband df with min,max for each arm
cyto.bed$arm=paste(cyto.bed$chr,substr(cyto.bed$band,0,1),sep="")
for(i in unique(cyto.bed$arm)) {subset=cyto.bed[cyto.bed$arm==i,]; 
min=min(subset$start)
max=max(subset$stop)
cyto <- rbind(cyto, data.frame(chrom=unique(subset$chr), start=min, stop=max,arm=i))
}

bin.size=10000
#Create dataframe with bins of bin.size into df called bins.df
bins.df <- data.frame(chrom=character(), 
                      start=numeric(), 
                      stop=numeric(),
                      mid=numeric(),
                      bin.num=character(), 
                      stringsAsFactors=FALSE)
bin=0
sa=0

for(i in unique(cyto$arm)){
  sub=cyto[cyto$arm==i,]
  chrom=sub$chrom
  start=sub$start
  stop=sub$stop
  
  
  starts=seq(start, stop, by=bin.size)
  stops=starts+bin.size
  chroms=rep(as.character(chrom),length(starts))
  arms=rep(as.character(i),length(starts))
  bins.df=rbind(bins.df, data.frame(chrom=chroms, start=starts, stop=stops, mid=starts,bin.num=arms))
}

bins.df$bin.num=paste(bins.df$bin.num,bins.df$start,bins.df$stop,sep="_")


trials = seg
colnames(trials) = c("sample","chr","start","stop","ntot","cell.frac","maj","min","clonal")
trials[is.na(trials)] <- 0

gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop), bin.num=bins.df$bin.num)
gr.sv = GRanges(seqnames=Rle(as.character(trials$chr)), IRanges(trials$start, trials$stop))

overlapGenes <- findOverlaps(gr.bins, gr.sv)
df.plot = data.frame(trials[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])

df.plot1 = df.plot %>% dplyr::filter(sample=="I-H-106917-T2-1-D1-2" | sample=="I-H-106917-T2-2-D1-2" | sample=="I-H-106917-T2-3-D1-2"| sample=="I-H-106917-T2-4-D1-2")
df.plot1$order = factor(df.plot1$chrom, levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))
df.plot1 = (df.plot1) %>% dplyr::mutate(value=(min+maj)*cell.frac) #%>% dplyr::mutate(value = ifelse(clonal=="subclonal",value*-1,value))

ggplot(df.plot1, aes(start.1,(value),color=clonal))+geom_line()+facet_grid(sample~order,scales = "free_x",space="free")+
  theme_bw(base_size = 22)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position="right", strip.background  = element_blank(), 
        legend.title = element_blank(),
        axis.ticks = element_line(size = 0), 
        panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.5, colour = "#666666", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_color_manual(values=c("dodgerblue","brown2"))+
  scale_fill_manual(values=c("#007B22"))+
  theme(legend.title=element_blank(),plot.margin = unit(c(0,0,0,0.1), "lines"))+
  xlab("")+ylab("Frequency")#+
  #coord_cartesian(ylim=c(-55,78)) 














