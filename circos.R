library(RCircos)
library(circlize)
library(dplyr)

sample="I-H-130720"
somp="I-H-130720-T1"

#sample="I-H-106917"
#somp="I-H-106917-T2"
text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)

setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data")
s=read.table(sprintf("%s.trunk",somp),sep="\t") %>% arrange(V1,V2)
t1=as.vector(s$V2)
s$row=c(t1[2:length(t1)],0)
s$idist=s$row-s$V2
s = s %>% dplyr::filter(idist>0) %>% dplyr::mutate(key=paste(V1,V2,sep=":"))
s$idist = log10(s$idist)
see=read.table(sprintf("%s_CA.txt",somp),sep="\t",header = F)
b1=read.table(sprintf("%s_CA.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b2=read.table(sprintf("%s_CG.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b3=read.table(sprintf("%s_CT.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b4=read.table(sprintf("%s_TA.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b5=read.table(sprintf("%s_TC.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b6=read.table(sprintf("%s_TG.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)


setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/bbg")
bbg=read.table(sprintf("%s_shared.txt",sample),header = T, sep = "\t") %>% 
  dplyr::select(chr,startpos, endpos, ntot) %>% dplyr::rename(value1 = ntot) %>% dplyr::mutate(chr=paste0("chr",chr))
bbg.max = max(bbg$value1)
bbg.min = min(bbg$value1)


sv = read.table("~/Downloads/all_brass.txt",header=T,sep="\t") %>% 
  dplyr::select(-c("contributing_read_names","contributing_read_number")) %>% 
  mutate(PCAWG_class=ifelse(PCAWG_class=="unbalanced_translocation","Unbalanced_translocation",as.character(PCAWG_class))) %>%
  mutate(PCAWG_class=ifelse(PCAWG_class=="Chromothripsis","chromothripsis",as.character(PCAWG_class))) %>%
  filter(PCAWG_class!="CSR" & shared!="ART")

sv=sv[grep(sample,sv$sample),] %>% filter(shared=="SHARED")

ig.tra = sv

ig.tra = ig.tra %>% mutate(color.ig = ifelse(PCAWG_class=="reciprocal_translocation", "#FB8072",
                                       ifelse(PCAWG_class=="chromoplexy", "#add8e6",
                                       ifelse(PCAWG_class=="chromothripsis", "#984EA3",
                                       ifelse(PCAWG_class=="complex", "#FF7F00",
                                       ifelse(PCAWG_class=="unclassified translocation", "black",
                                       ifelse(PCAWG_class=="local_n_jumps", "#377EB8",
                                       ifelse(PCAWG_class=="templated_insertion", "#A65628",
                                       ifelse(PCAWG_class=="Unbalanced_translocation", "#0a5d00",
                                       ifelse(PCAWG_class=="single", "#377EB8","#A9A9A9"))))))))))


sv1=ig.tra[,c("chr1.x","start1.x","end1.x","color.ig")]
sv1$chr1.x=paste0("chr",sv1$chr1.x)


sv2=ig.tra[,c("chr2.x","start2.x","end2.x","color.ig")]
sv2$chr2.x=paste0("chr",sv2$chr2.x)


colnames(sv1)=c("chr","start","end","value1")
colnames(sv2)=c("chr","start","end","value1")


pdf(sprintf("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/plots/%s_circos_shared.pdf",sample), useDingbats=FALSE,width = 7, height = 7)

circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
#circos.genomicLabels(text, labels.column = 4, cex=.6,side = "outside")

#circos.initializeWithIdeogram(plotType = NULL,ideogram.height=convert_height(3,"mm"))
circos.initializeWithIdeogram(plotType = c("ideogram", "axis", "labels"),ideogram.height=convert_height(1.75,"mm"),track.height =convert_height(1,"mm") )

bed_list = list(b1, b2,b3, b4,b5,b6)
col=c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C","#4DAF4A","#E41A1C")
#palette(c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C"))




circos.genomicTrack(bed_list, track.height = 0.125,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = i, ...)
                    })


#circos.genomicRainfall(bed_list,pch = 16, cex = 0.4, ylim=c(0,6.8),col =c("#add8e6","#984EA3","#FF7F00","#A65628","#FF7F00","#A65628"),track.height = 0.1)
circos.genomicRainfall(bed_list,pch = 16, cex = 0.5,ylim=c(0,7.75),col =c("#0000FF","#000000","#FF0000","#808080","#00ff00","#ffa8b7"),track.height = 0.2)

bed1 = generateRandomBed(nr = 100)
#z1 = read.table("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/I-H-106917_shared.dels.txt",header=F,sep="\t")
#z2 = read.table("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/I-H-106917_shared.ins.txt",header=F,sep="\t")
setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/indels")
z1 = read.table(list.files(pattern=sprintf("%s.*shared.del.txt$",sample)),header=F,sep="\t")
z2 = read.table(list.files(pattern=sprintf("%s.*shared.ins.txt$",sample)),header=F,sep="\t")

colnames(z1) = colnames(bed1)
colnames(z2)= colnames(bed1)
bed_list = list(z1,z2)
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("forestgreen", "black", "red"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                # i = getI(...)
                                #cell.xlim = get.cell.meta.data("cell.xlim")
                                #circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              }, track.height=0.05)

bbg1 = bbg %>% filter(value1>2) %>% mutate(value1=ifelse(value1>2.5,4,2)) %>% filter(value1!=2)
bbg2 = bbg %>% filter(value1<=2) %>% mutate(value1=ifelse(value1<1.5,0.5,2)) %>% filter(value1!=2)

bed_list = list(bbg1,bbg2)
f = colorRamp2(breaks = c(bbg.min,(bbg.min+0.45), 2, (bbg.max-2.5),bbg.max), colors = c("firebrick1","firebrick1", "white", "dodgerblue", "dodgerblue"))
f = colorRamp2(breaks = c(bbg.min,1.6, 2, 2.4,bbg.max), colors = c("firebrick1","firebrick1", "white", "dodgerblue", "dodgerblue"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                # i = getI(...)
                                #cell.xlim = get.cell.meta.data("cell.xlim")
                                #circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              }, track.height=0.075)



circos.genomicLink(sv1, sv2,col = as.vector(sv1$value1))


circos.clear()
dev.off()








#######################
##### Unique ##########
library(RCircos)
library(circlize)
library(dplyr)


text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)

setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data")
s=read.table(sprintf("%s_uniqued.temp",somp),sep="\t") %>% arrange(V1,V2)
t1=as.vector(s$V2)
s$row=c(t1[2:length(t1)],0)
s$idist=s$row-s$V2
s = s %>% dplyr::filter(idist>0) %>% dplyr::mutate(key=paste(V1,V2,sep=":"))
s$idist = log10(s$idist)

setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/unique")
b1=read.table(sprintf("%s_CA.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b2=read.table(sprintf("%s_CG.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b3=read.table(sprintf("%s_CT.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b4=read.table(sprintf("%s_TA.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b5=read.table(sprintf("%s_TC.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)
b6=read.table(sprintf("%s_TG.txt",somp),sep="\t",header = F) %>% dplyr::mutate(key=paste(V1,V2,sep=":")) %>% 
  left_join(s,by="key") %>% dplyr::mutate(chr=paste0("chr",V1.x), start=V2.x,end=V2.x+1) %>% dplyr::select(chr,start,end,idist)

setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/bbg")
bbg=read.table(sprintf("%s_uniqued.txt",sample),header = T, sep = "\t") %>% 
  dplyr::select(chr,startpos, endpos, ntot) %>% dplyr::rename(value1 = ntot) %>% dplyr::mutate(chr=paste0("chr",chr))
bbg.max = max(bbg$value1)
bbg.min = min(bbg$value1)


sv = read.table("~/Downloads/all_brass.txt",header=T,sep="\t") %>% 
  dplyr::select(-c("contributing_read_names","contributing_read_number")) %>% 
  mutate(PCAWG_class=ifelse(PCAWG_class=="unbalanced_translocation","Unbalanced_translocation",as.character(PCAWG_class))) %>%
  mutate(PCAWG_class=ifelse(PCAWG_class=="Chromothripsis","chromothripsis",as.character(PCAWG_class))) %>%
  filter(PCAWG_class!="CSR" & shared!="ART")

sv=sv[grep(sample,sv$sample),] %>% filter(shared=="UNIQUE")

ig.tra = sv

ig.tra = ig.tra %>% mutate(color.ig = ifelse(PCAWG_class=="reciprocal_translocation", "#FB8072",
                                             ifelse(PCAWG_class=="chromoplexy", "#add8e6",
                                                    ifelse(PCAWG_class=="chromothripsis", "#984EA3",
                                                           ifelse(PCAWG_class=="complex", "#FF7F00",
                                                                  ifelse(PCAWG_class=="unclassified translocation", "black",
                                                                         ifelse(PCAWG_class=="local_n_jumps", "#377EB8",
                                                                                ifelse(PCAWG_class=="templated_insertion", "#A65628",
                                                                                       ifelse(PCAWG_class=="Unbalanced_translocation", "#0a5d00",
                                                                                              ifelse(PCAWG_class=="single", "#377EB8","#A9A9A9"))))))))))


sv1=ig.tra[,c("chr1.x","start1.x","end1.x","color.ig")]
sv1$chr1.x=paste0("chr",sv1$chr1.x)


sv2=ig.tra[,c("chr2.x","start2.x","end2.x","color.ig")]
sv2$chr2.x=paste0("chr",sv2$chr2.x)


colnames(sv1)=c("chr","start","end","value1")
colnames(sv2)=c("chr","start","end","value1")


pdf(sprintf("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/plots/%s_circos_unique.pdf",sample), useDingbats=FALSE,width = 7, height = 7)

circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
#circos.genomicLabels(text, labels.column = 4, cex=.6,side = "outside")

#circos.initializeWithIdeogram(plotType = NULL,ideogram.height=convert_height(3,"mm"))
circos.initializeWithIdeogram(plotType = c("ideogram", "axis", "labels"),ideogram.height=convert_height(1.75,"mm"),track.height =convert_height(1,"mm") )

bed_list = list(b1, b2,b3, b4,b5,b6)
col=c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C","#4DAF4A","#E41A1C")
#palette(c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C"))




circos.genomicTrack(bed_list, track.height = 0.125,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = i, ...)
                    })


#circos.genomicRainfall(bed_list,pch = 16, cex = 0.4, ylim=c(0,6.8),col =c("#add8e6","#984EA3","#FF7F00","#A65628","#FF7F00","#A65628"),track.height = 0.1)
circos.genomicRainfall(bed_list,pch = 16, cex = 0.5,ylim=c(0,7.75),col =c("#0000FF","#000000","#FF0000","#808080","#00ff00","#ffa8b7"),track.height = 0.2)

bed1 = generateRandomBed(nr = 100)
#z1 = read.table("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/I-H-106917_shared.dels.txt",header=F,sep="\t")
#z2 = read.table("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/I-H-106917_shared.ins.txt",header=F,sep="\t")
setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/indels")
z1 = read.table(list.files(pattern=sprintf("%s.*unique.del.txt$",sample)),header=F,sep="\t")
z2 = read.table(list.files(pattern=sprintf("%s.*unique.ins.txt$",sample)),header=F,sep="\t")

colnames(z1) = colnames(bed1)
colnames(z2)= colnames(bed1)
bed_list = list(z1,z2)
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("forestgreen", "black", "red"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                # i = getI(...)
                                #cell.xlim = get.cell.meta.data("cell.xlim")
                                #circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              }, track.height=0.05)

bbg1 = bbg %>% filter(value1>2) %>% mutate(value1=ifelse(value1>2.5,4,2)) %>% filter(value1!=2)
bbg2 = bbg %>% filter(value1<=2) %>% mutate(value1=ifelse(value1<1.5,0.5,2)) %>% filter(value1!=2)

bed_list = list(bbg1,bbg2)
f = colorRamp2(breaks = c(bbg.min,(bbg.min+0.45), 2, (bbg.max-2.5),bbg.max), colors = c("firebrick1","firebrick1", "white", "dodgerblue", "dodgerblue"))
f = colorRamp2(breaks = c(bbg.min,1.6, 2, 2.4,bbg.max), colors = c("firebrick1","firebrick1", "white", "dodgerblue", "dodgerblue"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                # i = getI(...)
                                #cell.xlim = get.cell.meta.data("cell.xlim")
                                #circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              }, track.height=0.075)



circos.genomicLink(sv1, sv2,col = as.vector(sv1$value1))


circos.clear()

dev.off()

library(BSgenome)
library(reshape2)
library(ggplot2)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/vcf")
vcf_files=c(sprintf("%s_s.vcf",somp),sprintf("%s_u.vcf",somp))
sample_names=c(sprintf("%s Trunk",sample),sprintf("%s Unique",sample))
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat)

m=as.data.frame(mut_mat)
#> View(m)
m=t(as.data.frame(mut_mat))
m=(as.data.frame(mut_mat))
m$id=row.names(m)
m=melt(m)
m$nt=substr(m$id,3,5)

pdf(sprintf("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/plots/%s_96nt.pdf",sample), useDingbats=FALSE,width = 7, height = 5)
ggplot(m,aes(id,value,fill=nt))+geom_bar(stat="identity")+
  facet_grid(variable~nt,scales = "free_x")+
theme_bw(base_size = 15)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), #axis.text.y = element_blank(),
        legend.position="none", 
        strip.background  = element_rect(colour="black", fill="white", size=0.1, linetype="solid"), 
        strip.text.x = element_text(size=12, color="black",face="plain"),
        strip.text.y = element_text(size=12, color="black",face="plain"),
        legend.title=element_blank(), 
        axis.ticks = element_line(size = 0), 
        #panel.grid.major = element_blank(), 
        panel.spacing = unit(0.0005, "lines"),
        panel.border = element_rect(size = 0.001, colour = "#000000", fill = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
        #strip.text.x = element_blank()
        )+ 
  #theme(strip.text.x = element_text(size=12, color="black",face="bold"),)+
  theme(legend.title=element_blank(), plot.margin = unit(c(0,0,0.25,2.61), "lines"))+
  scale_fill_manual(values=c("#0000FF","#000000","#FF0000","#808080","#00ff00","#ffa8b7"))+
  ylab("")+xlab("")#+coord_cartesian(ylim=c(0,1))
dev.off()