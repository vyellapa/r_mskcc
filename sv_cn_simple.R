#circ<- kat_final
#circ$chr<- paste("chr", kat_final$chrom, sep="")

library(RCircos)
library(circlize)
library(dplyr)
h=read.table("~/Desktop/try_circlise.tsv",sep="\t",header = T)
colnames(h)=c("chr","start","end","value1")
head(h)


text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)
a1=read.table("/Users/yellapav/Downloads/chromoplexy_hotspots.csv",sep=",",header = T)
a2=read.table("/Users/yellapav/Downloads/chromothripsis_hotspots.csv",sep=",",header = T)
a3=read.table("/Users/yellapav/Downloads/complex_hotspots.csv",sep=",",header = T)
a4=read.table("/Users/yellapav/Downloads/templated_insertion_hotspots.csv",sep=",",header = T)
b1=read.table("/Users/yellapav/Downloads/del_hotspots.csv",sep=",",header = T)
b2=read.table("/Users/yellapav/Downloads/dup_hotspots.csv",sep=",",header = T)
b3=read.table("/Users/yellapav/Downloads/inv_hotspots.csv",sep=",",header = T)
b4=read.table("/Users/yellapav/Downloads/tra_hotspots.csv",sep=",",header = T)




h = h %>% mutate(value1 = ifelse(value1 > 40, 40,value1))


circos.par("start.degree" = 90)
circos.genomicLabels(text, labels.column = 4, cex=.6,side = "outside")

#circos.initializeWithIdeogram(plotType = NULL,ideogram.height=convert_height(3,"mm"))
circos.initializeWithIdeogram(ideogram.height=convert_height(2,"mm"))


#bed = generateRandomBed(nr = 3, nc = 2)
bed1 = generateRandomBed(nr = 500)
bed2 = generateRandomBed(nr = 500)
bed3 = generateRandomBed(nr = 500)
bed4 = generateRandomBed(nr = 500)
bed_list = list(b1, b2,b3, b4)
col=c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C")
#palette(c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C"))


#circos.genomicLabels(text, labels.column = 4, cex=.4,side = "outside")
#circos.genomicIdeogram(species = "hg19")

circos.genomicRainfall(bed_list,pch = 16, cex = 0.4, ylim=c(0,6.8),col =c("#E41A1C","#377EB8","#4DAF4A","#000000"),track.height = 0.1)
#circos.genomicTrack(bed_list, track.height = 0.125,
#                    panel.fun = function(region, value, ...) {
#                      i = getI(...)
#print(i)
#                      colors(c("#000000","#377EB8","#4DAF4A","#E41A1C"))
#col=rep(c("#000000","#377EB8","#4DAF4A","#E41A1C"),8)
#col=c(rep("#000000",nrow(a1)),rep("#377EB8",nrow(a2)),rep("#4DAF4A",nrow(a3)),rep("#E41A1C",nrow(a4)))
#if(i==3) {i=137}
#if(i==2) {i=139}
#if(i==1) {i=135}
#if(i==4) {i=1}
#c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C")
#                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = i, ...)
#                    })

bed_list = list(a1, a2,a3, a4)

#circos.genomicTrack(bed_list, track.height = 0.125,
#                    panel.fun = function(region, value, ...) {
#                      i = getI(...)
#                      circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = i, ...)
#                    })


circos.genomicRainfall(bed_list,pch = 16, cex = 0.4, ylim=c(0,6.8),col =c("#add8e6","#984EA3","#FF7F00","#A65628"),track.height = 0.1)

circos.genomicTrack(h, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h")
                    })

circos.clear()












################################################################################################################
library(RCircos)
library(circlize)
library(dplyr)
h=read.table("~/Desktop/try_circlise.tsv",sep="\t",header = T)
colnames(h)=c("chr","start","end","value1")
head(h)


text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)
a1=read.table("/Users/yellapav/Downloads/Chromoplexy_hotspots.csv",sep=",",header = T)
a2=read.table("/Users/yellapav/Downloads/Chromothripsis_hotspots.csv",sep=",",header = T)
a3=read.table("/Users/yellapav/Downloads/Complex_hotspots.csv",sep=",",header = T)
a4=read.table("/Users/yellapav/Downloads/Templated_insertion_hotspots.csv",sep=",",header = T)
b1=read.table("/Users/yellapav/Downloads/DEL_hotspots.csv",sep=",",header = T)
b2=read.table("/Users/yellapav/Downloads/DUP_hotspots.csv",sep=",",header = T)
b3=read.table("/Users/yellapav/Downloads/INV_hotspots.csv",sep=",",header = T)
b4=read.table("/Users/yellapav/Downloads/TRA_hotspots.csv",sep=",",header = T)


h = h %>% mutate(value1 = ifelse(value1 > 40, 40,value1))


library(circlize)

par(lwd = 0.5)

circos.par("cell.padding" = c(0, 0, 0, 0),"start.degree" = 90)
#circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)
circos.initializeWithIdeogram(plotType = NULL)

posTransform.fun = function(region) {
  return(region)
}

circos.genomicTrackPlotRegion(xxx, ylim = c(6.9, 7.9), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 0.6,font=2, posTransform = posTransform.fun,niceFacing = TRUE, track.margin=c(0.01,0.01))
}, track.height = 0.05, bg.border = NA)

#circos.genomicTrackPlotRegion(xx, ylim = c(0, 1), panel.fun = function(region, value, ...) {
#    circos.genomicText(region, value, y = 0, labels.column = 4, facing = "clockwise", adj = c(0, 0.5), cex = 0.5, posTransform = posTransform.fun,niceFacing = TRUE)
#  }, track.height = 0.05, bg.border = NA)

#circos.genomicPosTransformLines(text, posTransform = posTransform.fun, track.height = 0.07, col = "orange")

#textxx=head(text)
#textxx$V2=textxx$V2-5000
#textxx$V3=textxx$V3-5000

#text_list = list(text,textxx)

circos.genomicLabels(text, labels.column = 4, cex=.9,side = "outside",padding = 0.1, connection_height = convert_height(4, "mm"),col=c(rep("red",nrow(text)-10),rep("red",10)))
#circos.genomicLabels(text, labels.column = 4, cex=.6,side = "outside",padding = 0.001, connection_height = convert_height(7, "mm"), col="red")
#circos.genomicLabels(xx, labels.column = 1, cex=.6,side = "outside",padding = 0.02, connection_height = convert_height(7, "mm"))

cytoband = read.cytoband()$df
circos.genomicTrackPlotRegion(cytoband, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = cytoband.col(value[, 2]), border = NA, ...)
  cell.xlim = get.cell.meta.data("cell.xlim")
  cell.ylim = get.cell.meta.data("cell.ylim")
  circos.rect(cell.xlim[1], cell.ylim[1], cell.xlim[2], cell.ylim[2], border = "black")
  major.at = seq(0, cell.xlim[2], by = 50000000)
  major.labels = major.at/1000000
  l = major.at %% 50000000 == 0
  major.labels[l] = ""
  
  major.at.t = seq(0, cell.xlim[2], by = 25000000)
  major.labels.t = major.at/500000
  l.t = major.at.t %% 25000000 == 0
  major.labels.t[l.t] = ""
  circos.axis("top", major.at = major.at.t, labels = major.labels.t, labels.facing = "clockwise", labels.cex = 0.4, major.tick.percentage = 0.9)
  circos.text(major.at[l], rep(1.7, sum(l)), paste0(major.at[l]/1000000, ""), cex = 0.3, facing = "clockwise", adj = c(0, 0.5), niceFacing = TRUE)
}, bg.border = NA, track.height = 0.05)


bed_list = list(b1, b2,b3, b4)
col=c("#OOOOOO","#377EB8","#4DAF4A","#E41A1C")

circos.genomicRainfall(bed_list,pch = 16, cex = 0.4, ylim=c(0,7.5),col =c("#E41A1C","#377EB8","#4DAF4A","#000000"),track.height = 0.12)

bed_list = list(a1, a2,a3, a4)
circos.genomicRainfall(bed_list,pch = 16, cex = 0.4, ylim=c(0,7.5),col =c("#add8e6","#984EA3","#FF7F00","#A65628"),track.height = 0.12)

#circos.genomicDensity(b1, col = c("#E41A1C","#377EB8","#4DAF4A","#000000"),track.height=0.1)
circos.genomicDensity(bz, col = c("#4DAF4A"),track.height=0.1)
#circos.genomicTrack(h, 
#                    panel.fun = function(region, value, ...) {
#                      circos.genomicLines(region, value, type="h") #area=TRUE, col = "green", border="black"
#                    },track.height = 0.1)


circos.clear()

for(i in unique(cytoband$V1)) {minbp=min(cytoband[cytoband$V1==i,]$V2); maxbp=max(cytoband[cytoband$V1==i,]$V3);new=round((minbp+maxbp)/2); x=rbind(x,c(i,new,new+10,i))}
xx=data.frame(x)
xx=xx[grep("chr",xx$X1),]
xx$X4=gsub("chr","",xx$X4)
xx=head(xx,n=24)
write.table(xx,"~/Downloads/chr.pos.txt",sep="\t", row.names = F,col.names = F,quote=F)
xxx=read.table("~/Downloads/chr.pos.txt",sep="\t",header = F)

xxx = xxx %>% mutate(V2 = ifelse(V1=="chr11", V2-39000000,
                                 ifelse(V1=="chr15", V2-12000000,
                                        ifelse(V1=="chr4", V2+19000000,
                                               ifelse(V1=="chr16", V2-9000000,
                                                      ifelse(V1=="chr7", V2+9000000,
                                                             ifelse(V1=="chr1", V2+15000000,
                                                                    ifelse(V1=="chr13", V2+9000000,
                                                                           ifelse(V1=="chr20", V2-11000000,
                                                                                  ifelse(V1=="chr17", V2-9000000,V2))))))))))
xxx$V3=xxx$V2+10


############################################################################################




################################################################################################################
library(tidyr)
library(dplyr)



xxx=read.table("~/Downloads/chr.pos.txt",sep="\t",header = F)

xxx = xxx %>% mutate(V2 = ifelse(V1=="chr11", V2-55000000,
                                 ifelse(V1=="chr15", V2-12000000,
                                        ifelse(V1=="chr4", V2+25000000,
                                               ifelse(V1=="chr16", V2-900000,
                                                      ifelse(V1=="chr7", V2+9000000,
                                                             ifelse(V1=="chr1", V2+29000000,
                                                                    ifelse(V1=="chr13", V2+9000000,
                                                                           ifelse(V1=="chr20", V2-11000000,
                                                                                  ifelse(V1=="chr17", V2-9000000,V2))))))))))
xxx$V3=xxx$V2+10



h=read.table("~/Desktop/try_circlise.tsv",sep="\t",header = T)
colnames(h)=c("chr","start","end","value1")

#gnames=read.table("Downloads/sv_genes.txt",sep="\t",header=F)



text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)
ig.tra=read.table("/Users/yellapav/Downloads/IG_TRA_noncanonical_041219.txt",sep="\t",header = T)
#ig.tra=ig.tra[ig.tra$PCAWG_class=="Templated insertion",]
ig.tra=ig.tra[ig.tra$PCAWG_class!="Templated insertion",]

ig.tra = ig.tra %>% mutate(V3 = ifelse(chrom1=="22", pos1+5000000,
                                       ifelse(chrom2=="22", pos1+5000000,
                                              ifelse(chrom1=="14", pos1+5000000,
                                                     ifelse(chrom2=="14", pos1+5000000, pos1+2000000)))))

ig.tra = ig.tra %>% mutate(V4 = ifelse(chrom2=="22", pos2+5000000,
                                       ifelse(chrom1=="22", pos2+5000000,
                                              ifelse(chrom1=="14", pos2+5000000,
                                                     ifelse(chrom2=="14", pos2+5000000, pos2+2000000)))))


ig.tra = ig.tra %>% mutate(color.ig = ifelse(PCAWG_class=="Reciprocal translocation", "#FB8072",
                                             ifelse(PCAWG_class=="Chromoplexy", "#add8e6",
                                                    ifelse(PCAWG_class=="Chromothripsis", "#984EA3",
                                                           ifelse(PCAWG_class=="Complex", "#FF7F00",
                                                                  ifelse(PCAWG_class=="Unclassified translocation", "black",
                                                                         ifelse(PCAWG_class=="Local n distant jumps", "#377EB8",
                                                                                ifelse(PCAWG_class=="Templated insertion", "#A65628",
                                                                                       ifelse(PCAWG_class=="Unbalanced translocation", "#0a5d00",
                                                                                              ifelse(PCAWG_class=="Single", "#377EB8","‘#A9A9A9"))))))))))

#ig.tra$V3=ig.tra$pos1+5000000
b1=ig.tra[,c("chrom1","pos1","V3","color.ig")]
#b1$V3=b1$pos1+200000
b1$chrom1=paste0("chr",b1$chrom1)



b2=ig.tra[,c("chrom2","pos2","V4","color.ig")]
b2$chrom2=paste0("chr",b2$chrom2)

colnames(b1)=c("chr","start","end","value1")
colnames(b2)=c("chr","start","end","value1")

#h = h %>% mutate(value1 = ifelse(value1 > 40, 40,value1))



gnames=data.frame(grep("[A-Z]",unique(unlist(strsplit(as.vector(gsub("[/, ]","!",((ig.tra$gene)))),"!"))),value=T))
colnames(gnames)=c("V1")


bed=read.table("~/local/resources/GRCh37.e75.gene_boundaries.bed",sep="\t",header=F)

bed=separate(bed,V4,sep=";",c("V4","V5","V6"))
head(bed)
bed=left_join(gnames,bed,c("V1"="V6"))
bed=bed[,c(2,3,4,1)]
colnames(bed)=colnames(text)
bed$V1=paste0("chr",bed$V1)
bed$color="black"

bed1=rbind(bed,c("chr14",106032614,106032624,"IGH","red"))
bed1=rbind(bed1,c("chr22",23265085,23265095,"IGL","red"))
bed1=rbind(bed1,c("chr2",89156874,89156884,"IGK","red"))
bed1$V2=as.numeric(bed1$V2)
bed1$V3=as.numeric(bed1$V3)

bed1$chrome=as.numeric(gsub("chr","",bed1$V1))
bed1=arrange(bed1, chrome,V2)




library(circlize)

cnv=read.table("~/Downloads/commpass_cnv_new_nov_2019.txt",sep="\t",header = T) %>%
  dplyr::filter(sample=="MMRF_1890_1_BM")

cnvz = cnv %>% dplyr::mutate(chrom=paste0("chr",chrom)) %>% 
  dplyr::select(chrom,start,end,total) %>% dplyr::rename(value=total)

cnvz1 = cnv %>% dplyr::mutate(chrom=paste0("chr",chrom)) %>% 
  mutate(color = ifelse(total>2, "#FB8072","#add8e6")) %>% dplyr::select(chrom,start,end,total,color) %>% 
  dplyr::rename(value=total)

cnvz2 = cnv %>% dplyr::mutate(chrom=paste0("chr",chrom)) %>% 
  mutate(color = ifelse(total>2, "#FB8072","#add8e6")) %>% dplyr::select(chrom,start,end,minor,color) %>% 
  dplyr::rename(value=minor)

cnvtry = cnv %>% dplyr::mutate(chrom=paste0("chr",chrom)) %>% dplyr::mutate(value3=2) %>%
   dplyr::select(chrom,start,end,total,minor,value3) %>% 
  dplyr::rename(value1=total,value2=minor)

sv=read.table("~/Downloads/191112_delly_mapq_60_filtered_baseline.txt",sep="\t",header = T) %>%
  dplyr::filter(sample=="MMRF_1890_1_BM")


sv = sv %>% mutate(color = ifelse(PCAWG_class=="reciprocal_translocation", "#FB8072",
                                    ifelse(PCAWG_class=="chromoplexy", "#add8e6",
                                           ifelse(PCAWG_class=="chromothripsis", "#984EA3",
                                                  ifelse(PCAWG_class=="complex", "#FF7F00",
                                                         ifelse(PCAWG_class=="Unclassified translocation", "black",
                                                                ifelse(PCAWG_class=="local_n_jumps", "#377EB8",
                                                                       ifelse(PCAWG_class=="templated_insertion", "#A65628",
                                                                              ifelse(PCAWG_class=="unbalanced_inversion", "#0a5d00",
                                                                                     ifelse(PCAWG_class=="single", "#377EB8","#A9A9A9"))))))))))

zsv1 = sv %>% dplyr::mutate(chrom1=paste0("chr",chrom1), end=pos1+200000) %>% dplyr::select(chrom1,pos1,end,color)
zsv2 = sv %>% dplyr::mutate(chrom2=paste0("chr",chrom2), end=pos2+200000) %>% dplyr::select(chrom2,pos2,end,color)


  colnames(zsv1)=c("chr","start","end","value1")
  colnames(zsv2)=colnames(zsv1)
  
  
  

  
  head(zsv1)
  head(zsv2)
#par(lwd = 0.5)
#par(mar = c(-1, -1, -1, -1))

circos.par("cell.padding" = c(0, 0, 0, 0),canvas.xlim=c(-1,1),canvas.ylim=c(-1,1),"start.degree" = 90,"track.height" = 0.1)

circos.initializeWithIdeogram(plotType = NULL)

posTransform.fun = function(region) {
  return(region)
}
##Chrom numbers
circos.genomicTrackPlotRegion(xxx, ylim = c(2.1, 3.1), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 0.7,font=2, posTransform = posTransform.fun,niceFacing = TRUE, track.margin=c(0.01,0.01))
}, track.height = 0.05, bg.border = NA)


#gene labels
circos.genomicLabels(bed1, labels.column = 4, cex=0.7,side = "outside",padding = 0.2, connection_height = convert_height(2, "mm"),col=bed1$color)
#circos.genomicLabels(text, labels.column = 4, cex=.6,side = "outside",padding = 0.001, connection_height = convert_height(7, "mm"), col="red")
#circos.genomicLabels(xx, labels.column = 1, cex=.6,side = "outside",padding = 0.02, connection_height = convert_height(7, "mm"))

cytoband = read.cytoband()$df
circos.genomicTrackPlotRegion(cytoband, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = cytoband.col(value[, 2]), border = NA, ...)
  cell.xlim = get.cell.meta.data("cell.xlim")
  cell.ylim = get.cell.meta.data("cell.ylim")
  circos.rect(cell.xlim[1], cell.ylim[1], cell.xlim[2], cell.ylim[2], border = "black")
  major.at = seq(0, cell.xlim[2], by = 50000000)
  major.labels = major.at/1000000
  l = major.at %% 50000000 == 0
  major.labels[l] = ""
  
  major.at.t = seq(0, cell.xlim[2], by = 25000000)
  major.labels.t = major.at/500000
  l.t = major.at.t %% 25000000 == 0
  major.labels.t[l.t] = ""
  circos.axis("top", major.at = major.at.t, labels = major.labels.t, labels.facing = "clockwise", labels.cex = 0.4, major.tick.percentage = 0.9)
  circos.text(major.at[l], rep(1.7, sum(l)), paste0(major.at[l]/1000000, ""), cex = 0.3, facing = "clockwise", adj = c(0, 0.5), niceFacing = TRUE)
}, bg.border = NA, track.height = 0.05)


bed_list=c(cnvz1,cnvz2)
bed = generateRandomBed(nr = 100)
#circos.genomicTrackPlotRegion(bed_list, stack=TRUE,
#                              panel.fun = function(region, value, ...) {
#                                i=getI(...)
#  circos.genomicLines(region, value, col = i, ...)
#})

#circos.genomicTrackPlotRegion(cnvtry,ylim=c(0.25,2),panel.fun = function(region, value, ...) {
#  circos.genomicLines(region, value, col = 1:3,type="segment", ...)
#})

bed1= cnvz1 %>% dplyr::select(-color)
bed_list = list(bed1)
f = colorRamp2(breaks = c(0, 2, 3), colors = c("green", "black", "red"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                i = getI(...)
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              })


cnvtrr = cnvtry %>% dplyr::select(chrom,start,end,value1,value3)
circos.genomicTrackPlotRegion(cnvtrr, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type="segment",col = 1:2, ...)
})

circos.genomicLink(zsv1, zsv2,col = as.vector(zsv1$value1))

circos.clear()


circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 100, nc = 4)
cnvtrr = cnvtry %>% dplyr::select(chrom,start,end,value1,value3)
circos.genomicTrackPlotRegion(cnvtrr, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type="segment",col = 1:2, ...)
})
circos.clear()

circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram()

bed1 = generateRandomBed(nr = 100)
bed2 = generateRandomBed(nr = 100)
bed1= cnvz1 %>% dplyr::select(-color)
bed_list = list(bed1)
f = colorRamp2(breaks = c(0, 2, 3), colors = c("green", "black", "red"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                i = getI(...)
                                cell.xlim = get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              })

circos.clear()

circos.par("track.height" = 0.1)
circos.initializeWithIdeogram(plotType = NULL)

bed = generateRandomBed(nr = 100)
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", ...)
})

bed1 = generateRandomBed(nr = 100)
bed2 = generateRandomBed(nr = 100)
bed_list = list(bed1, bed2)

circos.genomicTrackPlotRegion(bed_list, panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicLines(region, value, col = i, ...)
})

circos.genomicTrackPlotRegion(bed_list, stack = TRUE, 
                              panel.fun = function(region, value, ...) {
                                i = getI(...)
                                circos.genomicLines(region, value, col = i, ...)
                              })

bed = generateRandomBed(nr = 100, nc = 4)
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, col = 1:4, ...)
})

circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicLines(region, value, col = i, ...)
})

bed = generateRandomBed(nr = 100)
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "segment", lwd = 2, ...)
})

circos.clear()



