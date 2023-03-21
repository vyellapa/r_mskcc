library(RCircos)
library(circlize)
library(dplyr)
library(tidyr)

text=read.table("~/Downloads/circos_sv_text.tsv",sep="\t",header = F)
xxx=read.table("~/Downloads/chr.pos.txt",sep="\t",header = F)



cnv=read.table("~/Downloads/commpass_cnv_new_nov_2019.txt",sep="\t",header = T) %>%
  dplyr::filter(sample=="MMRF_1890_1_BM")

cnv.tot = cnv %>% dplyr::mutate(chrom=paste0("chr",chrom)) %>% 
  mutate(color = ifelse(total>2, "#FB8072","#add8e6")) %>% 
  dplyr::select(chrom,start,end,total,color) %>% 
  dplyr::rename(value=total)

cnv.min = cnv %>% dplyr::mutate(chrom=paste0("chr",chrom)) %>% 
  mutate(color = ifelse(total>2, "#FB8072","#add8e6")) %>% 
  dplyr::select(chrom,start,end,minor,color) %>% 
  dplyr::rename(value=minor)


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

sv1 = sv %>% dplyr::mutate(chrom1=paste0("chr",chrom1), end=pos1+200000) %>% dplyr::select(chrom1,pos1,end,color)
sv2 = sv %>% dplyr::mutate(chrom2=paste0("chr",chrom2), end=pos2+200000) %>% dplyr::select(chrom2,pos2,end,color)

colnames(sv1)=c("chr","start","end","value1")
colnames(sv2)=colnames(sv1)
  

circos.par("cell.padding" = c(0, 0, 0, 0),canvas.xlim=c(-1,1),canvas.ylim=c(-1,1),"start.degree" = 90,"track.height" = 0.08)

circos.initializeWithIdeogram(plotType = NULL)

posTransform.fun = function(region) {
  return(region)
}

##Chrom numbers
circos.genomicTrackPlotRegion(xxx, ylim = c(1.8, 2.8), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 0.8,font=2, posTransform = posTransform.fun,niceFacing = TRUE, track.margin=c(0.01,0.01))
}, track.height = 0.05, bg.border = NA)


#gene labels
circos.genomicLabels(bed1, labels.column = 4, cex=0.7,side = "outside",padding = 0.2, connection_height = convert_height(2, "mm"),col=bed1$color)

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


bed.amp= cnv.tot %>% dplyr::select(-color) %>% filter(value>=2 )
bed.del= cnv.tot %>% dplyr::select(-color) %>% filter(value<=2)
bed_list = list(bed.amp,bed.del)
f = colorRamp2(breaks = c(0,1, 2,3,3.5), colors = c("#B8DC3C","#B8DC3C", "white", "#CB2402", "#CB2402"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE,
                              panel.fun = function(region, value, ...) {
                                
                                circos.genomicRect(region, value, col = f(value[[1]]), 
                                                   border = NA, ...)
                                i = getI(...)
                                # For a dashed line
                               # cell.xlim = get.cell.meta.data("cell.xlim")
                                # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#000000")
                              })


circos.genomicLink(sv1, sv2,col = as.vector(zsv1$value1))

circos.clear()


