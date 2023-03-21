library(dplyr)
library(GenomicRanges)
library("biomaRt")
library("RColorBrewer")
library(reshape2)
library(data.table)

setwd("/work/isabl/home/yellapav/projects/polyA/bams/output/ttests")

bb = read.table("../../../data/ipa.non.exonUTR.125bp.bed",header=F,sep="\t")
gr.seg = GRanges(seqnames=Rle(bb$V1), IRanges(as.numeric(as.character(bb$V2)), as.numeric(as.character(bb$V3))))
all.ipa = bb %>% dplyr::select(V1,V2,V3) %>% distinct()

i = "/work/isabl/home/yellapav/projects/polyA/bams/output/SRR4189473_IPAdepth.txt"
  print(i)
  r = fread(i,header=F,sep="\t")
  #r = distinct(r)
  r$sp = r$V2+1
  r$name = gsub("_IPAdepth.txt","",basename(i))
  name = gsub("_IPAdepth.txt","",basename(i))
  colnames(r) = c("chrom","st","dp","sp","name")
  
  gr.cyto = GRanges(seqnames=Rle(r$chrom), IRanges(r$st, r$sp))
  overlapGenes <- findOverlaps(gr.cyto, gr.seg)
  df1 = data.frame(bb[subjectHits(overlapGenes),], r[queryHits(overlapGenes),]) %>% distinct() 
  df1$pre = df1$V2+125
  df1$post = df1$V3-125
  
  
  pvall = list()
  for(j in 1:dim(all.ipa)[1]) {
    cat(c(j," "))
    
    temp.df = df1 %>% dplyr::filter(V1==all.ipa[j,1] & V2==all.ipa[j,2])
    pre = temp.df %>% dplyr::filter(st<pre & st>V2+25) 
    post = temp.df %>% dplyr::filter(st>post & st<V3-25) 
    if(nrow(pre)>5 & nrow(post)>5) {
    pval.temp = cbind(all.ipa[j,],as.data.frame(matrix(data.frame(unlist(t.test(pre$dp,post$dp)))[c(3,6,7),],nrow = 1)), name) 
    pvall[[j]] = pval.temp
    }
  }
  pval.df = do.call(rbind, pvall)
  




write.table(pval.df,"./SRR4189473_out.txt",col.names = F, row.names = F, append = F, quote = F, sep = "\t")

