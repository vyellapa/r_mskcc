#### 2 cohorts:
# Sanger cohort (PD as ID)
# chapman et al nature 2011 


#### upload libraries
library(RColorBrewer)
library(stringr)
library("GenomicRanges")
library(IRanges)
library(mclust)
library(readr) 
library(tidyr)
library(splitstackshape)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)
library(magrittr)
library(lmerTest)
library(knitr)
library(reshape2)
library(MASS)
library(lme4)
library(survival)
library(deconstructSigs)
library(VariantAnnotation) 
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(sjmisc)
library(IRanges)
library("GenomicRanges")
require(MASS)
library("emmeans")
library("seqinr")                # comp()
library("MutationalPatterns")    # cos_sim_matrix()
library(pheatmap)


#### Upload battenberg output for all samples

setwd("/Users/yellapav/Desktop/record/clock_paper/data")
ascat_final<- read.delim(file.path("copy_number.txt"),sep="\t", header=T, stringsAsFactors = F) 
cna2<- ascat_final[ascat_final$class!="normal",]

#### Upload clinical information for all samples

clin<- read.delim("age_all_chap_fm6.txt", sep="\t", stringsAsFactors = F)
first<- clin[clin$code=="first",]


#### upload all project 1212 and Chapman mutations with VAF and select only the first sample

myelo <- read.delim("Project1212Caveman.txt", sep="\t", stringsAsFactors = F, skip=66)
myelo2<- myelo[myelo$ASMD >140 & myelo$CLPM==0,]
myelo3<- myelo2[,c("Sample","Chrom","Pos","Ref","Alt","PM.Tum")]
# myelo3$Sample = substr(myelo3[,1],1,nchar(myelo3[,1])-1)
# sanger<- unique(myelo3[,1:5])
sanger<- myelo3[myelo3$Sample %in% first$Sanger.sample.ID, ]

chap <- read.delim("cave_pass_chapman.txt", sep="\t", stringsAsFactors = F)
chap2<- chap[chap$ASMD >84 & chap$CLPM==0,]
chap3<- chap2[,c("Sample","Chrom","Pos","Ref","Alt","PM.Tum")]
colnames(chap3)<- colnames(sanger)
cave<- rbind.data.frame(chap3, sanger)

# setwd("~/Desktop/clock/chap_fm6/")
# cave_all<- read.delim(file.path("cave_fm6_chap.txt"),sep="\t", header=T, stringsAsFactors = F) 
# sanger<- cave_all[grep("PD", cave_all$sample),]
# sanger$paz = substr(sanger[,1],1,nchar(sanger[,1])-1)
# mat<- unique(sanger[,c(6,2:5)])
# colnames(mat)[1]<-"sample"
# dana<- cave_all[-grep("PD", cave_all$sample),]
# cave<- rbind.data.frame(dana, mat)

#### Upload dirichlet data containing the trunk assigned to the "trunk" for each Sanger patients 

shared_cave<- read.delim("shared_clonal_mutations.txt", sep="\t", stringsAsFactors = F)

#### txt file with the NUMBER id CONVERSION FOR dp Chapman files

code_file<- read.delim("code_files.txt", sep="\t", stringsAsFactors = F)

#### create vector with the same file order

#setwd("/Volumes/GoogleDrive/My Drive/clock/data/files/intermediate/")
f<-c("1_DP_and_cluster_info.txt","2_DP_and_cluster_info.txt","3_DP_and_cluster_info.txt","4_DP_and_cluster_info.txt","5_DP_and_cluster_info.txt","6_DP_and_cluster_info.txt",
     "7_DP_and_cluster_info.txt","8_DP_and_cluster_info.txt","9_DP_and_cluster_info.txt","10_DP_and_cluster_info.txt","11_DP_and_cluster_info.txt","12_DP_and_cluster_info.txt",
     "13_DP_and_cluster_info.txt","14_DP_and_cluster_info.txt","15_DP_and_cluster_info.txt","16_DP_and_cluster_info.txt","17_DP_and_cluster_info.txt","18_DP_and_cluster_info.txt",
     "19_DP_and_cluster_info.txt","20_DP_and_cluster_info.txt","21_DP_and_cluster_info.txt")


#setwd("/Volumes/GoogleDrive/My Drive/clock/data/RDS_files///")
mut_l<- c("1_DPoutput_1000iters_300burnindataset.RData","2_DPoutput_1000iters_300burnindataset.RData","3_DPoutput_1000iters_300burnindataset.RData","4_DPoutput_1000iters_300burnindataset.RData",
          "5_DPoutput_1000iters_300burnindataset.RData","6_DPoutput_1000iters_300burnindataset.RData","7_DPoutput_1000iters_300burnindataset.RData","8_DPoutput_1000iters_300burnindataset.RData",
          "9_DPoutput_1000iters_300burnindataset.RData","10_DPoutput_1000iters_300burnindataset.RData","11_DPoutput_1000iters_300burnindataset.RData","12_DPoutput_1000iters_300burnindataset.RData",
          "13_DPoutput_1000iters_300burnindataset.RData","14_DPoutput_1000iters_300burnindataset.RData","15_DPoutput_1000iters_300burnindataset.RData","16_DPoutput_1000iters_300burnindataset.RData",
          "17_DPoutput_1000iters_300burnindataset.RData","18_DPoutput_1000iters_300burnindataset.RData","19_DPoutput_1000iters_300burnindataset.RData","20_DPoutput_1000iters_300burnindataset.RData",
          "21_DPoutput_1000iters_300burnindataset.RData")

#### create from DP intermdiate file the final file for the mol_time 
# compared to Kevin output the dan output doesn't have the rownames so you need this 
#  intermediate step to find the correct corrispondence

all_chap<- list()
for(i in (1:length(f)))
{
  ### get the file with each mutation annotated for each cluster
  setwd("/Users/yellapav/Desktop/record/clock_paper/data/files/intermediate/")
  sam<- read.delim(f[i], sep="\t", stringsAsFactors = F)
  
  ### reference mutation lost in the previous file
  setwd("/Users/yellapav/Desktop/record/clock_paper/data/files/RDS_files/")
  load(mut_l[i])
  mut_diri<-  as.data.frame.matrix(cbind(dataset[[1]], dataset[[2]]))
  diri_sam<- cbind.data.frame(mut_diri, sam)
  colnames(diri_sam)[1:2]<-c("Chrom","Pos")
  diri_sam$Pos<- as.numeric(as.character(diri_sam$Pos))
  diri_sam$Pos<- diri_sam$Pos +1
  diri_sam<- diri_sam[,c(1:2,ncol(diri_sam))]
  
  ### reference mutation lost in the previous file
  code_file_sam<- code_file[code_file$code == i, ] 
  cave_sam<- cave[cave$Sample == code_file_sam$sample, ]
  def_sam<- merge(diri_sam, cave_sam, by=c("Chrom","Pos"))
  def_sam<- def_sam[,c(4,1,2,5,6,7,3)]
  all_chap[[i]]<- def_sam
  
}
all_chap2<- do.call("rbind", all_chap)

### V0D92Z not possibile to reconstruct the tree


setwd("/Users/yellapav/Desktop/record/clock_paper/data/")
cluster_clin<- read.delim("sample_names_DIRIC.txt", sep="\t", stringsAsFactors = F, header=F)
cluster2<-cluster_clin[cluster_clin$V2=="trunk",]
# out <- str_split_fixed((cluster2[,1]),'_',2) 
all_chap2$code<- paste(all_chap2$Sample, all_chap2$most.likely.cluster, sep="_")
shared_chap<- all_chap2[all_chap2$code %in% cluster2$V1,]
shared_cave2<- shared_cave[grep("PD",shared_cave$Sample),]
shared_all<- rbind.data.frame(shared_chap[,c(1:6)],shared_cave2[,c(1:5,7)]) #### final clonal shared file for mol_time

### remove IGH regions

shared_all_igh<- shared_all[shared_all$Chrom == 14 & shared_all$Pos >106000000 &  shared_all$Pos< 107400000 | 
                              shared_all$Chrom == 22 & shared_all$Pos >22300000 &  shared_all$Pos< 23270000 | 
                              shared_all$Chrom == 2 & shared_all$Pos >89880000 &  shared_all$Pos< 90300000 |
                              shared_all$Chrom == 2 & shared_all$Pos >89100000 &  shared_all$Pos< 89700000,]
code<- rownames(shared_all_igh)
shared<- shared_all[!rownames(shared_all) %in% code,]


### V0D92Z and PD26434 no dirichlet data for mol_time (tot=49 pts)


#### filter battenberg for: 2:0. 2:1, 3:1, 4:1 status
ascat_final_cn2<- cna2[cna2$code_batt == "clonal",]
ascat_final_cn<- ascat_final_cn2[ascat_final_cn2$class_sub !="exclude",]
batt2<- ascat_final_cn[ascat_final_cn$class != "deletion" & ascat_final_cn$class !="Whole_chromosome_duplication" & ascat_final_cn$class !="LOH_3" &
                         ascat_final_cn$class != "LOH_5" & ascat_final_cn$class != "LOH_2"& ascat_final_cn$class != "bi_del" &
                         ascat_final_cn$class != "gain_6" & ascat_final_cn$class != "gain_7"& ascat_final_cn$class != "bi_del" & 
                         ascat_final_cn$class != "LOH_4" & ascat_final_cn$class != "LOH_High" & ascat_final_cn$class != "LOH_8" 
                       & ascat_final_cn$class != "LOH_9" & ascat_final_cn$class != "gain_4" & ascat_final_cn$class != "LOH_6" 
                       & ascat_final_cn$class != "-",]
batt<-batt2[(batt2$end-batt2$start)>1000000,]


clin<- read.delim("age_all_chap_fm6.txt", sep="\t", stringsAsFactors = F)
first_samples<- clin[clin$code=="first",]

### remove pts without gains 

sample_list<- sort(unique(first_samples$Sanger.sample.ID)[!unique(first_samples$Sanger.sample.ID) %in% 
                                                            c("PD26434c","PD26402a","V0D92Z","PD26425e", "PD26425f",
                                                              "PD26427a","PD26428a")])


ccf<- read.delim("ccf_samples.txt", sep="\t", stringsAsFactors = F)


####
#### functions for mol_time
####
  
spit = function (dbg, mess, ...) {
  if (dbg) { print( sprintf(mess,...) ) }
}

spat = function (dbg, name, var) {
  if (dbg) {print(name); print(var)}
}

oppstrand = function(x) {   # comp() complements dna sequences
  trin1 = paste(comp(rev(strsplit(substr(x,1,3),"")[[1]]), forceToLower=F), collapse="")
  trin2 = paste(comp(rev(strsplit(substr(x,5,7),"")[[1]]), forceToLower=F), collapse="")
  return(sprintf("%s>%s",trin1,trin2))
}

em_signatures = function(sigts.defn,mut.freqs,max.iter,dbg) {
  nof.sigts.defn <- ncol(sigts.defn)
  alpha = runif(nof.sigts.defn); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (i in 1:max.iter) {
    contr = t(array(alpha, dim=c(nof.sigts.defn,96))) * sigts.defn
    probs = contr/array(rowSums(contr), dim=dim(contr))
    probs[is.na(probs)] = 0
    probs = probs * mut.freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  spit(dbg, "em: exit iteration: %d", i)
  
  return( alpha )
}



clusters_plot <- function(sample_code, myMclust, first_def, cn_code, grey.thresh, eps=FALSE) {
  code<- unique(first_def$code)
  nof.clusts <- length(code)     #: nicos,    print(first_def$code); print(code);
  # here 17.09.04: resuffle the colours according to the mean ccf value of each cluster
  centroids <- numeric()
  for (t in 1:nof.clusts) {
    first_def_code <- first_def[first_def$code==code[t],]
    centroids[t] <- mean( first_def_code$ccf )
  }
  col<-c("brown3","forestgreen","cornflowerblue","blueviolet")[1:nof.clusts]
  col<-col[order(centroids)]
  col<-c(col,"grey50")
  max.clust <- max(code)
  grey.clust <- max.clust + 1
  
  par(mar=c(9,6,5,5), xpd=F)
  kk<- paste(sample_code,cn_code,sep="_")
  clust.prbs <- myMclust$z[,code]   # nicos: Mclust can have hidden clusters (remove them as they often are sub-clusts, with high co-prb
  ext.code <- c(code,grey.clust)   # fixme:
  
  if (is.matrix(clust.prbs) ) {   # fixme:
    for ( i in (1:nrow(clust.prbs)) ) { 
      row.i <- clust.prbs[i,]
      if (length(row.i[row.i>grey.thresh])>1 ) first_def$code[i] <- grey.clust;
    } 
  }
  
  pdf(sprintf("%s_clusters.pdf",kk), width=10, height=5) # was: "%s_evolution_plot_.pdf"
  clusters_plot_disp( nof.clusts, first_def, ext.code, col, cn_code )
  dev.off()
  
  if (eps) {
    setEPS()
    postscript(sprintf("%s_clusters.eps",kk),paper="special",family="NimbusSan") # also can try: ComputerModern
    clusters_plot_disp( nof.clusts, first_def, ext.code, col, cn_code )
    dev.off()
  }
  return(grey.clust)
}

clusters_plot_disp <- function( nof.clusts, first_def, ext.code, col, cn_code) {
  for ( z in (1:(nof.clusts+1)))  {
    first_def_col<- first_def[first_def$code == ext.code[z],]
    if (length(first_def_col$Pos) > 0) {
      plot(first_def_col$Pos, first_def_col$ccf, col=col[z],yaxt="n",xaxt="n", main=cn_code,ylab="", xlab="", ylim=c(0,150), pch=20,cex=1.5)
      par(new=TRUE)
    }
  }
  axis(2, at = seq(0,150, by=10),  col.axis = "black", labels=seq(0,150, by=10), lwd.ticks =1, cex.axis = 1.5, las=1)
}

# LOH
tetraploid.2.2.or.UPD.2.0.est <- function(ploidy.1, ploidy.2, bootstrap = TRUE, iter = 1000) {
  retv <- list()
  point.est <- ploidy.2 / (ploidy.2 + ploidy.1 / 2)
  if (bootstrap) {
    total.muts <- ploidy.1 + ploidy.2
    prob.2 <- ploidy.2 / total.muts
    resamps <- rbinom(n=iter, size=total.muts, p=prob.2)
    ests.bs <- resamps / (resamps + (total.muts - resamps)/2)
    ests.95 <- quantile(ests.bs, c(0.025,0.975))
    # return(c(point.est, ests.95))
    retv$pests <- c(point.est,ests.95)
    retv$noe   <- 1
    retv$boot1 <- ests.bs
  }
  else {retv$pests <- point.est}
  
  return(retv)
}

# 1 extra gain
triploid.2.1.est <- function(ploidy.1, ploidy.2, bootstrap = TRUE, iter=1000) {
  retv <- list()
  point.est <- ploidy.2 / (ploidy.2 + (ploidy.1 - ploidy.2)/3)
  if (bootstrap) {
    total.muts <- ploidy.1 + ploidy.2
    prob.2 <- ploidy.2 / total.muts
    resamps <- rbinom(n=iter, size=total.muts, p=prob.2)
    ests.bs <- resamps / (resamps + (total.muts - 2*resamps)/3)
    ests.95 <- quantile(ests.bs, c(0.025,0.975))
    # return(c(point.est, ests.95))
    retv$pests <- c(point.est,ests.95)
    retv$noe   <- 1
    retv$boot1 <- ests.bs
  }
  else {retv$pests <- point.est} # return(point.est)}
  
  return(retv)
}

# 2 extra gain
tetraploid.3.1.est <- function(ploidy.1, ploidy.2, ploidy.3, bootstrap=TRUE, iter=1000) {
  retv <- list()
  point.est.t1 <- ploidy.3 / (ploidy.3 + ploidy.2 + (ploidy.1 - 2*ploidy.2 - ploidy.3) / 4)
  point.est.t2 <- (ploidy.3 + ploidy.2) / (ploidy.3 + ploidy.2 + (ploidy.1 - 2*ploidy.2 - ploidy.3) / 4)
  if (bootstrap) {
    total.muts <- ploidy.1 + ploidy.2 + ploidy.3
    prob.2 <- ploidy.2 / total.muts
    prob.3 <- ploidy.3 / total.muts
    resamps <- rmultinom(n=iter, size=total.muts, p=c(1-prob.2-prob.3, prob.2, prob.3))
    ests.bs.t1 <- resamps[3,] / (resamps[3,] + resamps[2,] + (resamps[1,] - 2*resamps[2,] - resamps[3,])/4)
    ests.bs.t2 <- (resamps[3,] + resamps[2,]) / (resamps[3,] + resamps[2,] + (resamps[1,] - 2*resamps[2,] - resamps[3,])/4)
    ests.95.t1 <- quantile(ests.bs.t1, c(0.025,0.975))
    ests.95.t2 <- quantile(ests.bs.t2, c(0.025,0.975))
    retv$pests <- c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2)
    retv$noe   <- 2
    retv$boot1 <- ests.bs.t1
    retv$boot2 <- ests.bs.t2
    # return(c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2))
  }
  # else {return(c(point.est.t1, point.est.t2))}
  else {retv$pests <- c(point.est.t1, point.est.t2)}
  
  return(retv)
}

# 3 extra gain
tetraploid.4.1.est <- function(ploidy.1, ploidy.2, ploidy.3, ploidy.4,bootstrap=TRUE, iter=1000) {
  retv <- list()
  point.est.t1 <- (ploidy.4 * 5) / (ploidy.1 + ploidy.2 *2+ ploidy.3*3 +ploidy.4*4)
  point.est.t2 <- ((ploidy.4 +  ploidy.3)* 5) / (ploidy.1 + ploidy.2 *2+ ploidy.3*3 +ploidy.4*4)
  point.est.t3 <- ((ploidy.4 +  ploidy.3 + ploidy.2)* 5) / (ploidy.1 + ploidy.2 *2+ ploidy.3*3 +ploidy.4*4)
  if (bootstrap) {
    total.muts <- ploidy.1 + ploidy.2 + ploidy.3 + ploidy.4
    prob.2 <- ploidy.2 / total.muts
    prob.3 <- ploidy.3 / total.muts
    prob.4 <- ploidy.4 / total.muts
    resamps <- rmultinom(n=iter, size=total.muts, p=c(1-prob.2-prob.3- prob.4, prob.2, prob.3, prob.4))
    
    ests.bs.t1 <- ((resamps[4,]) * 5) / (resamps[1,] + resamps[2,] *2+ resamps[3,]*3 +resamps[4,]*4)
    ests.bs.t2 <- ((resamps[4,] +  resamps[3,])* 5) / (resamps[1,] + resamps[2,] *2+ resamps[3,]*3 +resamps[4,]*4)
    ests.bs.t3 <- ((resamps[4,] +  resamps[3,]+ resamps[2,])* 5) / (resamps[1,] + resamps[2,] *2+ resamps[3,]*3 +resamps[4,]*4)
    
    ests.95.t1 <- quantile(ests.bs.t1, c(0.025,0.975))
    ests.95.t2 <- quantile(ests.bs.t2, c(0.025,0.975))
    ests.95.t3 <- quantile(ests.bs.t3, c(0.025,0.975))
    
    retv$pests <- c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2,point.est.t3, ests.95.t3)
    retv$noe   <- 3
    retv$boot1 <- ests.bs.t1
    retv$boot2 <- ests.bs.t2
    retv$boot3 <- ests.bs.t3
    
    # return(c(point.est.t1, ests.95.t1, point.est.t2, ests.95.t2,point.est.t3, ests.95.t3))
  }
  # else {return(c(point.est.t1, point.est.t2,point.est.t3))}
  else {retv$pests <- c(point.est.t1,point.est.t2,point.est.t3)}
  
  return(retv)
}


mol_time_ic_gen <- function(all_code, iter=1000) {
  ###### generate molecular time and IC for gain
  pc<-list()
  bs<-list()  # nicos: bootstraps
  bs.code <- character()
  bs.segm <- character()
  bs.cnvn <- integer()
  bs.chrm <- character()
  bs.chvn <- character()
  j <- 1      # nicos: bootstraps index (accelerates faster than i)
  all_code_sgms <- all_code$segment
  ac_sgms_map <- list()
  for (i in 1:length(all_code_sgms)) {
    j = 1
    tchr = str_split_fixed(all_code_sgms[i],'-',4)[1]
    pt = paste0(tchr,letters[j])
    while (! is.null(ac_sgms_map[[pt]])) {
      j = j + 1  # let 's hope we don't run out of letters
      pt = paste0(tchr,letters[j])
      if (j>10) print(err2)
    }
    ac_sgms_map[[pt]] <- all_code_sgms[i]
  }
  all_code_gain<- all_code[all_code$cyto == "gain",]
  if(nrow(all_code_gain)>0)
  {
    for(i in (1:nrow(all_code_gain)))
    {
      CODE<-c(all_code_gain$cyto[i],all_code_gain$segment[i])
      est.res <- triploid.2.1.est(all_code_gain$cn1_clon[i],all_code_gain$cn2[i],iter=iter)
      pc_ic_row<- est.res$pests
      pc[[i]]<-c(CODE,pc_ic_row)
      
      for(k in(1:est.res$noe))
      { 
        bootN <- paste0("boot",k)
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,sep="_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    pc_plot_parameter_gain <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    gain_1<- as.data.frame.matrix(cbind(pc_plot_parameter_gain, matrix(, nrow = nrow(pc_plot_parameter_gain), ncol = 6)))
    colnames(gain_1)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {gain_1<-NULL}
  
  ###### generate molecular time and IC for gain_2
  
  pc<-list()
  all_code_gain_2<- all_code[all_code$cyto == "gain_2",]
  if(nrow(all_code_gain_2)>0)
  {
    for(i in (1:nrow(all_code_gain_2)))
    {
      CODE<-c(all_code_gain_2$cyto[i],all_code_gain_2$segment[i])
      est.res <- tetraploid.3.1.est(all_code_gain_2$cn1_clon[i],all_code_gain_2$cn2[i],all_code_gain_2$cn3[i],iter=iter)
      pc_ic_row<- est.res$pests
      pc[[i]]<-c(CODE,pc_ic_row)
      
      for(k in(1:est.res$noe))
      { 
        bootN <- paste("boot",k,sep="")
        # bs[[j]] <- c(CODE,est.res[bootN])
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,sep="_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    pc_plot_parameter_gain_2 <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    gain_2<- as.data.frame.matrix(cbind(pc_plot_parameter_gain_2, matrix(, nrow = nrow(pc_plot_parameter_gain_2), ncol = 3)))
    colnames(gain_2)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {gain_2<-NULL}
  ###### generate molecular time and IC for gain_3
  
  pc<-list()
  all_code_gain_3<- all_code[all_code$cyto == "gain_3",]
  if(nrow(all_code_gain_3)>0)
  {
    for(i in (1:nrow(all_code_gain_3)))
    {
      CODE<-c(all_code_gain_3$cyto[i],all_code_gain_3$segment[i])
      est.res   <- tetraploid.4.1.est(all_code_gain_3$cn1_clon[i],all_code_gain_3$cn2[i],all_code_gain_3$cn3[i],all_code_gain_3$cn4[i],iter=iter)
      pc_ic_row <- est.res$pests
      pc[[i]]   <- c(CODE,pc_ic_row)
      
      for(k in(1:est.res$noe))
      { 
        bootN <- paste("boot",k,sep="")
        # bs[[j]] <- c(CODE,est.res[bootN])
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,"_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    gain_3 <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    colnames(gain_3)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {gain_3<-NULL}
  
  ###### generate molecular time and IC for LOH
  pc<-list()
  all_code_gain_LOH<- all_code[all_code$cyto == "LOH",]
  if(nrow(all_code_gain_LOH)>0)
  {
    for(i in (1:nrow(all_code_gain_LOH)))
    {
      CODE<-c(all_code_gain_LOH$cyto[i],all_code_gain_LOH$segment[i])
      est.res   <- tetraploid.2.2.or.UPD.2.0.est(all_code_gain_LOH$cn1_clon[i],all_code_gain_LOH$cn2[i],iter=iter)
      pc_ic_row <- est.res$pests
      pc[[i]]<-c(CODE,pc_ic_row)
      
      for(k in(1:est.res$noe))
      { 
        bootN <- paste("boot",k,sep="")
        # bs[[j]] <- c(CODE,est.res[bootN])
        bs.code[j] <- CODE[1]
        bs.segm[j] <- CODE[2]
        bs.cnvn[j] <- k
        out3 <- str_split_fixed((CODE[2]),'-',4)
        bs.chrm[j] <- out3[,1]
        # bs.chvn[j] <- paste(out3[,1],k,sep="_")
        bs.chvn[j] <- paste(names(which(ac_sgms_map==CODE[2])),k,sep="_")
        bs[[j]] <- est.res[[bootN]]
        j <- j + 1
      }
    }
    pc_plot_parameter_loh <- as.data.frame.matrix(do.call("rbind", pc), stringsAsFactors = F) # merge results from all chromosomes  
    loh<- as.data.frame.matrix(cbind(pc_plot_parameter_loh))
    loh_1<- as.data.frame.matrix(cbind(pc_plot_parameter_loh, matrix(, nrow = nrow(pc_plot_parameter_loh), ncol = 6)))
    colnames(loh_1)<-c("code","segment","plot_dot1","IC_down_cn1","IC_up_cn1","plot_dot2","IC_down_cn2","IC_up_cn2","plot_dot3","IC_down_cn3","IC_up_cn3")
  }else
  {loh_1<-NULL}
  
  #### create the final file
  
  # all_code_final<- as.data.frame.matrix(rbind(gain_1,gain_2,gain_3,loh_1))
  retv <- list()
  # rboots <- as.data.frame.matrix(do.call("rbind", bs), stringsAsFactors = F) # cast the list to a matrix
  # retv$boots <- as.data.frame.matrix(cbind(rboots))
  boots <- data.frame( code=bs.code, chrom=bs.chrm, segment=bs.segm, cnvn=bs.cnvn, chrom_cnvn=bs.chvn )
  boots<- boots[complete.cases(boots),]
  for( k in (1:iter) ) {
    acc <- numeric()
    for(m in (1:(j-1))) {
      acc <- c(acc,bs[[m]][k])
    }
    boots <- cbind( boots, acc )
  }
  retv$boots <- boots
  retv$pests <- rbind(gain_1,gain_2,gain_3,loh_1)
  return(retv)
}

##### assign cluster to CN1, CN2, CN3, CN4 etc
clust_assign <- function(matrix, first_def, cn_code) {
  
  if(unique(first_def$class) == "gain") ###### single gain
  {
    if(nrow(matrix) == 2)
    {
      DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],0,0,as.numeric(sum(matrix$tot_mut)),cn_code)))
      colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      ## cyto[[i]] <- DF
    }else
    { if(nrow(matrix) == 3)
    {
      if(matrix$median[2]<((matrix$median[1] +matrix$median[3])/2)) ##### avarage between the two extreme cluster to assign the one in the middle
      {
        DF <- t(as.matrix(c(matrix$tot_mut[1] + matrix$tot_mut[2],matrix$tot_mut[3],0,0,as.numeric(sum(matrix$tot_mut)),cn_code)))
        colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      }
      else
      {
        DF <- t(as.matrix(c(matrix$tot_mut[1] , matrix$tot_mut[2] + matrix$tot_mut[3],0,0,as.numeric(sum(matrix$tot_mut)),cn_code)))
        colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      }
      ## cyto[[i]] <- DF
    }}}
  else
  {
    if(unique(first_def$class) == "gain_2") ###### 2 extra gains
    { if(nrow(matrix) == 2)
    {
      DF <- t(as.matrix(c(matrix$tot_mut[1],0,matrix$tot_mut[2],0,as.numeric(sum(matrix$tot_mut)),cn_code)))
      colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      ## cyto[[i]] <- DF
    }else
    {if(nrow(matrix) == 3)
    {
      DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],matrix$tot_mut[3],0,as.numeric(sum(matrix$tot_mut)),cn_code)))
      colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      ## cyto[[i]] <- DF
    }}}
    else
    { if(unique(first_def$class) == "gain_3") ###### 3 extra gains
    { if(nrow(matrix) == 3)
    {
      if((matrix$median[2]>(matrix$median[3]/4))&(matrix$median[2]<(matrix$median[3]/2)))
      {DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],0,matrix$tot_mut[3],as.numeric(sum(matrix$tot_mut)),cn_code)))
      colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      }else
      {if((matrix$median[2]>(matrix$median[3]/2))&(matrix$median[2]<(matrix$median[3]*3/4)))
      {
        DF <- t(as.matrix(c(matrix$tot_mut[1],0,matrix$tot_mut[2],matrix$tot_mut[3],as.numeric(sum(matrix$tot_mut)),cn_code)))
        colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      }
        else {stop("gain_3.point1") }
      }
      ## cyto[[i]] <- DF
    }else
    {if(nrow(matrix) == 4)
    {
      DF <- t(as.matrix(c(matrix$tot_mut[1],matrix$tot_mut[2],matrix$tot_mut[3],matrix$tot_mut[4],as.numeric(sum(matrix$tot_mut)),cn_code)))
      colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
      ## cyto[[i]] <- DF
    }}}else
    {
      if(unique(first_def$class) == "LOH")
      {if(nrow(matrix) == 2)
      {
        DF <- t(as.matrix(c(as.numeric(as.character(matrix$tot_mut[1])),as.numeric(as.character(matrix$tot_mut[2])),0,0,sum(as.numeric(as.character(matrix$tot_mut))),cn_code)))
        colnames(DF)<-c("cn1_clon","cn2","cn3", "cn4","tot","segment")
        ## cyto[[i]] <- DF
      }else
      {if(nrow(matrix) == 1) # fixme:
        stop("less than 2 clusters in LOH")
      }}
    }
    }}
  
  return(DF)
}

mol_time_warn <- function(warn,mess,arg) {
  if (warn) {
    print(sprintf(mess,arg))
  }
  return(TRUE)
}

mol_time_plot <- function (sample_code,all_code_final) {
  pdf(sprintf("%s_mol_time_estimate.pdf",sample_code), height=6, width=13)
  par(mar=c(10,7,3,3), xpd=F)
  par(mfrow=c(1,1))
  name_all2<-paste(all_code_final$chr, all_code_final$type, sep=" ")
  
  plot(c(1:nrow(all_code_final)),all_code_final$plot_dot1, ylim=c(0,2), xlim=c(0,(length(all_code_final$plot_dot1)+1)), xaxt="n", ylab="", 
       xlab="", pch=16, col=("dodgerblue"), 
       las=2, cex.axis=2, cex.lab = 1.5, cex=1.5, main=sample_code, cex.main=2.5)
  axis(1, at = c(1:nrow(all_code_final)),  col.axis = "black", labels=name_all2 , lwd.ticks =1, cex.axis = 2, las=2)
  mtext("Molecular Time", side = 2, line = 4, cex=3)
  
  par(new=TRUE)
  for(i in (1:nrow(all_code_final))) {
    segments(i,all_code_final$IC_down_cn1[i], i,all_code_final$IC_up_cn1[i],lwd=1, lty = 2, col="grey50")
    segments(i-0.2,all_code_final$IC_down_cn1[i],i+0.2,all_code_final$IC_down_cn1[i],lwd=1, lty = 1, col="grey50")
    segments(i-0.2,all_code_final$IC_up_cn1[i],i+0.2,all_code_final$IC_up_cn1[i],lwd=1, lty = 1, col="grey50")
  }
  
  second<-all_code_final
  second[is.na(second)]<-0
  par(new=TRUE)
  plot(c(1:nrow(all_code_final)) +0.1, all_code_final$plot_dot2, ylim=c(0,2), xlim=c(0,(length(all_code_final$plot_dot1)+1)), xaxt="n", ylab="", xlab="", pch=16, col=("red"), 
       las=2, cex.axis=2, cex.lab = 2, cex=1.5)
  
  num_gain_2<-which(all_code_final$plot_dot2!=0)
  
  for(w in (num_gain_2)) {
    segments(w+0.1,second$IC_down_cn2[w], w+0.1, second$IC_up_cn2[w], lwd=1, lty = 2, col="coral1")
    segments(w-0.1,second$IC_down_cn2[w],w+0.3,second$IC_down_cn2[w],lwd=1, lty = 1, col="coral1")
    segments(w-0.1,second$IC_up_cn2[w],w+0.3,second$IC_up_cn2[w],lwd=1, lty = 1, col="coral1")
    par(new=TRUE)
  }
  
  third<-all_code_final
  third[is.na(third)]<-0
  par(new=TRUE)
  plot(c(1:nrow(all_code_final)) +0.2, all_code_final$plot_dot3, ylim=c(0,2), xlim=c(0,(length(all_code_final$plot_dot3)+1)), xaxt="n", ylab="", xlab="", pch=16, col=("forestgreen"), 
       las=2, cex.axis=2, cex.lab = 2, cex=1.5)
  
  num_gain_2<-which(all_code_final$plot_dot3!=0)
  
  for(w in (num_gain_2)) {
    segments(w +0.2,third$IC_down_cn3[w], w+0.2, third$IC_up_cn3[w], lwd=1, lty = 2, col="forestgreen")
    segments(w-0.2+0.2,third$IC_down_cn3[w],w+0.2+0.2,third$IC_down_cn3[w],lwd=1, lty = 1, col="forestgreen")
    segments(w-0.2+0.2,third$IC_up_cn3[w],w+0.2+0.2,third$IC_up_cn3[w],lwd=1, lty = 1, col="forestgreen")
    par(new=TRUE)
  }
  
  dev.off()
  return(TRUE)
}



#######
####### clusterization part
#######

grey.thresh = 0.2
hclust.thresh=0.4
seed=3
iter=1000
warn=F
eps=F
annotated_mut_all2<- list()
shared$PM.Tum<- as.numeric(as.character(shared$PM.Tum))

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_18.01.2018/plot_final_19.03.11//")
for(j in (1:length(sample_list)))
{
  # j=11
  sample_code<- sample_list[j]
  print( sprintf("sample: %s", sample_code) )
  conc2<-ccf[ccf$Sample == sample_code,]
  aberrant<- as.numeric(conc2$CCF) ### extract the purity of the selected sample
  cave_sample<- shared[shared$Sample == sample_code,]
  ascat_sample_spec_type<- batt[batt$sample == sample_code,]
  cave_sample$chr<- paste0("chr",cave_sample$Chrom)
  ascat_sample_spec_type$chr<- paste0("chr",ascat_sample_spec_type$Chrom)
  gr0 = with(cave_sample, GRanges(chr, IRanges(start=(Pos), end=(Pos))))
  gr1 = with(ascat_sample_spec_type, GRanges(chr, IRanges(start=start, end=end)))
  ranges2 <- merge(as.data.frame(gr0),as.data.frame(gr1),by="seqnames",suffixes=c("A","B"))
  ranges2 <- ranges2[with(ranges2, startB <= startA & endB >= endA),]
  ascat_brass_second<- ranges2[,c("seqnames","startA","startB","endB")]
  colnames(ascat_brass_second)<- c("chr", "Pos","start","end")
  int2<- merge(ascat_sample_spec_type, ascat_brass_second, by=c("chr","start","end"))
  cave_ascat<- merge(int2, cave_sample, by=c("chr","Pos","Chrom"))
  cave_ascat<- unique(cave_ascat)
  if(nrow(cave_ascat)>0) {
    ##### create a code specific for each gain in order to filter out all gains with less than 50 SNV
    cave_ascat$ascat_code<-paste(cave_ascat$chr,cave_ascat$start, cave_ascat$end,cave_ascat$class,sep="-")
    alfa_sample<-as.data.frame(table(cave_ascat$ascat_code))
    alfa_sample2<-alfa_sample[alfa_sample$Freq>=50,]
    gain_mut<-as.character(alfa_sample2$Var1)
    cave_ascat_filt <- subset(cave_ascat, ascat_code %in% gain_mut)
  } else { mol_time_warn(warn,"no cave_ascat %s",sample_code)
    next  }     # goes to the next value of for loop
  if(nrow(cave_ascat_filt)>0) {
    all_code =NULL
  } else { mol_time_warn(warn,"no cave_ascat_filt %s",sample_code)
    next  }     # goes to the next value of for loop
  
  # setwd("/nfs/team78pc18/fm6/relapse_MM/molecular_time/plots")
  # setwd("/nfs/users/nfs_n/na11/ac/wtsi/proj/myeloma/mol_time")
  # setwd("/home/na11/ac/wtsi/proj/myeloma/mol_time/res.thres.0.1.wG2")
  DF<-NULL
  cyto<-list()
  mclust_mm<-list()
  chr_segment_list <-unique(cave_ascat_filt$ascat_code)
  
  # print( sprintf("for sample: %s", sample_code) )
  # print( sprintf("chr_segement_list: %i", length(chr_segment_list)) )
  
  
  annotated_mut_sample<- list()
  for(i in (1:length(chr_segment_list)))
  {
    cn_code<- chr_segment_list[i]
    chr1<- cave_ascat_filt[cave_ascat_filt$ascat_code ==chr_segment_list[i],]
    first<- chr1[,c("Sample", "Chrom","Pos","Ref","Alt", "PM.Tum","class")]
    first$Pos<- as.numeric(as.character(first$Pos))
    first$ccf<- (first$PM.Tum /aberrant)*100
    X<- first$ccf
    #mod5 = densityMclust(X)
    #plot(mod5, what = "density", type = "image", col = "dodgerblue3")  ##### add this if you want to see the density plot of the ccf
    myData <- data.frame(x=X)
    # get parameters for most optimal model. nicos: with at least 2 clusters
    
    if(unique(first$class) == "LOH"){ myMclust <- Mclust(myData,G=2,verbose=FALSE)}else
    { myMclust <- Mclust(myData,G=2:4,verbose=FALSE)}
    # add a column in myData CLUST with the cluster.
    myData$CLUST <- myMclust$classification
    colnames(myData)<-"ccf"
    first_def<- as.data.frame.matrix(cbind(first, myData[-1]))
    colnames(first_def)[ncol(first_def)]<- "code"
    grey.clust <- clusters_plot(sample_code, myMclust, first_def, cn_code, grey.thresh, eps=eps)
    first_def <- first_def[first_def$code != grey.clust,]
    
    
    ####### CN/allele assignment part
    
    number_clust<-sort(unique(first_def$code))
    median=NULL
    
    ##### estimate the median allelic fraction for each cluster
    
    # print( sprintf("for cluster sample: %s", sample_code) )
    # # print( sprintf("for number_clust: %i", length(number_clust)) )
    for(h in (1:length(number_clust))) {
      cluster<- first_def[first_def$code == number_clust[h],]
      alfa<- median(cluster$ccf)
      median<- c(median, alfa)
    } # for h
    code_up<-which.max(median)
    a<- table(first_def$code)
    
    matrix<- as.data.frame.matrix(cbind(median, rownames(a), a))
    colnames(matrix)<-c("median","cluster_number","tot_mut")
    matrix$tot_mut<-as.numeric(as.character(matrix$tot_mut))
    matrix$median<-as.numeric(as.character(matrix$median))
    matrix<-matrix[order(matrix$median),]
    
    cyto[[i]] <- clust_assign(matrix,first_def,cn_code)
    
    mclust_mm[[i]]<-first_def
    
    first_def_tab<- first_def
    first_def_tab$median<- NULL
    for(ww in (1:length(median)))
    {
      first_def_tab$median[first_def_tab$code ==ww]<- median[ww]
    }
    
    annotated_mut_sample[[i]]<- first_def_tab
    
    
  }   # for i
  
  annotated_mut_all<- do.call("rbind", annotated_mut_sample)
  
  # merge results from all chromosomes  and create a data frame with all SNVs and clusters
  
  mclust_mm_final<-     as.data.frame.matrix(do.call("rbind", mclust_mm), stringsAsFactors = F)
  #write.table(mclust_mm_final,sprintf("%s_with_cluster.txt",sample_code), quote=F, sep="\t", row.names=F, col.names = T)
  
  all_code <- as.data.frame.matrix(do.call("rbind", cyto), stringsAsFactors = F) # merge results from all chromosomes
  all_code[is.na(all_code)]<-0
  
  out2 <- str_split_fixed((all_code$segment),'-',4)
  all_code$cyto <- out2[,4]
  all_code$chr<-out2[,1]
  all_code$chrom <- gsub("chr","",all_code$chr)
  all_code<- all_code[order(all_code$chrom),]
  
  cols = c(1:5,ncol(all_code))
  all_code[,cols] = apply(all_code[,cols], 2, function(x) as.numeric(as.character(x)))
  
  file_code<- paste0(sample_code,"_cluster_absolute_number_summary.txt")
  #write.table(all_code,file_code,sep="\t", quote=F)
  
  
  ####### calculation part
  
  ####### Peter's script part bootstrap to generate IC and dots: one function for each possible gain combination
  # all_code_final<- as.data.frame.matrix(rbind(gain_1,gain_2,gain_3,loh_1))
  # all_code_final <- mol_time_ic_gen(all_code)$pests
  
  # nicos: introducing _comb as variables that contain the point estimates ($pests) and bootstrap runs ($bootN)
  
  all_code_final_comb <- mol_time_ic_gen(all_code,iter=iter)
  all_code_final <- all_code_final_comb$pests
  
  cols = c(3:ncol(all_code_final))
  all_code_final[,cols] = apply(all_code_final[,cols], 2, function(x) as.numeric(as.character(x)))
  
  out <- str_split_fixed((all_code_final$segment),'-',4)
  all_code_final$chr<- (out[,1])
  
  type <- str_split_fixed((out[,4]),'_',2)
  all_code_final$type<- type[,1]
  all_code_final$chrom<-gsub("chr","",all_code_final$chr)
  all_code_final$chrom[all_code_final$chrom=="X"]<-23
  all_code_final$chrom<-as.numeric(as.character(all_code_final$chrom))
  all_code_final<-all_code_final[order(all_code_final$chrom),]
  all_code_final$chrom[all_code_final$chrom==23]<-"X"
  
  ###### summary of pre and post gain mutation loads for each chromosome
  
  file_code<- paste0(sample_code,"_cluster_summary.txt")
  write.table(all_code_final,file_code,sep="\t", quote=F)
  
  mol_time_plot(sample_code,all_code_final);
  
  #all_code_final_hc<- all_code_final # nicos: not used anywhere ?!
  #all_code_final_hc[is.na(all_code_final_hc)]<-0
  
  if(nrow(all_code_final)>3){
    
    second<-all_code_final[ !is.na(all_code_final$plot_dot2),]
    third<-all_code_final[ !is.na(all_code_final$plot_dot3),]
    
    k<- as.matrix(c(all_code_final$plot_dot1, second$plot_dot2, third$plot_dot3))
    rownames(k)<- c(all_code_final$chr, second$chr, third$chr)
    write.table(k,paste0(sample_code,"_h_clust_tbl.txt"),sep="\t", quote=F)
    all_code2<-  (dist(k))
    pdf(sprintf("%s_h_clust.pdf",sample_code), height=6, width=8)
    plot(hclust(all_code2))
    dev.off()
    write.table(all_code_final_comb$boots,paste0(sample_code,"_boots.txt"),sep="\t", quote=F)
    hcls <- character()
    
    for( i in (1:iter)){
      
      hmx <- as.matrix( all_code_final_comb$boots[[5+i]] )
      rownames(hmx) <- all_code_final_comb$boots[[5]]
      hds <- dist(hmx)
      hcl <- hclust(hds)
      hct <- cutree(hcl, h=hclust.thresh)
      hct <- hct[order(names(hct))]   # align along lexico
      # we probably need to standarise cluster naming ?
      hcls[i] <- paste( hct, collapse="_" )
    }   # for i
    
    thcl <- table(hcls)
    
    max.len <- max( length(hct), length(thcl) )
    hnas <- rep(NA,max.len-length(thcl))
    
    edf <- data.frame( seg_cnvn=c(names(hct),rep(NA,max.len - length(hct))), hcls=c(names(thcl),hnas), hfreq=c(as.vector(thcl),hnas) )
    write.table(edf,paste0(sample_code,"_hfreqs.txt"),sep="\t", quote=F)
    
    nms_thcl <- names(thcl)
    cncl <- matrix(nrow=length(nms_thcl),ncol=length(as.integer(strsplit(nms_thcl[1],split="_")[[1]])) )
    # cncl <- matrix()
    
    for (i in 1:length(nms_thcl)) {
      cncl[i,] <- as.integer(strsplit(nms_thcl[i],split="_")[[1]])
    }
    
    colnames(cncl) <- names(hct)
    max.len.2 <- max( length(cncl), length(thcl) )
    hnas.2 <- rep(NA,max.len.2-length(thcl))
    hfdf <- data.frame(cncl, hfreq=c(as.vector(thcl),hnas.2))
    write.table(hfdf,paste0(sample_code,"_chrom_clust_freqs.txt"),sep="\t", quote=F)
    
  }  # IF THEN
  
  annotated_mut_all2[[j]]<- annotated_mut_all
}  # for j

duple<- do.call("rbind", annotated_mut_all2)
duple$sample_code<- paste(duple$Sample, duple$Chrom, sep="_")

### manual correction

duple<- duple[duple$Sample != "MMRC02G7BM" & duple$Chrom != 15,] ### too much uncertainty
duple<- duple[duple$Sample != "MMRC02G7BM" & duple$Chrom != 12,] ### too much uncertainty


# ### PD26400a chr 1has the first cluster with 0 mutations. This means that for some reason everything is moved up...
# duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 1"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 1"] == 2]<- 1
# duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 1"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 1"] == 3]<- 2
# duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 1"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 1"] == 4]<- 3


# ### PD26400a chr 18 non duplicated are splitted in 2. This means that for some reason everything is moved up...
duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 18"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 18"] == 1]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 18"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 18"] == 2]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 18"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26400a 18"] == 3]<- 2

### PD26410a chr 11 has the first cluster with 0 mutations. This means that for some reason everything is moved up...
duple$code[paste(duple$Sample, duple$Chrom) == "PD26410d 11"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26410d 11"] == 2]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26410d 11"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26410d 11"] == 3]<- 2
duple$code[paste(duple$Sample, duple$Chrom) == "PD26410d 11"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26410d 11"] == 4]<- 3

### PD26419a chr 3 has the first cluster with 0 mutations. This means that for some reason everything is moved up...
duple$code[paste(duple$Sample, duple$Chrom) == "PD26419a 3"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26419a 3"] == 1]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26419a 3"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26419a 3"] == 2]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26419a 3"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26419a 3"] == 3]<- 2

### PD26432c chr 3 has the first cluster with 0 mutations. This means that for some reason everything is moved up...
duple$code[paste(duple$Sample, duple$Chrom) == "PD26432c 7"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26432c 7"] == 1]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26432c 7"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26432c 7"] == 2]<- 1
duple$code[paste(duple$Sample, duple$Chrom) == "PD26432c 7"][duple$code[paste(duple$Sample, duple$Chrom) == "PD26432c 7"] == 3]<- 2



# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# 
# write.table(duple, "mol_time_output_march.txt", sep="\t", quote=F)
#setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
duple<- read.delim("mol_time_output_march.txt", sep="\t", stringsAsFactors = F, header=T)



###########################################################################################################################################################
### combine all files
############################################################################################################################################################


#### this is the manually curated output from mol_time function - it contains mol_time and IC fo each CN gains

#setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
summary_all2<- read.delim("order_matrix.txt", sep="\t", stringsAsFactors = F) 

gr0 = with(summary_all2, GRanges(chr, IRanges(start=(start), end=(end))))
values(gr0) <- DataFrame(code_ti = summary_all2$segment, 
                         cnv = summary_all2$cnv,
                         plot_dot1 = summary_all2$plot_dot1,
                         Sample = summary_all2$Sample,
                         type=summary_all2$code,
                         order = summary_all2$order,
                         double_chr<- summary_all2$code_extra_gains)

duple$chr<- paste0("chr",duple$Chrom)
duple$Pos<- as.numeric(as.character(duple$Pos))
gr1 = with(duple, GRanges(chr, IRanges(start=Pos, end=Pos)))
values(gr1) <- DataFrame(Sample = duple$Sample, 
                         Ref = duple$Ref,
                         Alt = duple$Alt,
                         PM.Tum = duple$PM.Tum,
                         ccf = duple$ccf,
                         class = duple$class,
                         sample_code = duple$sample_code,
                         median = duple$median,
                         code=duple$code)

ranges2 <- merge(as.data.frame(gr1),as.data.frame(gr0),by="seqnames",suffixes=c("A","B"))
ranges2 <- ranges2[with(ranges2, startB <= startA & endB >= endA),]
ranges3<- ranges2[which(ranges2$SampleA == ranges2$SampleB),]
ranges3$mol_order<- paste(ranges3$SampleA, gsub("chr","",ranges3$seqnames), sep="_") #### sample ID pasted to the the chromosome
colnames(ranges3)[6]<-"Sample"

fin2<- ranges3[,c("Sample","seqnames","startA","Ref","Alt","PM.Tum","code","ccf","class","median","order","double_chr....summary_all2.code_extra_gains")]
head(fin2)
colnames(fin2)[1:5]<-c("Sample","Chrom","Pos","Ref","Alt")

#########################################################################################################
### test signature for time window or for chromosome
#########################################################################################################
### upload trinucleotide context for each mutation for each dirichlet cluster in patients with more than one samples
#setwd("~/Desktop/clock/MM_tree_structure//")
mutations<-readRDS("all_dirichlet_MM_cave_trinucl_for_extraction.RDS")
head(mutations)
out <- str_split_fixed((mutations$sampleID),'_',2) 
mutations$Sample<- out[,1]

mutations2<- mutations[,c(1:5,10)]

all_chap2_5<- all_chap2[,c(8,2:5,1)]
colnames(all_chap2_5)<-colnames(mutations2)[1:6]

diri<- rbind(mutations2, all_chap2_5) #### combine dirichlet data from chapman and fm6 to get also the subclonal fraction
diri<- unique(diri[,1:6])
diri$chr<- paste0("chr",diri$chr)

colnames(fin2)[2:5]<- colnames(diri)[2:5]

# myelo3$Sample = substr(myelo3[,1],1,nchar(myelo3[,1])-1)

fin3<- fin2
colnames(fin3)[1:5]<-c("Sample","chr","pos","ref","mut")

fin3$Sample[(grep("PD", fin3$Sample))] = substr(fin3$Sample[(grep("PD", fin3$Sample))],1,nchar(fin3$Sample[(grep("PD", fin3$Sample))])-1)
require(plyr)
int_diri<- join(diri, fin3, by=c("Sample","chr","pos","ref","mut"))

##################################################################################################################################################################
### 96 classes for each patient DP considering each cluster as independent 
############################################################################################################################################################

int_diri_DP<- int_diri[,c(1:5)]
colnames(int_diri_DP)<-  c("sample","chr","pos","ref","alt")
DP<- mut.to.sigs.input(mut.ref = int_diri_DP,
                       sample.id = "sample",
                       chr = "chr",
                       pos = "pos",
                       ref = "ref",
                       alt = "alt",
                       bsg = BSgenome.Hsapiens.UCSC.hg19)

mm_col <- c(rep("lightskyblue",16), rep("black",16), rep("firebrick2",16),rep("gray88",16),rep("darkolivegreen3",16),
            rep("lightpink1",16))


# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(DP, "96_dirichlet_per_samples.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# DP<- read.delim("96_dirichlet_per_samples.txt", sep="\t", stringsAsFactors = F)
#  

##### only pre post gainand subclonal signature contribution where you have a gain (clonal and subclonal SNVs)

summary_all_u<- summary_all2
summary_all_u$Sample[(grep("PD", summary_all_u$Sample))] = substr(summary_all_u$Sample[(grep("PD", summary_all_u$Sample))],
                                                                  1,nchar(summary_all_u$Sample[(grep("PD", summary_all_u$Sample))])-1)

gr0 = with(summary_all_u, GRanges(chr, IRanges(start=(start), end=(end))))
values(gr0) <- DataFrame(code_ti = summary_all_u$segment,
                         cnv = summary_all_u$cnv,
                         plot_dot1 = summary_all_u$plot_dot1,
                         Sample = summary_all_u$Sample)

gr1_pre_gain = with(int_diri, GRanges(chr, IRanges(start=pos, end=pos)))
values(gr1_pre_gain) <- DataFrame(Sample = int_diri$Sample,
                                  duple_ID= int_diri$sampleID,
                                  ref = int_diri$ref,
                                  mut = int_diri$mut,
                                  PM.Tum = int_diri$PM.Tum,
                                  ccf = int_diri$ccf,
                                  class = int_diri$class,
                                  median = int_diri$median,
                                  code=int_diri$code)

ranges2_gain2 <- merge(as.data.frame(gr1_pre_gain),as.data.frame(gr0),by="seqnames",suffixes=c("A","B"))
ranges2_gain <- ranges2_gain2[with(ranges2_gain2, startB <= startA & endB >= endA),]
ranges3_clon_sub<- ranges2_gain[which(ranges2_gain$SampleA == ranges2_gain$SampleB),]

ranges3_clon_sub$code[is.na(ranges3_clon_sub$code)]<-"subclonal"
ranges3_clon_sub$code[ranges3_clon_sub$code == "3"]<-"2"
ranges3_clon_sub$gain_code<- paste(ranges3_clon_sub$SampleA, ranges3_clon_sub$code, sep="_")

mol_time_sig<- ranges3_clon_sub[,c("gain_code","seqnames","startA","ref","mut")]
mol_time_sig<-(unique(mol_time_sig))

colnames(mol_time_sig)<-  c("sample","chr","pos","ref","alt")
mol_time_sig_96<- mut.to.sigs.input(mut.ref = mol_time_sig,
                                    sample.id = "sample",
                                    chr = "chr",
                                    pos = "pos",
                                    ref = "ref",
                                    alt = "alt",
                                    bsg = BSgenome.Hsapiens.UCSC.hg19)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(mol_time_sig_96, "96_pre_post_gain_per_samples.txt", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(mol_time_sig_96, "96_pre_post_gain_per_samples.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# mol_time_sig_96<- read.delim("96_pre_post_gain_per_samples.txt", sep="\t", stringsAsFactors = F)


##### look for intergains
##### 

ranges3_clon_sub<- ranges2_gain[which(ranges2_gain$SampleA == ranges2_gain$SampleB),]
ranges3_clon_sub$code[is.na(ranges3_clon_sub$code)]<-"subclonal"
ranges3_clon_sub$gain_code_int<- paste(ranges3_clon_sub$SampleA, ranges3_clon_sub$code, sep="_")
mol_time_sig<- ranges3_clon_sub[,c("gain_code_int","seqnames","startA","ref","mut")]
colnames(mol_time_sig)<-  c("sample","chr","pos","ref","alt")
mol_time_sig_96_int<- mut.to.sigs.input(mut.ref = mol_time_sig,
                                        sample.id = "sample",
                                        chr = "chr",
                                        pos = "pos",
                                        ref = "ref",
                                        alt = "alt",
                                        bsg = BSgenome.Hsapiens.UCSC.hg19)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(mol_time_sig_96_int, "96_pre_post_gain_per_samples_int_samples.txt", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(mol_time_sig_96_int, "96_pre_post_gain_per_samples_int_samples.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# mol_time_sig_96_int<- read.delim("96_pre_post_gain_per_samples_int_samples.txt", sep="\t", stringsAsFactors = F)




##### 
##### look for LOH only
##### 

ranges_loh<- ranges3[ranges3$cnv == "LOH",]
ranges_loh$gain_code_int<- paste(ranges_loh$Sample, ranges_loh$seqnames, ranges_loh$code, sep="_")
mol_time_sig_LOH<- ranges_loh[,c("gain_code_int","seqnames","startA","Ref","Alt")]
colnames(mol_time_sig_LOH)<-  c("sample","chr","pos","ref","alt")
mol_time_sig_LOH<- mut.to.sigs.input(mut.ref = mol_time_sig_LOH,
                                     sample.id = "sample",
                                     chr = "chr",
                                     pos = "pos",
                                     ref = "ref",
                                     alt = "alt",
                                     bsg = BSgenome.Hsapiens.UCSC.hg19)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(mol_time_sig_LOH, "96_pre_post_gain_LOH_per_samples_int_samples.txt", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# # write.table(mol_time_sig_LOH, "96_pre_post_gain_LOH_per_samples_int_samples.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# mol_time_sig_LOH<- read.delim("96_pre_post_gain_LOH_per_samples_int_samples.txt", sep="\t", stringsAsFactors = F)
# 

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### 
#####  pre post gain with chromosome order per time window (only clonal SNVs)
##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
fin2<- ranges3[,c("Sample","seqnames","startA","Ref","Alt","PM.Tum","code","ccf","class","median","order","code_ti")]
head(fin2)
colnames(fin2)[1:5]<-c("Sample","Chrom","Pos","Ref","Alt")
single_chrom<- fin2
single_chrom$Sample[(grep("PD", single_chrom$Sample))] = substr(single_chrom$Sample[(grep("PD", single_chrom$Sample))],
                                                                1,nchar(single_chrom$Sample[(grep("PD", single_chrom$Sample))])-1)

single_chrom$sampleID<- paste(single_chrom$Sample, single_chrom$code , single_chrom$Chrom, single_chrom$code_ti, sep="_") ## code is the pre-post gain flag / order is the time window order per chr
single_chrom2<- single_chrom[,c("sampleID","Chrom","Pos","Ref","Alt")]
colnames(single_chrom2)<-  c("sample","chr","pos","ref","alt")
mol_time_single_chomr<- mut.to.sigs.input(mut.ref = single_chrom2,
                                          sample.id = "sample",
                                          chr = "chr",
                                          pos = "pos",
                                          ref = "ref",
                                          alt = "alt",
                                          bsg = BSgenome.Hsapiens.UCSC.hg19)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(mol_time_single_chomr, "96_pre_post_gain_per_chromosome_Sample.txt", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(mol_time_single_chomr, "96_pre_post_gain_per_chromosome_Sample.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# mol_time_single_chomr<- read.delim("96_pre_post_gain_per_chromosome_Sample.txt", sep="\t", stringsAsFactors = F)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### 
##### all mutations to define signature yes or no
##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

myelo2<- myelo[myelo$ASMD >140 & myelo$CLPM==0,]
myelo3<- myelo2[,c("Sample","Chrom","Pos","Ref","Alt")]
myelo3$Sample[myelo3$Sample == "WGA_PD26405c"]<-"PD26405c"
myelo4<- myelo3[myelo3$Sample %in% clin$Sanger.sample.ID,]
chap2<- chap[chap$ASMD >84 & chap$CLPM==0,]
chap3<- chap2[,c("Sample","Chrom","Pos","Ref","Alt")]
colnames(chap3)<- colnames(myelo3)
cave_all_mut<- rbind.data.frame(chap3, myelo4)

colnames(cave_all_mut)<-  c("sample","chr","pos","ref","alt")
cave_all_mut$chr<- paste0("chr", cave_all_mut$chr)
all_sig_96<- mut.to.sigs.input(mut.ref = cave_all_mut,
                               sample.id = "sample",
                               chr = "chr",
                               pos = "pos",
                               ref = "ref",
                               alt = "alt",
                               bsg = BSgenome.Hsapiens.UCSC.hg19)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(all_sig_96, "96_all_mutations_per_samples.txt", sep="\t", quote=F)
#
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(all_sig_96, "96_all_mutations_per_samples.txt", sep="\t", quote=F)
#
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# all_sig_96<- read.delim("96_all_mutations_per_samples.txt", sep="\t", stringsAsFactors = F)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### 
#####  signatures investigation using fitting approach 
##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# 
# 


### FIRST FITTING IS DONE ON ALL CASES WITH ALL SIGNATURES EXTRACTED BY NNMF AND HDP; 
### FROM THIS WE WILL CREATE THE SIGNATURE CATALOGUE FOR EACH SAMPLES THAT WILL BE USED FOR ALL NEXT ANALYSIS

dbg=FALSE
consigts.defn <- read.csv("mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(all_sig_96)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total

samples.muts<-samples.muts[,colSums(samples.muts)>100]
samples<- colnames(samples.muts)



consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)

sample.sigt.profs=NULL

if (is.null(sample.sigt.profs)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  sigt.profs <- list()
  for (i in 1:length(samples)) {
    sigt.profs[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}

sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000



for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  # if(length(grep("PD",samples[j]  ))==1){
  #   code_sample= substr(samples[j],1,nchar(samples[j])-1)
  # }else
  #   {code_sample<- samples[j] }
  code_sample<- samples[j]
  sample.sigts <- unique(sigt.profs[[ code_sample]])
  
  # we don't need the following anymore .... but it also sorts the vector
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  
  
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    # if(!is.na(rem.names))
    names(cosReduction) <- rem.names
    if (min(cosReduction) < 0.01) {
      spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
      rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
      rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      reducing = TRUE
    }else {spit(dbg, "exiting while...")
      reducing = FALSE
    }
  }
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)
sig<- sig[order(rownames(sig)),]
#   
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(sig, "all_patients_contribution.txt", sep="\t", quote=F)
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(sig, "all_patients_contribution.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# sig<- read.delim("all_patients_contribution.txt", sep="\t", stringsAsFactors = F)
head(sig)
sig_catalogue<- sig[,-ncol(sig)]  
sig_catalogue$sample<- rownames(sig_catalogue)
sig_catalogue$paz<- rownames(sig_catalogue)
sig_catalogue$paz[grep("PD",sig_catalogue$paz)] = substr(sig_catalogue$sample[grep("PD",sig_catalogue$paz)]
                                                         ,1,nchar(sig_catalogue$sample[grep("PD",sig_catalogue$paz)])-1)

#### create list with annotated which signatures are active for each sample - this will be used for the extraction on different contexts

int_data_framne<- NULL
sample_list<- unique(sig_catalogue$paz)
for(i in (1:length(sample_list)))
{
  sig_catalogue_sing<- sig_catalogue[sig_catalogue$paz == sample_list[i], ]
  if(nrow(sig_catalogue_sing)==1){
    int_data_framne<- rbind(int_data_framne, sig_catalogue_sing[,-9])
  }else{
    int_data_framne<- rbind(int_data_framne, c(colSums(sig_catalogue_sing[,1:8]), unique(sig_catalogue_sing[,10])))
  }
}

int_data_framne2<- as.data.frame.matrix(int_data_framne)
rownames(int_data_framne2)<- int_data_framne2$paz

sig_catalogue_fin<- int_data_framne2[,-ncol(int_data_framne2)]
for(i in (1:ncol(sig_catalogue_fin)))
{
  sig_catalogue_fin[,i][sig_catalogue_fin[,i] !=0]<- colnames(sig_catalogue_fin)[i]
  sig_catalogue_fin[,i][sig_catalogue_fin[,i] ==0]<- NA
}


sample.sigt.profs<- list()
for(i in (1:nrow(sig_catalogue_fin)))
{
  kk<- as.character(as.character((sig_catalogue_fin[i,])))
  sample.sigt.profs[[rownames(sig_catalogue_fin)[i]]]<- kk[!is.na(kk)]
}


##################################################################################################
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### 
#####  signatures investigation for unique case
##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

myelo2<- myelo[myelo$ASMD >140 & myelo$CLPM==0,]
myelo3<- myelo2[,c("Sample","Chrom","Pos","Ref","Alt")]
myelo3$Sample[myelo3$Sample == "WGA_PD26405c"]<-"PD26405c"
myelo4<- myelo3[myelo3$Sample %in% clin$Sanger.sample.ID,]
myelo4[,1] = substr(myelo4[,1],1,nchar(myelo4[,1])-1)

chap2<- chap[chap$ASMD >84 & chap$CLPM==0,]
chap3<- chap2[,c("Sample","Chrom","Pos","Ref","Alt")]
colnames(chap3)<- colnames(myelo3)
cave_all_mut<- rbind.data.frame(chap3, myelo4)
cave_all_mut<- unique(cave_all_mut[,1:5])
colnames(cave_all_mut)<-  c("sample","chr","pos","ref","alt")
cave_all_mut$chr<- paste0("chr", cave_all_mut$chr)
all_sig_96_unique<- mut.to.sigs.input(mut.ref = cave_all_mut,
                                      sample.id = "sample",
                                      chr = "chr",
                                      pos = "pos",
                                      ref = "ref",
                                      alt = "alt",
                                      bsg = BSgenome.Hsapiens.UCSC.hg19)


# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(all_sig_96_unique, "96_all_mutations_per_unique_paz.txt", sep="\t", quote=F)
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(all_sig_96_unique, "96_all_mutations_per_unique_paz.txt", sep="\t", quote=F)
# 

consigts.defn <- read.csv("mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(all_sig_96_unique)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total

samples.muts<-samples.muts[,colSums(samples.muts)>100]
samples<- colnames(samples.muts)

consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)
sample.sigt.profs2=NULL
if (is.null(sample.sigt.profs2)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  sigt.profs2 <- list()
  for (i in 1:length(samples)) {
    sigt.profs2[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}

sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000

for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  ### split smapleID pre-post gain data to identify the correct sample
  out <- str_split_fixed(( samples[j]),'_',2)
  samples_list<- out[,1]
  ###
  sample.sigts <- unique(sigt.profs2[[ samples_list ]])
  sample.sigts <- c("Signature.Subs.01",sample.sigts)
  
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    
    ### test the cos-similarities of each signature but SBS1 and SBS5
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    # if(!is.na(rem.names))
    names(cosReduction) <- rem.names
    if (min(cosReduction) < 0.01) {
      spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
      rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
      rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      reducing = TRUE
    }else {spit(dbg, "exiting while...")
      reducing = FALSE
    }
  }
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig_u <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(sig_u, "all_unique_patients", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(sig_u, "all_unique_patients", sep="\t", quote=F)
# 
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# sig_u<- read.delim( "all_unique_patients", sep="\t", stringsAsFactors = F)
# # 
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# sig<- read.delim("all_patients_contribution.txt", sep="\t", stringsAsFactors = F)
sig_all<- sig_u
sig_all$V2<- rownames(sig_all)

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/transplanat/") ### CLINICAL DATA FOR ALL CASES
trans<- read.delim("sct_data_2.txt", sep="\t", stringsAsFactors = F)

fin_all<- merge(sig_all,trans, by="V2")
fin_all$col_code<- paste(fin_all$stage, fin_all$transplant)
fin_all$col_code2<- fin_all$col_code

fin_all$col_code2[fin_all$col_code2 == "RR-MM NA"] <-"dodgerblue3"
fin_all$col_code2[fin_all$col_code2 == "ND MM "] <-"brown2"
fin_all$col_code2[fin_all$col_code2 =="SMM no_SCT"] <- "forestgreen"
fin_all$col_code2[fin_all$col_code2 =="RR-MM SCT"] <-"grey20"
fin_all$col_code2[fin_all$col_code2 == "RR-MM no_SCT"] <-"purple"
fin_all$col_code2[fin_all$col_code2 == " "] <- "grey50"
fin_all$col_code2[fin_all$col_code2 == "RR-MM "] <- "grey50"
fin_all<-fin_all[order(fin_all$col_code2, fin_all$mutations),]

name_1<- gsub("Signature.Subs.0","SBS",colnames(sig)[-ncol(sig)])
name_2<-gsub<- gsub("Signature.Subs.","SBS",name_1)
name_3<- gsub<- gsub("MM1","SBS-MM1",name_2)

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/transplanat/")
# pdf("barplot_signatures_figure_1.pdf", width = 15, height = 17)
par(mfrow=c(2,1), mar=c(8,10,5,10),xpd=T)
x<- barplot(as.matrix(t(fin_all[,2:9])), col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
            las=2, cex.names = 0.7, space=rep(0, nrow(fin_all)), border = NA, names.arg = rep(NA, nrow(fin_all)), cex.axis=2.5)

# axis(1, at=x, labels=fin_all$V2, col=fin_all$col_code2, las=2, cex=1)
mtext(side=1, text = fin_all$V2, at=x, line = 0.5,  col=fin_all$col_code2, las=2, cex=1)
legend("topright",legend=name_3,bty="n", pch=15,
       col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
       cex=1.5, pt.cex=1.5, inset=c(-0.12,0.0),x.intersp = 1,y.intersp = 1)
x<- barplot(as.matrix(t(fin_all[,2:9]*fin_all$mutations)), col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
            las=2, cex.names = 0.7, space=rep(0, nrow(fin_all)), border = NA,, names.arg = rep(NA, nrow(fin_all)), cex.axis=2.5)
mtext(side=1, text = fin_all$V2, at=x, line = 0.5,  col=fin_all$col_code2, las=2, cex=1)


legend("topright",legend=name_3,bty="n", pch=15,
       col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
       cex=1.5, pt.cex=1.5, inset=c(-0.12,0.0),x.intersp = 1,y.intersp = 1)

# dev.off()

##########################################################################################################################################
###  run pre-post gain considering only catalogue of signature extracted by MM-FIT for each patient
################################################################################################################################################
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# mol_time_sig_96<- read.delim("96_pre_post_gain_per_samples.txt", sep="\t", stringsAsFactors = F)

consigts.defn <- read.csv("mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(mol_time_sig_96)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total

samples.muts<-samples.muts[,colSums(samples.muts)>100]
samples<- colnames(samples.muts)

consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)

if (is.null(sample.sigt.profs)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  sigt.profs <- list()
  for (i in 1:length(samples)) {
    sigt.profs[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}

sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000


for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  ### split smapleID pre-post gain data to identify the correct sample
  out <- str_split_fixed(( samples[j]),'_',2)
  samples_list<- out[,1]
  ###
  sample.sigts <- unique(sigt.profs[[ samples_list ]])
  sample.sigts <- c("Signature.Subs.01",sample.sigts)
  
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    
    ### test the cos-similarities of each signature but SBS1 and SBS5
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    # if(!is.na(rem.names))
    names(cosReduction) <- rem.names
    if (min(cosReduction) < 0.01) {
      spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
      rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
      rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      reducing = TRUE
    }else {spit(dbg, "exiting while...")
      reducing = FALSE
    }
  }
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig_gain <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)


### multiple gain intermediate mut collpased with the pre-gain
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# sig_gain<- read.delim("all_patients_pre_gain_post_gain_contribution.txt", sep="\t", stringsAsFactors = F) ### multiple gain intermediate mut collpased with the pre-gain

head(sig_gain)
sig_gain[sig_gain$Signature.Subs.09 !=0,][,-1]

par(mfrow=c(1,1))
barplot(as.matrix(t(sig_gain[,-c( ncol(sig_gain))])), col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
        las=2, cex.names = 0.7, names.arg = sig_gain$sampleID, space=rep(0, nrow(sig_gain)), border = NA)
legend("topright",legend=(colnames(sig_gain)[-c( ncol(sig_gain))]),bty="n", pch=15, 
       col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
       cex=1, pt.cex=1, inset=c(-0.1,0.0),x.intersp = 1,y.intersp = 1)

sig2<- sig_gain
sig2$sample<- rownames(sig2)  
out <- str_split_fixed((sig2$sample),'_',2) 
sig2$sampleID<- out[,1]
sig2$code<- out[,2]
sig2$code[sig2$code == 1]<- 5
sig2$code[sig2$code == 2]<- 1
sig2$code[sig2$code == "subclonal"]<- 3
sig2$code[sig2$code == 5]<- 2
sig2$code<- as.numeric(as.character(sig2$code))
sam<- unique(sig2$sampleID)
par(mfrow=c(3,3))
signature_list<- colnames(sig2)[1:8]

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_signatures/")
for(j in (1:8))
{
  pdf(paste0("mol_time_sig",j,".pdf"), width = 7.5, height=10)
  par(mar=c(10,8,5,5))
  for(i in (1:length(sam)))
  {
    sig2_sam<- sig2[sig2$sampleID == sam[i],]
    
    plot(sig2_sam$code, sig2_sam[,j], col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3"))[j],
         pch=16, ylim=c(0,1), xlim=c(0,3.5), xaxt="n", xlab="", ylab="", yaxt="n", bty="n")
    
    par(new=T)
    if(nrow(sig2_sam)>1){
      for(i in (1:(nrow(sig2_sam)-1)))
      {
        segments(sig2_sam$code[i],  sig2_sam[i,j], sig2_sam$code[i+1],  sig2_sam[i+1,j],
                 col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3"))[j], lty = 2)
        par(new=T)
      }
      par(new=T)
    }
  }
  axis(1, at=c(1:3), label=c("pre-gain","post-gain","subclonal"), cex.axis=2, las=2)
  axis(2, at=seq(0,1,by=0.2), label=seq(0,1,by=0.2), las=2, cex.axis=2.5)
  title(signature_list[j])
  print(signature_list[j])
  print(pairwise.wilcox.test(sig2[,j],sig2$code))
  print(summary(lm(sig2[,j]~sig2$code)))
  par(new=F)
  dev.off()
}



cluster_clin<-sig2[,10:12]
rownames(cluster_clin)<- cluster_clin$sample
cluster_clin$code[cluster_clin$code == 1]<-"pre-gain"
cluster_clin$code[cluster_clin$code == 2]<-"post-gain"
cluster_clin$code[cluster_clin$code == 3]<-"subclonal"
vec<- cluster_clin$code
names(vec)<- rownames(cluster_clin)
ann_colors = list(vec=c( "pre-gain"="red","post-gain"="orange","subclonal"="green"))
sig2<- sig2[order(sig2$code),]
colnames(sig2)<- gsub("Signature.Subs.0","SBS", colnames(sig2))
colnames(sig2)<- gsub("Signature.Subs.","SBS", colnames(sig2))
colnames(sig2)<- gsub("MM1","SBS-MM1", colnames(sig2))

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_signatures/")
pdf("heatmap_pre_post.pdf", width = 10, height = 8)
library(pheatmap)
pheatmap((as.matrix(t(sig2[,1:8]))), annotation_col=as.data.frame(vec), annotation_colors = ann_colors, show_colnames = F, 
         annotation_names_col = F, space=rep(0,nrow(sig2) ), border_color = F, fontsize_row = 20)
dev.off()


#### explain that MM1 was detected in clonal fraction without mulktiple samples where likely the trunk don't reflect the real trunk


################################################################################################################################################
################################################################################################################################################
###
###  run DP  considering only catalogue of signature extracted by MM-FIT for each patient
###
################################################################################################################################################
################################################################################################################################################
# 
consigts.defn <- read.csv("mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(DP)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total

samples.muts<-samples.muts[,colSums(samples.muts)>100]
samples<- colnames(samples.muts)

consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)

if (is.null(sample.sigt.profs)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  sigt.profs <- list()
  for (i in 1:length(samples)) {
    sigt.profs[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}

sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000

for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  ### split smapleID pre-post gain data to identify the correct sample
  out <- str_split_fixed(( samples[j]),'_',2)
  samples_list<- out[,1]
  ###
  sample.sigts <- unique(sigt.profs[[ samples_list ]])
  sample.sigts <- c("Signature.Subs.01",sample.sigts)
  
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    
    ### test the cos-similarities of each signature but SBS1 and SBS5
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    # if(!is.na(rem.names))
    names(cosReduction) <- rem.names
    if (min(cosReduction) < 0.01) {
      spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
      rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
      rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      reducing = TRUE
    }else {spit(dbg, "exiting while...")
      reducing = FALSE
    }
  }
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig_DP <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)


# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(sig, "all_patients_DP_contribution.txt", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(sig, "all_patients_DP_contribution.txt", sep="\t", quote=F)
# 
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# sig_DP<- read.delim("all_patients_DP_contribution.txt", sep="\t", stringsAsFactors = F)

example2<- sig_DP[grep("PD26403", rownames(sig)),]
colnames(example2)<- gsub("Signature.Subs.0","SBS",colnames(example2))
colnames(example2)<- gsub("Signature.Subs.","SBS",colnames(example2))
colnames(example2)<- gsub("MM1","SBS-MM1",colnames(example2))


#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/transplanat/")
for(i in (1:nrow(example2)))
{
  pdf(paste0("cluster",i,".pdf"), width = 5, height = 5 )
  par(mar=c(4,4,4,4))
  pie(as.numeric(as.character(example2[i,c(1,3:5,8)])), col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3"))[c(1,3:5,8)],
      cex=1, labels = colnames(example2)[c(1,3:5,8)], radius = 0.9, main=example2$mutations[i])
  dev.off()
}

par(mfrow=c(1,1))
barplot(as.matrix(t(sig_DP[,-c( ncol(sig_DP))])), col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
        las=2, cex.names = 0.7, names.arg = sig_DP$sampleID, space=rep(0, nrow(sig_DP)), border = NA)
legend("topright",legend=(colnames(sig_DP)[-c( ncol(sig_DP))]),bty="n", pch=15, 
       col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
       cex=1, pt.cex=1, inset=c(-0.1,0.0),x.intersp = 1,y.intersp = 1)


pheatmap(as.matrix(t(sig_DP[,1:8]))) ### with all samples

#### select only cases with multiple samples

multiple_samples<- c("PD26400", "PD26401","PD26402","PD26403",  "PD26404", "PD26405",
                     "PD26406", "PD26407","PD26408","PD26409",  "PD26411", "PD26412",
                     "PD26414", "PD26415","PD26416","PD26418",  "PD26419", "PD26420",
                     "PD26422", "PD26423","PD264024","PD26425",  "PD26427", "PD26428",
                     "PD26432", "PD26435")

#setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
cluster_clin<- read.delim("sample_names_DIRIC.txt", sep="\t", stringsAsFactors = F, header=F)

cluster_clin$V1<- gsub("PD26406_06","PD26406_6", cluster_clin$V1)

count_kk<- sig_DP[grep("PD", rownames(sig_DP)),1:9]
count_kk$V1<- rownames(count_kk)
com_dp_count<- merge(count_kk, cluster_clin, by="V1")



kk<- sig_DP[grep("PD", rownames(sig_DP)),1:8]
kk2<- kk[-grep("NA", rownames(kk)),]
kk2<- kk2[-grep("26429", rownames(kk2)),]
kk2<- kk2[-grep("26410", rownames(kk2)),]
kk2<- kk2
kk2$V1<- rownames(kk2)
com_dp<- merge(kk2, cluster_clin, by="V1")

out <- str_split_fixed((com_dp$V1),'_',2) 
rownames(com_dp)<- (com_dp$V1)
vec<- com_dp$V2
names(vec)<- rownames(com_dp)
ann_colors = list(vec=c( "trunk"="red","branches"="green"))
pheatmap(as.matrix( t(kk2[,-ncol(kk2)])), annotation_col=as.data.frame(vec), annotation_colors = ann_colors)

################################################
####
#### check AID in the branches
####
################################################

head(cluster_clin)

#setwd("~/Desktop/clock/2019_dirichlet_dl8/")
code_file<- read.delim("code_files.txt", sep="\t", stringsAsFactors = F)

setwd("/Users/yellapav/Desktop/record/clock_paper/data/files/optimal_info")
ff<-list.files()

chap_list<- list()
for(i in (1:length(code_file)))
{
  tree<- read.delim(ff[i], stringsAsFactors = F, sep="\t")
  num<- code_file[code_file$code == gsub("_optimaInfo.txt","",ff[i]),]
  tree$sample<- num$sample
  tree$clon_code<- tree$location
  tree$clon_code[tree$clon_code>=0.8]<-"clonal"
  tree$clon_code[tree$clon_code<0.8]<-"subclonal"
  chap_list[[i]]<- tree[,c("sample", "no.of.mutations","cluster.no","clon_code")]
}
chap_list2<- do.call("rbind", chap_list)

setwd("~/Desktop/record/clock_paper/optimal_info_kevin//")
ffk<-list.files()[grep("txt", list.files())][-c(11,24,29,27 )]
kevin_list<- list()
for(i in (1:length(ffk)))
{
  tree<- read.delim(ffk[i], stringsAsFactors = F, sep="\t")
  
  tree$sample<- gsub("_optimaInfo_0.01.txt","", ffk[i])
  tree$clon_code<- NA
  columns_num<- grep("PD", colnames(tree))
  for(j in (1:length(columns_num)))
  {
    tree$clon_code[tree[,columns_num[j]]>=0.8]<-"clonal"
  }
  tree$clon_code[is.na( tree$clon_code)]<-"subclonal"
  tree$cluster.no<- rownames(tree)
  kevin_list[[i]]<- tree[,c("sample", "no.of.mutations.assigned","cluster.no","clon_code")]
}
kevin_list2<- do.call("rbind", kevin_list)


#setwd("~/Desktop/record/clock_paper/optimal_info_kevin/")
ffk_single<-list.files()[grep("txt", list.files())][c(11,24,29,27 )]
kevin_list_single<- list()

for(i in (1:length(ffk_single)))
{
  tree<- read.delim(ffk_single[i], stringsAsFactors = F, sep="\t")
  
  tree$sample<- gsub("_optimaInfo_0.01.txt","", ffk_single[i])
  tree$clon_code<- NA
  tree$clon_code[tree[,columns_num[j]]>=0.8]<-"clonal"
  tree$clon_code[is.na( tree$clon_code)]<-"subclonal"
  tree$cluster.no<- rownames(tree)
  kevin_list_single[[i]]<- tree[,c("sample", "no.of.mutations","cluster.no","clon_code")]
}
kevin_list_single2<- do.call("rbind", kevin_list_single)

colnames(chap_list2)<- colnames(kevin_list2)
colnames(kevin_list_single2)<- colnames(kevin_list2)
all_trees<- rbind.data.frame(kevin_list2, kevin_list_single2, chap_list2)
all_trees2<- all_trees[all_trees$no.of.mutations.assigned >100, ]
all_trees2$sampleID<- paste(all_trees2$sample, all_trees2$cluster.no, sep="_")

# colnames(com_dp)[11]<- "paz"
colnames(com_dp)[1]<- "sampleID"
def_dp<- merge(com_dp, all_trees2, by="sampleID")

aid_branch<- def_dp[def_dp$Signature.Subs.09!=0 & def_dp$V2=="branches",]


mm_col <- c(rep("lightskyblue",16), rep("black",16), rep("firebrick2",16),rep("gray88",16),rep("darkolivegreen3",16),
            rep("lightpink1",16))

aid_branc<- com_dp_count[com_dp_count$Signature.Subs.09!=0 & com_dp_count$V2=="branches",]

rownames(def_dp)<- def_dp$sampleID
vec<- paste(def_dp$V2, def_dp$clon_code)
names(vec)<- paste(def_dp$sampleID)
ann_colors = list(vec=c( "trunk clonal"="forestgreen","branches clonal"="darkolivegreen2","branches subclonal"="goldenrod1"))

def_dp_h<- def_dp
colnames(def_dp_h)<- gsub("Signature.Subs.0", "SBS", colnames(def_dp_h))
colnames(def_dp_h)<- gsub("Signature.Subs.", "SBS", colnames(def_dp_h))
colnames(def_dp_h)[colnames(def_dp_h)=="MM1"]<-"SBS-MM1"

setwd("/Users/yellapav/Desktop/record/clock_paper/figures/")
pdf("Heat_map_DP.pdf", width = 20, height = 10)
library(pheatmap)
pheatmap(as.matrix( t(def_dp_h[,2:9])), annotation_col=as.data.frame(vec), annotation_colors = ann_colors, 
         show_colnames =  F, fontsize_row = 20, annotation_legend =F, space=rep(0, nrow(def_dp_h)),border_color =F)
dev.off()

# pdf("Heat_map_DP_Legend.pdf", width = 20, height = 10)
# plot(1,1)
# legend("topright", bty = "n", col=c("forestgreen", "darkolivegreen2", "goldenrod1"),
# legend=c("Early Clonal","Late Clonal", "Subclonal"), pch=15, cex=4)
# dev.off()

kk<- def_dp[def_dp$Signature.Subs.09!=0 & def_dp$V2=="branches",]
unique(kk$sample)

def_dp_h$vio<- paste(def_dp_h$V2, def_dp_h$clon_code, sep=" - ")


def_dp_h$vio[def_dp_h$vio == "trunk - clonal" ]<- "Early Clonal"
def_dp_h$vio[def_dp_h$vio == "branches - subclonal"]<- "Subclonal"
def_dp_h$vio[def_dp_h$vio == "branches - clonal"]<-"Late Clonal"


# median_order<- median_order[order(median_order$x, decreasing=T),]
palette<- brewer.pal(n = 8, name = "RdBu")
county_list<- colnames(def_dp_h)[2:9]
setwd("/Users/yellapav/Desktop/record/clock_paper/figures/")
for (i in c(1:length((county_list)))) { 
  
  def_dp2<- def_dp_h[,c(i+1,ncol(def_dp_h))]  
  colnames(def_dp2)<-c("var","cat")
  median_order<- as.data.frame(aggregate(def_dp_h[, i+1], list(def_dp_h$vio), median))
  median_order$Group.1<- as.character(median_order$Group.1)
  
  par(mar=c(8,8,2,2))
  p<- ggplot(def_dp2, aes(x=as.factor(cat), y=as.numeric(var), fill=cat)) +
    geom_violin(trim=FALSE)+ guides(fill=FALSE) + 
    scale_fill_manual(values=c("forestgreen", "darkolivegreen2", "goldenrod1")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.text=element_text(size=20)) +
    ylim(0,1) +  geom_jitter(shape=16, position=position_jitter(0.2), size=2) + 
    stat_summary(fun.y=median, geom="point", size=4, color="brown2") +
    scale_x_discrete(limits=(median_order$Group.1)) +
    xlab("") + ylab("Contribution") +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) + 
    ggtitle(colnames(def_dp_h)[i+1]) + 
    theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))
  p      
  ggsave(paste0(colnames(def_dp_h)[i+1],".pdf"), width = 7, height = 10)
}

### trunk contribution and median
sample_list<- unique(def_dp_h$sample)
median_tree=NULL
for(i in (1:length(sample_list)))
{
  def_dp_h_sam<-  def_dp_h[def_dp_h$sample==sample_list[i],]
  median_tree<- rbind(median_tree,c(sample_list[i],def_dp_h_sam$no.of.mutations.assigned[def_dp_h_sam$vio =="Early Clonal"]/sum(def_dp_h_sam$no.of.mutations.assigned)))
}
median_tree2<- as.data.frame.matrix(median_tree)
median_tree2$V2<- as.numeric(as.character(median_tree2$V2))

def_dp_h[def_dp_h$vio =="Early Clonal" & def_dp_h$SBS2 ==0,]
def_dp_h[def_dp_h$vio !="Early Clonal" & def_dp_h$SBS2 ==0,]

def_dp_h[def_dp_h$vio =="Early Clonal" & def_dp_h$SBS9 !=0,]

#### check 96 classes of branches and late clonal with AID
aid_late<- def_dp_h[def_dp_h$vio !="Early Clonal" & def_dp_h$SBS9 ==0,][,c(1:10,14,15)]

aid_late<- def_dp_h[def_dp_h$vio !="Early Clonal" & def_dp_h$SBS9 !=0,][,c(1:10,14,15)]
jj<- as.matrix(rowSums(DP))
jj2<- as.data.frame(jj)
jj2$sampleID <- rownames(jj2)

fin_aid<- merge(aid_late,jj2, by="sampleID")
fin_aid$V1<- as.numeric(as.character(fin_aid$V1))
fin_aid$aid_burden<- fin_aid$SBS9*fin_aid$V1
fin_aid2<- fin_aid[fin_aid$aid_burden >100,]




#setwd("/Users/yellapav/Desktop/record/clock_paper/figures/clock/analysis/dirichlet_kevin/")
pdf("aid_in_branches.pdf", width = 30, height = 20)
par(mfrow=c(4,2), mar=c(5,10,5,5))
DP_aid_lat<- DP[rownames(DP) %in% fin_aid2$sampleID,]
for(i in (1:nrow(DP_aid_lat)))
{
  barplot(as.numeric(DP_aid_lat[i,]), col=mm_col, names.arg ="", las=2, main=rownames(DP_aid_lat)[i], border = F, cex.axis=2.5,
          cex.main=3)
  mtext(side=2, cex=2, text = "Number of mutations", line=5)
}
dev.off()  

def_dp_h[def_dp_h$sample %in%fin_aid2$sample, ][,c(1:10,14,15)]

################################################################################################################################################
################################################################################################################################################
###
###  run only for LOH
###
################################################################################################################################################
# ################################################################################################################################################
# 
consigts.defn <- read.csv("/Users/yellapav/Desktop/record/clock_paper/data/mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(mol_time_sig_LOH)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total

samples.muts<-samples.muts[,colSums(samples.muts)>50]
samples<- colnames(samples.muts)

consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)

if (is.null(sample.sigt.profs)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  sigt.profs <- list()
  for (i in 1:length(samples)) {
    sigt.profs[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}

sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000

for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  ### split smapleID pre-post gain data to identify the correct sample
  out <- str_split_fixed(( samples[j]),'_',2)
  samples_list<- out[,1]
  ###
  if(length(grep("PD", samples_list))>0){
    samples_list<- substr(samples_list,1,nchar(samples_list)-1)
  }
  sample.sigts <- unique(sigt.profs[[ samples_list ]])
  sample.sigts <- c("Signature.Subs.01",sample.sigts)
  
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  # if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    
    ### test the cos-similarities of each signature but SBS1 and SBS5
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    
    if(length(rem.names)>0){
      
      names(cosReduction) <- rem.names
      if (min(cosReduction) < 0.01) {
        spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
        rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
        rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        reducing = TRUE
      }else {spit(dbg, "exiting while...")
        reducing = FALSE
      }
    }else{
      reducing = FALSE
    }
    spit(dbg,"... while exited")
    rem.alpha.names <- names(rem.alpha)
    for (n in 1:length(consigt.names)) {
      if (consigt.names[n] %in% rem.alpha.names) {
        sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
      }
      else {
        sigt.fraction[n,j] <- 0
      }
    }
  }
  
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig_LOH <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)
# 
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(sig, "all_patients_loh_PRE_POST_contribution.txt", sep="\t", quote=F)
# 
# 
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# sig<- read.delim("all_patients_loh_PRE_POST_contribution.txt", sep="\t", stringsAsFactors = F) ### multiple gain intermediate mut collpased with the pre-gain

par(mfrow=c(1,1))
barplot(as.matrix(t(sig_LOH[,-c( ncol(sig_LOH))])), col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
        las=2, cex.names = 0.7, names.arg = sig_LOH$sampleID, space=rep(0, nrow(sig_LOH)), border = NA)
legend("topright",legend=(colnames(sig_LOH)[-c( ncol(sig_LOH))]),bty="n", pch=15, 
       col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
       cex=1, pt.cex=1, inset=c(-0.1,0.0),x.intersp = 1,y.intersp = 1)

sig_LOH2<- sig_LOH[order(rownames(sig_LOH)),]
sig_LOH2[,c("Signature.Subs.09","mutations")]
sig_LOH2$sample<- rownames(sig_LOH2)  
out <- str_split_fixed((sig_LOH2$sample),'_',3) 
sig_LOH2$sampleID<- out[,1]
sig_LOH2$code<- out[,3]
sig_LOH2$code[sig_LOH2$code == 1]<- 5 ### reorder number for plotting
sig_LOH2$code[sig_LOH2$code == 2]<- 1
sig_LOH2$code[sig_LOH2$code == "subclonal"]<- 3
sig_LOH2$code[sig_LOH2$code == 5]<- 2
sig_LOH2$code<- as.numeric(as.character(sig_LOH2$code))
sam<- unique(sig_LOH2$sampleID)
sig_LOH2<- sig_LOH2[sig_LOH2$code!="subclonal",]

sig_LOH2_no_apo<- sig_LOH2[!sig_LOH2$sampleID %in% c("MMRC006BBM","PD26419"),] 

par(mfrow=c(3,3))
signature_list<- colnames(sig_LOH2)[1:8]
for(j in (1:8))
{
  for(i in (1:length(sam)))
  {
    sig_LOH2_sam<- sig_LOH2[sig_LOH2$sampleID == sam[i],]
    
    plot(sig_LOH2_sam$code, sig_LOH2_sam[,j], col="forestgreen", pch=16, ylim=c(0,1), xlim=c(0,3.5), xaxt="n", xlab="", ylab="", yaxt="n", bty="n")
    
    par(new=T)
    if(nrow(sig_LOH2_sam)>1){
      for(i in (1:(nrow(sig_LOH2_sam)-1)))
      {
        segments(sig_LOH2_sam$code[i],  sig_LOH2_sam[i,j], sig_LOH2_sam$code[i+1],  sig_LOH2_sam[i+1,j], col = "forestgreen", lty = 2)
        par(new=T)
      }
      par(new=T)
    }
  }
  axis(1, at=c(1:3), label=c("pre-gain","post-gain","subclonal"))
  axis(2, at=seq(0,1,by=0.2), label=seq(0,1,by=0.2), las=2)
  title(signature_list[j])
  par(new=F)
}


sig_LOH2$code[sig_LOH2$code==1]<-"pre-LOH"
sig_LOH2$code[sig_LOH2$code==2]<-"post-LOH"
sig_LOH2$code<- factor(sig_LOH2$code, levels = c("pre-LOH","post-LOH"))
median_order<- as.data.frame(aggregate(sig_LOH2$Signature.Subs.09, list(sig_LOH2$code), median))
median_order$Group.1<- as.character(median_order$Group.1)
sig_LOH2<- sig_LOH2[order(sig_LOH2$code),]

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_signatures/")
par(mar=c(8,8,2,2))
p<- ggplot(sig_LOH2, aes(x=as.factor(code), y=as.numeric(Signature.Subs.09), fill=code)) +
  geom_violin(trim=FALSE)+ guides(fill=FALSE) + 
  scale_fill_manual(values=c("darkolivegreen2", "forestgreen")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), axis.text=element_text(size=20)) +
  ylim(0,1) +  geom_jitter(shape=16, position=position_jitter(0.2), size=2) + 
  stat_summary(fun.y=median, geom="point", size=4, color="brown2") +
  scale_x_discrete(limits=(median_order$Group.1)) +
  xlab("") + ylab("Contribution") +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) 
p      
# ggsave("LOH_signatures.pdf", width = 7, height = 10)

head(mol_time_sig_LOH)
mol_time_sig_LOH2<- mol_time_sig_LOH
mol_time_sig_LOH2$sample<- rownames(mol_time_sig_LOH2)  
out <- str_split_fixed((mol_time_sig_LOH2$sample),'_',3) 
mol_time_sig_LOH2$sampleID<- out[,1]
mol_time_sig_LOH2$code<- out[,3]

mol_time_sig_LOH2_2<- mol_time_sig_LOH2[mol_time_sig_LOH2$sampleID %in% sig_LOH2$sampleID,]
mol_time_sig_LOH2_2<- mol_time_sig_LOH2_2[!mol_time_sig_LOH2_2$sampleID %in% c("MMRC006BBM","PD26419"),] ### exclude cases without AID to show differential AID contribution
code_list<- unique(mol_time_sig_LOH2_2$code)
tab_96<- list()
for(i in (1:3))
{
  int_file<-  mol_time_sig_LOH2_2[mol_time_sig_LOH2_2$code == code_list[i],]
  tab_96[[i]]<- c(colSums(int_file[,1:96]), unique(int_file$code))
}
tab_92_2<- do.call("rbind", tab_96)
tab_92_2<- as.data.frame.matrix(tab_92_2)
tab_92_2[,1:96]<- apply(tab_92_2[,1:96], 2, function(x){as.numeric(as.character(x))})

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_signatures/")
pdf("96_classes_pre_post.pdf", width = 10, height = 5)
par(mfrow=c(2,1))
barplot(as.numeric(tab_92_2[2,1:96]), col=mm_col, main="pre-LOH", las=2, border = F)
barplot(as.numeric(tab_92_2[1,1:96]), col=mm_col, main="post-LOH", las=2, border = F)
dev.off()


mol_time_sig_LOH2_2<- mol_time_sig_LOH2[mol_time_sig_LOH2$sampleID %in% sig_LOH2$sampleID,]
mol_time_sig_LOH2_2<- mol_time_sig_LOH2_2[mol_time_sig_LOH2_2$sampleID %in% c("MMRC006BBM","PD26419"),] ## include only cases without AID
code_list<- unique(mol_time_sig_LOH2_2$code)
tab_92_2_apo<- list()
for(i in (1:3))
{
  int_file<-  mol_time_sig_LOH2_2[mol_time_sig_LOH2_2$code == code_list[i],]
  tab_92_2_apo[[i]]<- c(colSums(int_file[,1:96]), unique(int_file$code))
}
tab_92_2_apo<- do.call("rbind", tab_96)
tab_92_2_apo<- as.data.frame.matrix(tab_92_2_apo)
tab_92_2_apo[,1:96]<- apply(tab_92_2_apo[,1:96], 2, function(x){as.numeric(as.character(x))})

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/no_aid_just_APOBEC//")
pdf("96_classes_pre_post_onlky_APOBEC.pdf", width = 10, height = 5)
par(mfrow=c(2,1))
barplot(as.numeric(tab_92_2_apo[2,1:96]), col=mm_col, main="pre-LOH", las=2, border = F)
barplot(as.numeric(tab_92_2_apo[1,1:96]), col=mm_col, main="post-LOH", las=2, border = F)
dev.off()


################################################################################################################################################
################################################################################################################################################
###
### Intermediate_gains
###
###
################################################################################################################################################
################################################################################################################################################
# 
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# mol_time_sig_96_int<- read.delim("96_pre_post_gain_per_samples_int_samples.txt", sep="\t", stringsAsFactors = F)

consigts.defn <- read.csv("/Users/yellapav/Desktop/record/clock_paper/data/mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(mol_time_sig_96_int)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total

samples.muts<-samples.muts[,colSums(samples.muts)>50]
samples<- colnames(samples.muts)

consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)

dbg=FALSE
if (is.null(sample.sigt.profs)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  sigt.profs <- list()
  for (i in 1:length(samples)) {
    sigt.profs[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}

sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000


for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  ### split smapleID pre-post gain data to identify the correct sample
  out <- str_split_fixed(( samples[j]),'_',2)
  samples_list<- out[,1]
  ###
  sample.sigts <- unique(sigt.profs[[ samples_list ]])
  sample.sigts <- c("Signature.Subs.01",sample.sigts)
  
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  # if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    
    ### test the cos-similarities of each signature but SBS1 and SBS5
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    
    if(length(rem.names)>0){
      
      names(cosReduction) <- rem.names
      if (min(cosReduction) < 0.01) {
        spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
        rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
        rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        reducing = TRUE
      }else {spit(dbg, "exiting while...")
        reducing = FALSE
      }
    }else{
      reducing = FALSE
    }
    spit(dbg,"... while exited")
    rem.alpha.names <- names(rem.alpha)
    for (n in 1:length(consigt.names)) {
      if (consigt.names[n] %in% rem.alpha.names) {
        sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
      }
      else {
        sigt.fraction[n,j] <- 0
      }
    }
  }
  
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig_int <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(sig, "all_patients_INT_PRE_POST_contribution.txt", sep="\t", quote=F)

# 
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# sig_int<- read.delim("all_patients_INT_PRE_POST_contribution.txt", sep="\t", stringsAsFactors = F) ### multiple gain intermediate mut collpased with the pre-gain
sig_int$sample<- rownames(sig_int)  
out <- str_split_fixed((sig_int$sample),'_',2) 
sig_int$sampleID<- out[,1]
sig_int$code<- out[,2]
sam_multi<- sig_int[sig_int$code==3,]


sig_int$code[sig_int$code=="subclonal"]<-5
sig_int$code[sig_int$code=="1"]<-4
sig_int$code[sig_int$code=="2"]<-2
sig_int$code[sig_int$code=="3"]<-1
sig_int$code[sig_int$code=="4"]<-3
sig_int$code[sig_int$code=="5"]<-4
sig_int$code<- as.numeric(as.character(sig_int$code))

sig_int_int<- sig_int[sig_int$sampleID %in% sam_multi$sampleID,]
sig_int_int[,1:8]<- apply(sig_int_int[,1:8],2, function(x){as.numeric(as.character(x))})

sam<- unique(sig_int_int$sampleID)
par(mfrow=c(3,3))
signature_list<- colnames(sig_int_int)[1:8]
for(j in (1:8))
{
  for(i in (1:length(sam)))
  {
    sig2_sam<- sig_int_int[sig_int_int$sampleID == sam[i],]
    
    plot(sig2_sam$code, sig2_sam[,j], col="forestgreen", pch=16, ylim=c(0,1), xlim=c(0,4.5), xaxt="n", xlab="", ylab="", yaxt="n", bty="n")
    
    par(new=T)
    sig2_sam<-sig2_sam[order(sig2_sam$code),]
    for(i in (1:(nrow(sig2_sam)-1)))
    {
      segments(sig2_sam$code[i],  sig2_sam[i,j], sig2_sam$code[i+1],  sig2_sam[i+1,j], col = "forestgreen", lty = 2)
      par(new=T)
    }
    par(new=T)
  }
  axis(1, at=c(1:4), label=c("pre-1st gain","pre-2nd gain","post-gain","subclonal"))
  axis(2, at=seq(0,1,by=0.2), label=seq(0,1,by=0.2), las=2)
  title(signature_list[j])
  par(new=F)
}

mol_time_sig_INT<-mol_time_sig_96_int
mol_time_sig_INT$sample<- rownames(mol_time_sig_INT)  
out <- str_split_fixed((mol_time_sig_INT$sample),'_',2) 
mol_time_sig_INT$sampleID<- out[,1]
mol_time_sig_INT$code<- out[,2]

mol_time_sig_INT_2<- mol_time_sig_INT[mol_time_sig_INT$sampleID %in% sam,]
mol_time_sig_INT_2<- mol_time_sig_INT_2[!mol_time_sig_INT_2$sampleID %in% c("MMRC006BBM","PD26419"),]
code_list<- unique(mol_time_sig_INT_2$code)
tab_96_int<- list()
for(i in (1:4))
{
  int_file<-  mol_time_sig_INT_2[mol_time_sig_INT_2$code == code_list[i],]
  tab_96_int[[i]]<- c(colSums(int_file[,1:96]), unique(int_file$code))
}
tab_92_2_int<- do.call("rbind", tab_96_int)
tab_92_2_int<- as.data.frame.matrix(tab_92_2_int)
tab_92_2_int[,1:96]<- apply(tab_92_2_int[,1:96], 2, function(x){as.numeric(as.character(x))})

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_signatures/")
pdf("INT_mutations_classes_pre_post.pdf", width = 10, height = 3)
par(mfrow=c(1,1))
barplot(as.numeric(tab_92_2_int[4,1:96]), col=mm_col, main="Intermediate", las=2, border = F)
# barplot(as.numeric(tab_92_2[3,1:96]), col=mm_col, main="post-LOH", las=2)
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### figure paper for inter-gains and LOH [XX]
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/mol_time_signatures/")
#pdf("96_classes_pre_post.pdf", width = 10, height = 10)
par(mfrow=c(3,1))
barplot(as.numeric(tab_92_2_loh[2,1:96]), col=mm_col, main="pre-LOH", las=2, border = F)
barplot(as.numeric(tab_92_2_loh[1,1:96]), col=mm_col, main="post-LOH", las=2, border = F)
barplot(as.numeric(tab_92_2_int[4,1:96]), col=mm_col, main="Intermediate", las=2, border = F)
# barplot(as.numeric(tab_92_2[3,1:96]), col=mm_col, main="post-LOH", las=2)
#dev.off()





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### figure supplementary just with APOBEC high cases
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
setwd("/Users/yellapav/Desktop/record/clock_paper/figures/")
mol_time_sig_INT<-mol_time_sig_96_int
mol_time_sig_INT$sample<- rownames(mol_time_sig_INT)  
out <- str_split_fixed((mol_time_sig_INT$sample),'_',2) 
mol_time_sig_INT$sampleID<- out[,1]
mol_time_sig_INT$code<- out[,2]

mol_time_sig_INT_2<- mol_time_sig_INT[mol_time_sig_INT$sampleID %in% sam,]
mol_time_sig_INT_2<- mol_time_sig_INT_2[mol_time_sig_INT_2$sampleID %in% c("MMRC006BBM","PD26419"),]
code_list<- unique(mol_time_sig_INT_2$code)
tab_96<- list()
for(i in (1:4))
{
  int_file<-  mol_time_sig_INT_2[mol_time_sig_INT_2$code == code_list[i],]
  tab_96[[i]]<- c(colSums(int_file[,1:96]), unique(int_file$code))
}
tab_92_2<- do.call("rbind", tab_96)
tab_92_2<- as.data.frame.matrix(tab_92_2)
tab_92_2[,1:96]<- apply(tab_92_2[,1:96], 2, function(x){as.numeric(as.character(x))})

###

########### figure paper

#setwd("/Volumes/GoogleDrive/My Drive/clock/analysis/no_aid_just_APOBEC/")
pdf("96_classes_pre_post.pdf", width = 10, height = 10)
par(mfrow=c(3,1))
barplot(as.numeric(tab_92_2[2,1:96]), col=mm_col, main="pre-LOH", las=2, border = F)
barplot(as.numeric(tab_92_2[1,1:96]), col=mm_col, main="post-LOH", las=2, border = F)
barplot(as.numeric(tab_92_2_int[4,1:96]), col=mm_col, main="int", las=2, border = F)
# barplot(as.numeric(tab_92_2[3,1:96]), col=mm_col, main="post-LOH", las=2)
dev.off()


################################################################################################################################################
################################################################################################################################################
###
###  run pre-post gain for each chromosome gain - ANALYSIS FOR CLOCK
###
###
################################################################################################################################################
################################################################################################################################################
# 
consigts.defn <- read.csv("/Users/yellapav/Desktop/record/clock_paper/data/mm_signature_definitions.csv")
consigts.defn<- consigts.defn[,-c(1,10)]
tcons.defn <- t(consigts.defn[,3:dim(consigts.defn)[2]])# 

mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),
                                        substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
muttype = rep(1:96,2)
mutlist_oppstrand = sapply(mutlist, oppstrand)
muttype = rep(1:96,2)
names(muttype) = c(mutlist,mutlist_oppstrand)

samples.muts <- as.data.frame(t(mol_time_single_chomr)) # rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total
#
# samples.muts<-samples.muts[,colSums(samples.muts)>100]
samples<- colnames(samples.muts)

consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric )  #n be cautious in the original the first two columns are included
consigt.names <- colnames(consigts.defn)

if (is.null(sample.sigt.profs)) {
  spit(dbg, "using mm prior signature profiles")
  mm.sigts <- sort(c("Signature.Subs.01","Signature.Subs.02","Signature.Subs.05","Signature.Subs.08","Signature.Subs.09",
                     "Signature.Subs.13","Signature.Subs.18","MM1"))
  
  
  ##### adjust list of signature for each case
  
  sigt.profs <- list()
  for (i in 1:length(samples)) {
    sigt.profs[[samples[i]]] <- mm.sigts
  }
} else if (is.character(sample.sigt.profs)) {
  library("R.filesets")            # loadRDS()
  spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
  sigt.profs <- loadRDS(sample.sigt.profs)
} else {
  spit(dbg, "using mm signature profiles from input argument")
  sigt.profs <- sample.sigt.profs
}



sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
rownames(sigt.fraction) <- consigt.names
colnames(sigt.fraction) <- samples
max.em.iter=2000

for (j in 1:length(samples)) {
  sample.mut.freqs = as.numeric(samples.muts[,j])
  sample.mut.freqs[is.na(sample.mut.freqs)] = 0
  
  ### split smapleID pre-post gain data to identify the correct sample
  out <- str_split_fixed(( samples[j]),'_',2)
  samples_list<- out[,1]
  ###
  sample.sigts <- unique(sigt.profs[[ samples_list ]])
  sample.sigts <- c("Signature.Subs.01",sample.sigts)
  
  sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
  sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
  
  spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
  
  alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
  spat(dbg, "alpha", alpha)
  sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
  sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]  #nicos: this aint doing nothing... sampleAlpha = alpha
  
  # if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
  spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
  reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs) ######### i think this is teh part where we should perfor it...
  sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
  spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
  
  rem.alpha <- sampleAlpha                     # this will hold the final result
  rem.sample.consigts.defn <- sample.consigts.defn
  spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
  reducing = TRUE
  while (reducing) {
    spat(dbg, "in the while, rem.alpha: ", rem.alpha)
    cosReduction <- NULL
    
    ### test the cos-similarities of each signature but SBS1 and SBS5
    rem.names <- setdiff(names(rem.alpha),c("Signature.Subs.01","Signature.Subs.05"))
    for(c in rem.names){
      spit(dbg, "doing c: %s", c)
      red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
      red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
      red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
      red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
      cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
    }
    
    if(length(rem.names)>0){
      
      names(cosReduction) <- rem.names
      if (min(cosReduction) < 0.01) {
        spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
        rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
        rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        reducing = TRUE
      }else {spit(dbg, "exiting while...")
        reducing = FALSE
      }
    }else{
      reducing = FALSE
    }
    spit(dbg,"... while exited")
    rem.alpha.names <- names(rem.alpha)
    for (n in 1:length(consigt.names)) {
      if (consigt.names[n] %in% rem.alpha.names) {
        sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
      }
      else {
        sigt.fraction[n,j] <- 0
      }
    }
  }
  
  spit(dbg,"... while exited")
  rem.alpha.names <- names(rem.alpha)
  for (n in 1:length(consigt.names)) {
    if (consigt.names[n] %in% rem.alpha.names) {
      sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
    }
    else {
      sigt.fraction[n,j] <- 0
    }
  }
}

spat(dbg, "sigt.fraction", sigt.fraction)
tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
colsums.samples.muts <- colSums(samples.muts)
sig_clock <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)


# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# write.table(sig, "all_patients_pre_gain_single_chrom_contribution.txt", sep="\t", quote=F)
# setwd("/Volumes/GoogleDrive/My Drive/clock/data/")
# write.table(sig, "all_patients_pre_gain_single_chrom_contribution.txt", sep="\t", quote=F)
# 
# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# sig<- read.delim("all_patients_pre_gain_single_chrom_contribution.txt", sep="\t", stringsAsFactors = F)

# par(mfrow=c(1,1))
# barplot(as.matrix(t(sig[,-c( ncol(sig))])), col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
#         las=2, cex.names = 0.7, names.arg = sig$sampleID, space=rep(0, nrow(sig)), border = NA)
# legend("topright",legend=(colnames(sig)[-c( ncol(sig))]),bty="n", pch=15,
#        col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
#        cex=1, pt.cex=1, inset=c(-0.1,0.0),x.intersp = 1,y.intersp = 1)
# 
# pheatmap(as.matrix(t(sig[,1:8])))

head(sig_clock)
sig_clock$case<- gsub("-","_", rownames(sig_clock))
out <- str_split_fixed((sig_clock$case),'_',7) 
sig_clock$pre_post<- out[,2]
# sig_clock$order_time_window<- out[,4]
sig_clock$cnv_code<- out[,7]
sig_clock$sample<- out[,1]
sig_clock$chr<- out[,3]
sig_clock$start<- as.numeric(as.character(out[,5]))
sig_clock$end<- as.numeric(as.character(out[,6]))
sig_clock$Signature.Subs.05_abs<- sig_clock$Signature.Subs.05 * sig_clock$mutations ################# select signature 5 for clock
sig_clock$sampleID<-paste(sig_clock$sample, sig_clock$chr, sig_clock$start, sig_clock$end, sig_clock$cnv_code,sep="_")

# which(code_list == "MMRC00GABM_chr15_20001774_102431166_gain_2")

sig_clock$cnv_code[sig_clock$cnv_code == "134944160_gain"]<-"gain"

#######################################################################################################################################################################################
########################################################################################################################################################################################
###
### the following chromosomes need manual curation beacuse mclust splitted in 2 groups the "clonal - cn1" variants. 
### MANUAL CURATION = When you have single gain and 3 clusters, collapse the first 2
###
########################################################################################################################################################################################
########################################################################################################################################################################################

final_clock__mol_time<- NULL
code_list<- sort(unique(sig_clock$sampleID ))
for(i in (1:length(code_list)))
{
  sig_list<- sig_clock[sig_clock$sampleID == code_list[i],]
  
  if(unique(sig_list$cnv_code) == "gain")
  {
    if(nrow(sig_list)==2){
      first<- sig_list[sig_list$pre_post== "1",]
      cn1_clon<- first$Signature.Subs.05_abs
      second<- sig_list[sig_list$pre_post== "2",]
      cn2<- second$Signature.Subs.05_abs
      cn3<- 0
      cn4<- 0
      tot<- sum(sig_list$Signature.Subs.05_abs)
      segment <- unique(sig_list$sampleID)
      cyto<- unique(sig_list$cnv_code)
      chr<- unique(sig_list$chr)
      chrom<- gsub("chr","", chr)
    }else{
      first<- sig_list[sig_list$pre_post== "1" | sig_list$pre_post== "2",]
      if(sig_list$pre_post== "1"){
        cn1_clon<- sum(first$Signature.Subs.05_abs)
        second<- sig_list[sig_list$pre_post== "3",]
        cn2<- 0
        cn3<- 0
        cn4<- 0
        tot<- sum(sig_list$Signature.Subs.05_abs)
        segment <- unique(sig_list$sampleID)
        cyto<- unique(sig_list$cnv_code)
        chr<- unique(sig_list$chr)
        chrom<- gsub("chr","", chr)
      }}
  }else{
    if(unique(sig_list$cnv_code) == "gain_2")
    {
      if(nrow(sig_list) == 2)
      {
        
        jj_code<- sort(sig_list$pre_post)
        first<- sig_list[sig_list$pre_post== jj_code[1],]
        cn1_clon<- first$Signature.Subs.05_abs
        second<- sig_list[sig_list$pre_post==jj_code[2],]
        cn2<- second$Signature.Subs.05_abs
        cn3<- 0
        cn4<- 0
        tot<- sum(sig_list$Signature.Subs.05_abs)
        segment <- unique(sig_list$sampleID)
        cyto<- unique(sig_list$cnv_code)
        chr<- unique(sig_list$chr)
        chrom<- gsub("chr","", chr) 
      }else{
        jj_code<- sort(sig_list$pre_post)
        first<- sig_list[sig_list$pre_post== jj_code[1],]
        cn1_clon<- first$Signature.Subs.05_abs
        second<- sig_list[sig_list$pre_post==jj_code[2],]
        cn2<- second$Signature.Subs.05_abs
        third<-  sig_list[sig_list$pre_post== jj_code[3],]
        cn3<- third$Signature.Subs.05_abs
        cn4<- 0
        tot<- sum(sig_list$Signature.Subs.05_abs)
        segment <- unique(sig_list$sampleID)
        cyto<- unique(sig_list$cnv_code)
        chr<- unique(sig_list$chr)
        chrom<- gsub("chr","", chr) 
      }
    }else{
      jj_code<- sort(sig_list$pre_post)
      first<- sig_list[sig_list$pre_post== jj_code[1],]
      cn1_clon<- first$Signature.Subs.05_abs
      second<- sig_list[sig_list$pre_post== jj_code[2],]
      cn2<- second$Signature.Subs.05_abs
      cn3<- 0
      cn4<- 0
      tot<- sum(sig_list$Signature.Subs.05_abs)
      segment <- unique(sig_list$sampleID)
      cyto<- unique(sig_list$cnv_code)
      chr<- unique(sig_list$chr)
      chrom<- gsub("chr","", chr)
      
    }
  }
  line_mol_time_5<-  c(cn1_clon, cn2, cn3, cn4, tot, segment, cyto, chr, chrom)
  print(paste(i,length(line_mol_time_5) ))
  final_clock__mol_time<- rbind(final_clock__mol_time,as.matrix(t(line_mol_time_5 )))
}


final_clock__mol_time2<- as.data.frame(final_clock__mol_time)
rownames(final_clock__mol_time2)<-   code_list
colnames(final_clock__mol_time2)<-  c("cn1_clon", "cn2", "cn3", "cn4", "tot", "segment", "cyto", "chr", "chrom")
final_clock__mol_time2[,1:5]<- apply(final_clock__mol_time2[,1:5], 2, function(x){as.numeric(as.character(x))})

final_clock__mol_time2$segment<- gsub("_","-",final_clock__mol_time2$segment)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### mol_time for signature 5 only with LOH
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

final_clock__mol_time2_loh<- final_clock__mol_time2[final_clock__mol_time2$cyto == "LOH",]
loh_sig_5<- list()
for(i in (1:nrow(final_clock__mol_time2_loh)))
{
  mol_time_loh<- tetraploid.2.2.or.UPD.2.0.est(as.integer(final_clock__mol_time2_loh$cn1_clon)[i], as.integer(final_clock__mol_time2_loh$cn2)[i], iter=iter)
  pc_ic_row <- mol_time_loh$pests
  loh_sig_5[[i]]<- c(as.numeric(as.character(final_clock__mol_time2_loh[i,1:5])), as.character((final_clock__mol_time2_loh[i,6:8])), pc_ic_row)
}

loh_sig_5_df<- do.call("rbind", loh_sig_5)
loh_sig_5_df2<- as.data.frame(loh_sig_5_df)
colnames(loh_sig_5_df2)[1:9]<- c(colnames(final_clock__mol_time2)[1:8], "mol_time")


loh_sig_5_df2[,c(1:5,9:11)]<- apply(loh_sig_5_df2[,c(1:5,9:11)], 2, function(x) {as.numeric(as.character(x))})
loh_sig_5_df2_chr_big<- loh_sig_5_df2[loh_sig_5_df2$cn1_clon + loh_sig_5_df2$cn2 > 50 & loh_sig_5_df2$cn2 > 5,]

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### mol_time for signature 5 only with trisomy - gain
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

final_clock__mol_time2_gain<- final_clock__mol_time2[final_clock__mol_time2$cyto == "gain",]
gain_sig_5<- list()
for(i in (1:nrow(final_clock__mol_time2_gain)))
{
  mol_time_gain<- triploid.2.1.est(as.integer(final_clock__mol_time2_gain$cn1_clon)[i], as.integer(final_clock__mol_time2_gain$cn2)[i], iter=iter)
  pc_ic_row <- mol_time_gain$pests
  gain_sig_5[[i]]<- c(as.numeric(as.character(final_clock__mol_time2_gain[i,1:5])), as.character(final_clock__mol_time2_gain[i,6:8]), pc_ic_row)
}

gain_sig_5_df<- do.call("rbind", gain_sig_5)
gain_sig_5_df2<- as.data.frame(gain_sig_5_df)
colnames(gain_sig_5_df2)[1:9]<- c(colnames(final_clock__mol_time2)[1:8], "mol_time")



gain_sig_5_df2[,c(1:5,9:11)]<- apply(gain_sig_5_df2[,c(1:5,9:11)], 2, function(x) {as.numeric(as.character(x))})
gain_sig_5_df2_chr_big<- gain_sig_5_df2[gain_sig_5_df2$cn1_clon + gain_sig_5_df2$cn2 > 50 & gain_sig_5_df2$cn2 > 5,]

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### mol_time for signature 5 only with tetrasomy - gain x2 - EXCLUDED FROM HERE
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# final_clock__mol_time2_gain_2<- final_clock__mol_time2[final_clock__mol_time2$cyto == "gain_2",]
# gain_2_sig_5<- list()
# for(i in (1:nrow(final_clock__mol_time2_gain_2)))
# {
#   mol_time_gain_2<- tetraploid.3.1.est(as.integer(final_clock__mol_time2_gain_2$cn1_clon)[i], as.integer(final_clock__mol_time2_gain_2$cn2)[i], as.integer(final_clock__mol_time2_gain_2$cn3)[i], iter=iter)
#   pc_ic_row <- mol_time_gain_2$pests
#   gain_2_sig_5[[i]]<- c(as.numeric(as.character(final_clock__mol_time2_gain_2[i,1:5])), as.character(final_clock__mol_time2_gain_2[i,6:8]), pc_ic_row)
# }
# 
# gain_2_sig_5_df<- do.call("rbind", gain_2_sig_5)
# gain_2_sig_5_df2<- as.data.frame(gain_2_sig_5_df)
# colnames(gain_2_sig_5_df2)[1:9]<- c(colnames(final_clock__mol_time2)[1:8], "mol_time")
# 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### create final file
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


alfa<- rbind.data.frame(gain_sig_5_df2_chr_big, loh_sig_5_df2_chr_big)

order_db<- unique(ranges3[,c("code_ti","order","Sample")] )


order_db$Sample[(grep("PD", order_db$Sample))] = substr(order_db$Sample[(grep("PD", order_db$Sample))],1,nchar(order_db$Sample[(grep("PD", order_db$Sample))])-1)
order_db$segment<- paste(order_db$Sample, order_db$code_ti, sep="-")

final<- merge(alfa, order_db, by="segment")

# setwd("~/Desktop/clock/mol_time_18.01.2018/aa_files/")
# tab_mrca<- read.delim("MRCA_multiple_sample.txt", sep="\t", stringsAsFactors = F)

setwd("/Users/yellapav/Desktop/record/clock_paper/data/")
clin<- read.delim("age_all_chap_fm6.txt", sep="\t", stringsAsFactors = F)
first<- clin[clin$code=="first",]
colnames(first)[3]<-"Sample"

#setwd("~/Desktop/clock/")
tab_mrca<- read.delim("MRCA_multiple_sample_updated.txt", sep="\t",stringsAsFactors = F)

tab_mrca<- unique(tab_mrca[,c("ID","Age","MRCA.est","Est.mut.rate.per.Gb.per.year")])
colnames(tab_mrca)[1]<-"Sample"

clock_int<- merge(tab_mrca, final, by="Sample")

clock<- merge(clock_int[,-2], first, by="Sample")
clock<- clock[!is.na(clock$Est.mut.rate.per.Gb.per.year),]
library(stringr)
out <- str_split_fixed((clock$segment),'-',5) 
callable.genome.size<- as.numeric(as.character(out[,4]))- as.numeric(as.character(out[,3]))
clock$size<- callable.genome.size/1000000000
clock$chr<- out[,2]
# number mut per year per callable genome

clock3<- clock[clock$cn1_clon + clock$cn2 >50 & clock$cn1_clon>10 & clock$cn2> 5,] #### mol_time crietria applied also to signature 5 (more tha  50 SNVs)
clock3<- clock3[order(clock3$Sample),]


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### clock time for each chromosome (high concordance with mol_time but not perfect - some chromosomes change time window)
### likely due to slightely different signature landscape and low mutation burden
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

clock2<-clock3
age_gain<- clock2$cn2/(clock2$size* clock2$Est.mut.rate.per.Gb.per.year)

setwd("/Users/yellapav/Desktop/record/clock_paper/figures/")
pdf("timing_each_chr.pdf", width = 20, height = 30)
par(mar=c(10,10,5,5), xpd=T)
plot(age_gain, 1:length(age_gain), pch=20, xlim=c(-10,(max(clock2$Age)+5)), 
     axes=FALSE, xlab="Estimated time Chrom Gain (years)", ylab="", cex=2)
segments(x0 = age_gain*(clock2$`2.5%`/clock2$mol_time), 
         x1 = age_gain*(clock2$`97.5%`/clock2$mol_time),
         y0 = 1:nrow(clock2))
axis(side=1, at=(0:90)*10, labels=(0:90)*10, cex.axis=2)

text(x = rep(0, nrow(clock2)),
     y = 1:nrow(clock2), labels = paste(clock2$segment, clock2$chr), adj = 1, cex=1)
par(new=T)

plot(clock2$Age, 1:nrow(clock2), pch=15, col="brown3",  xlab="", ylab="",
     axes=FALSE,  , xlim=c(-10,(max(clock2$Age)+5)), cex=2)

par(new=T)

plot(clock2$MRCA.est, 1:nrow(clock2), pch=15, col="blue",  xlab="", ylab="",
     axes=FALSE,  , xlim=c(-10,(max(clock2$Age)+5)), cex=2)
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 
### clock time for each time window - this is the test referred to the other code
### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
head(clock2)

clock2$MRCA.est/clock2$Est.mut.rate.per.Gb.per.year

clock2$code_order<- paste(clock2$Sample, clock2$order, sep="_")
code_order<- unique(clock2$code_order)
compress<- list()
for(i in (1:length(code_order)))
{
  clock2_sam<- clock2[clock2$code_order == code_order[i],]
  first<- clock2_sam[1,]
  first$cn1_clon<- sum(clock2_sam$cn1_clon)
  first$cn2<- sum(clock2_sam$cn2)
  first$cn3<- sum(clock2_sam$cn3)
  first$cn4<- sum(clock2_sam$cn4)
  first$tot<- sum(clock2_sam$tot)
  first$size<- sum(clock2_sam$size)
  compress[[i]]<- first
}

compress2<- do.call("rbind", compress)
compress2<- compress2[compress2$Sample !="PD26412",]
age_gain<- compress2$cn2/(compress2$size* compress2$Est.mut.rate.per.Gb.per.year)
compress2$age_gain<- age_gain

pdf("test.pdf", width = 15, height = 15)
par(mar=c(10,10,5,15), xpd=T)
compress2$age_gain<- as.numeric(as.character(compress2$age_gain))
compress2$MRCA.est<- as.numeric(as.character(compress2$MRCA.est))
compress2$Age<- as.numeric(as.character(compress2$Age))
sample_list<- unique(compress2$Sample)
for(i in (1:length(sample_list)))
{
  compress2_paz<- compress2[compress2$Sample == sample_list[i],]
  if(nrow(compress2_paz)==1)
  {
    plot(compress2_paz$age_gain,i, pch=20, xlim=c(-10,(max(compress2$Age)+5)), ylim=c(0,length(sample_list)), col="forestgreen",
         axes=FALSE, xlab="Estimated time Chrom Gain (years)", ylab="", cex=2)
    segments(x0 = compress2_paz$age_gain*(compress2_paz$`2.5%`/compress2_paz$mol_time), 
             x1 = compress2_paz$age_gain*(compress2_paz$`97.5%`/compress2_paz$mol_time),
             y0 = i, col="forestgreen")
    
    par(new=T)
    plot(compress2_paz$Age,i, pch=16, col="brown3",  xlab="", ylab="",ylim=c(0,length(sample_list)),
         axes=FALSE,  , xlim=c(-10,(max(compress2$Age)+5)), cex=2)
    
    par(new=T)
    plot(compress2_paz$MRCA.est, i, pch=16, col="blue",  xlab="", ylab="",ylim=c(0,length(sample_list)),
         axes=FALSE,  , xlim=c(-10,(max(compress2$Age)+5)), cex=2)
  }else{
    
    compress2_paz<- compress2_paz[order(compress2_paz$age_gain),] ### here i ordered for age and not for number of order - number of order means time window. 2 maybe the first (i.e PD26435)
    plot(compress2_paz$age_gain,rep(i, nrow(compress2_paz)), pch=20, xlim=c(-10,(max(compress2$Age)+5)), ylim=c(0,length(sample_list)), 
         axes=FALSE, xlab="Estimated time Chrom Gain (years)", ylab="", cex=2, col=c("forestgreen","darkolivegreen3"))
    segments(x0 = compress2_paz$age_gain[1]*(compress2_paz$`2.5%`[1]/compress2_paz$mol_time[1]), 
             x1 = compress2_paz$age_gain[1]*(compress2_paz$`97.5%`[1]/compress2_paz$mol_time[1]),
             y0 = i, col="forestgreen")
    segments(x0 = compress2_paz$age_gain[2]*(compress2_paz$`2.5%`[2]/compress2_paz$mol_time[2]), 
             x1 = compress2_paz$age_gain[2]*(compress2_paz$`97.5%`[2]/compress2_paz$mol_time[2]),
             y0 = i, col="darkolivegreen4")
    par(new=T)
    plot(unique(compress2_paz$Age),i, pch=16, col="brown3",  xlab="", ylab="",ylim=c(0,length(sample_list)),
         axes=FALSE,  xlim=c(-10,(max(compress2$Age)+5)), cex=2)
    
    par(new=T)
    plot(unique(compress2_paz$MRCA.est), i, pch=16, col="blue",  xlab="", ylab="",ylim=c(0,length(sample_list)),
         axes=FALSE,  xlim=c(-10,(max(compress2$Age)+5)), cex=2)
    
  }
  par(new=T)
}

axis(side=1, at=seq(0,80, by=10), labels=seq(0,80, by=10), cex.axis=2)
text(x = rep(0, length(sample_list)),
     y = 1:length(sample_list), labels = sample_list, adj = 1.5, cex=1)

legend("topright",legend=c("Diagnosis","MRCA","1nd Time Window","2nd Time Window"),bty="n", pch=16, 
       col=c("brown3","blue","forestgreen","darkolivegreen4"),
       cex=2, pt.cex=2, inset=c(-0.30,0),x.intersp = 1,y.intersp = 1)
dev.off()



compress2_first<- compress2[compress2$order==1,]
dim(compress2_first)
dim(compress2_first[compress2_first$age_gain<20,])
dim(compress2_first[compress2_first$age_gain>30,])

compress2_second<- compress2[compress2$order==2,]

compress2_mult<- compress2[compress2$Sample %in% compress2_second$Sample,]
sam_list_gain<- unique(compress2_mult$Sample)
all_second<- NULL
for(i in (1:length(sam_list_gain)))
{
  compress2_mult2<- compress2_mult[compress2_mult$Sample==sam_list_gain[i],]
  compress2_mult2<- compress2_mult2[order(compress2_mult2$order),]
  x<- compress2_mult2$age_gain[2]- compress2_mult2$age_gain[1]
  all_second<- rbind( all_second, c(sam_list_gain[i], x))
}

head(all_second)
all_second2<- as.data.frame.matrix(all_second)
all_second2$V2<- as.numeric(as.character(all_second2$V2))
median(all_second2$V2[complete.cases(all_second2$V2)])
range(all_second2$V2[complete.cases(all_second2$V2)])
