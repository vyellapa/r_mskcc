library(dplyr)
library(GenomicRanges)
library("biomaRt"); library("RColorBrewer"); library(reshape2) ;library(data.table); library(limma)


setwd("/Users/yellapav/Desktop/p220_2019/misc/files")

driver_p = fread("driver_genes/driver_fra_morgan.txt",sep="\t",header=T)
driver_c = fread("driver_genes/reference_cnv_prediction.txt",sep="\t",header=T)

NS = c("missense_variant","splice_acceptor_variant","splice_donor_variant","start_lost","stop_gained","stop_lost","stop_retained_variant")

s = list()
## SNVs #####
for(i in list.files(path="snvs/filtered/",pattern = ".txt$",full.names = T)) {
  print(i)
  temp = fread(i,sep="\t",header=F)
  s[[i]] = temp
}

snv = do.call(rbind,s) %>% dplyr::select(V2,V3,V4,V5,V6,V7,V8,V15,V20,V21,V22,V24,V26,V27,V28,V29,V32,V33) %>% distinct()
snv = snv[grep(paste(NS,collapse = "|"), snv$V24),]
snv = snv %>% dplyr::select(V2,V3,V4,V5,V6,V7,V8,V20,V24,V26,V32,V33) %>% distinct()

snv = snv[snv$V26 %in% unique(driver_p$Gene.Symbol),]
colnames(snv) = c("TARGET_NAME","REFERENCE_NAME","CHR","START","END","REF","ALT","VAF","VEP_Consequence","VEP_SYMBOL","VEP_HGVSc","VEP_HGVSp") 
snv.df.wgs = snv %>% dplyr::mutate(assay="WGS") 



#### Indels #####
s = list()
for(i in list.files(path="indels/",pattern = ".txt$",full.names = T)) {
  print(i)
  temp = fread(i,sep="\t",header=T)
  s[[i]] = temp
}

indels = do.call(rbind,s)
indels = indels[indels$VEP_SYMBOL %in% unique(driver_p$Gene.Symbol),]
NSI = c("frameshift_variant","splice_acceptor_variant","stop_retained_variant","splice_donor_variant","stop_gained","transcript_ablation","stop_lost","start_lost","start_retained_variant")
indels.df.wgs = indels[grep(paste(NSI,collapse = "|"), indels$VEP_Consequence),c(2,3,4,5,6,7,8,16,29,31,37,38)] %>% distinct() %>% dplyr::mutate(assay="WGS") 




######## CNV #########
s = list()
for(i in list.files(path="battenberg/",pattern = "_subclones.txt$",full.names = T)) {
  print(i)
  temp = fread(i,sep="\t",header=T)
  temp$sample = as.character(strsplit2(basename(i),"_")[1])
  s[[i]] = temp
}

cnv = do.call(rbind,s)
cnv$V1 = cnv$sample
cnv = cnv[,1:14] %>% dplyr::select(-pval)
cnv.wgs.bkup = cnv

gr.seg = GRanges(seqnames=Rle(driver_c$chrom), IRanges(as.numeric(as.character(driver_c$start)), as.numeric(as.character(driver_c$end))))
gr.cyto = GRanges(seqnames=Rle(cnv$chr), IRanges(as.numeric(as.character(cnv$startpos)), as.numeric(as.character(cnv$endpos))))
overlapGenes <- findOverlaps(gr.cyto, gr.seg)
cnv.df = data.frame(driver_c[subjectHits(overlapGenes),], cnv[queryHits(overlapGenes),]) %>% distinct() 

cnv.cols=c("V1","chr","startpos","endpos","BAF","LogR","ntot","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","chrom","start","end","gene","band","effect_original","code")  
cnv.df = cnv.df[,cnv.cols]

bl1 = cnv.df %>% dplyr::filter(nMin1_A==0 | nMin2_A==0) %>% dplyr::filter(effect_original=="loss")
bg1 = cnv.df %>% dplyr::filter(nMin1_A>1 | nMin2_A>1 | nMaj1_A > 2 | nMaj2_A>2 | (nMaj1_A > 1 & nMin1_A>0) | (nMaj2_A > 1 & nMin2_A>0)) %>% dplyr::filter(effect_original=="gain")

cnv.df.wgs = unique(rbind(bl1,bg1)) %>% dplyr::rename(Sample=V1) %>% dplyr::mutate(assay="WGS") 





######## SV #########
s = list()
for(i in list.files(path="brass/",pattern = "*.bedpe$",full.names = T)) {
  print(i)
  temp = read.table(i,sep="\t",header=F)
  #temp$sample = as.character(strsplit2(basename(i),"_")[1])
  temp = temp[,c(1,2,3,4,5,6,8,9,10,11,12,28,37)]
  s[[i]] = temp
}

sv = do.call(rbind,s)


sv1 = read.table("multi_sample/I-H-130720_consolidated.txt",sep="\t",header=F)
sv2 = read.table("multi_sample/I-H-130718_consolidated.txt",sep="\t",header=F)

sys.call("less /Users/yellapav/Desktop/p220_2019/misc/files/multi_sample/I-H-130720_consolidated.txt |awk -F'\t' '{if(($2==4 && $5==14) || ($1~/^id/)) print $0}'|cut -f2-10,12-17,24,32")
sys.call("less /Users/yellapav/Desktop/p220_2019/misc/files/multi_sample/I-H-130718_consolidated.txt |awk '{if(($2==4 && $5==14) || ($1~/id/)) print $0}'|cut -f2-10,12-19,28,37")

#less /Users/yellapav/Desktop/p220_2019/misc/files/multi_sample/I-H-130720_consolidated.txt |awk -F'\t' '{if(($2==4 && $5==14) || ($1~/^id/)) print $0}'|cut -f2-10,12-17,24,32 > sv_I130720.txt
#less /Users/yellapav/Desktop/p220_2019/misc/files/multi_sample/I-H-130718_consolidated.txt |awk '{if(($2==4 && $5==14) || ($1~/id/)) print $0}'|cut -f2-10,12-19,28,37 > sv_I130718.txt



######### Weinhold ############
######## CNV #########
s = list()
setwd("/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/weinhold")
for(i in list.files(path="bbg/",pattern = "_prof.txt$",full.names = T)) {
  print(i)
  temp = fread(i,sep="\t",header=T)
  temp$sample = as.character(strsplit2(basename(i),"_p")[1])
  s[[i]] = temp
}

wcnv = do.call(rbind,s)
wcnv$V1 = wcnv$sample
wcnv = wcnv[,1:14] %>% dplyr::select(-pval)
cnv.exome.bkup = wcnv

gr.seg = GRanges(seqnames=Rle(driver_c$chrom), IRanges(as.numeric(as.character(driver_c$start)), as.numeric(as.character(driver_c$end))))
gr.cyto = GRanges(seqnames=Rle(wcnv$chr), IRanges(as.numeric(as.character(wcnv$startpos)), as.numeric(as.character(wcnv$endpos))))
overlapGenes <- findOverlaps(gr.cyto, gr.seg)
wcnv.df = data.frame(driver_c[subjectHits(overlapGenes),], wcnv[queryHits(overlapGenes),]) %>% distinct() 

cnv.cols=c("V1","chr","startpos","endpos","BAF","LogR","ntot","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","chrom","start","end","gene","band","effect_original","code")  
wcnv.df = wcnv.df[,cnv.cols]

wcnv.df = wcnv.df %>% dplyr::filter(!(is.na(nMaj2_A))) %>% dplyr::rename(Sample=V1)
head(wcnv.df)
#l1 = wcnv.df%>% dplyr::filter(LogR<0 & effect_original=="loss")
l2 = wcnv.df%>% dplyr::filter(nMin2_A==0 & effect_original=="loss")
#loss=rbind(l1,l2) %>% distinct()
loss=l2 %>% distinct()

#g1 = wcnv.df  %>% dplyr::filter(LogR>0 & effect_original=="gain")
g2 = wcnv.df  %>% dplyr::filter(nMaj2_A>2 & effect_original=="gain")

#gain = rbind(g1,g2) %>% distinct()
gain = g2 %>% distinct()

wcnv.df.results = unique(rbind(loss,gain)) %>% dplyr::mutate(assay="WXS")

cnv.df.all = distinct(rbind(wcnv.df.results, cnv.df.wgs))


#### Indels #####
s = list()
for(i in list.files(path="snv_indel/",pattern = "indel.txt$",full.names = T)) {
  print(i)
  temp = fread(i,sep="\t",header=T)
  s[[i]] = temp
}

indels = do.call(rbind,s)
indels = indels[indels$VEP_SYMBOL %in% unique(driver_p$Gene.Symbol),]
NSI = c("frameshift_variant","splice_acceptor_variant","stop_retained_variant","splice_donor_variant","stop_gained","transcript_ablation","stop_lost","start_lost","start_retained_variant")

indels.df.wxs = indels[grep(paste(NSI,collapse = "|"), indels$VEP_Consequence),] %>% distinct() %>% dplyr::mutate(assay="WXS") %>% dplyr::rename(TARGET_VAF_MEAN=pindel_TARGET_VAF)
indels.df.all = distinct(rbind(indels.df.wgs,indels.df.wxs))



## SNVs #####
for(i in list.files(path="snv_indel/",pattern = "snv.txt$",full.names = T)) {
  print(i)
  snv = fread(i,sep="\t",header=T)
}
NS = c("missense_variant","splice_acceptor_variant","splice_donor_variant","start_lost","stop_gained","stop_lost","stop_retained_variant")

snv = snv[grep(paste(NS,collapse = "|"), snv$VEP_Consequence),] %>% distinct()

snv = snv[snv$VEP_SYMBOL %in% unique(driver_p$Gene.Symbol),] %>% dplyr::mutate(assay="WXS") 
colnames(snv) = c("TARGET_NAME","REFERENCE_NAME","CHR","START","END","REF","ALT","VAF","VEP_Consequence","VEP_SYMBOL","VEP_HGVSc","VEP_HGVSp","assay") 
snv.df.wxs = snv

snv.df.all = distinct(rbind(snv.df.wxs,snv.df.wgs))

write.table(snv.df.all,file="/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/autopsy_supplementary_snv.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
write.table(cnv.df.all,file="/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/autopsy_supplementary_cnv.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
write.table(indels.df.all,file="/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/autopsy_supplementary_indels.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)



###### IGH ########
isa.ex = read.table("/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/weinhold/samples_isabl.txt", header = T, sep = "\t") %>% 
  dplyr::filter(Sample.Category=="TUMOR") %>%
  dplyr::select(Experiment.System.ID,Individual.ID) %>% distinct() %>% 
  dplyr::rename(Sample = Experiment.System.ID,iid=Individual.ID)

igh.ex = read.table("/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/weinhold/igh.txt", header = F, sep = "\t") %>% dplyr::rename(iid = V1,igh = V2)
igh.ex = dplyr::left_join(isa.ex,igh.ex, by = "iid") %>% dplyr::select(-iid) %>% distinct() %>% dplyr::mutate(type="WXS") 

isa.wg = read.table("/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/weinhold/samples_wgs.txt", header = T, sep = "\t") %>%
  dplyr::filter(Sample.Category=="TUMOR" & Sequencing.Technique=="DNA|WG") %>%
  dplyr::select(Experiment.System.ID,Individual.ID) %>% distinct() %>% 
  dplyr::rename(Sample = Experiment.System.ID,iid=Individual.ID) %>% dplyr::select(Sample) %>% dplyr::mutate(igh=0,type="WGS") 
isa.wg[grep("130718|130720",isa.wg$Sample),]$igh=1

igh.wg = isa.wg
sv.df.all = rbind(igh.wg, igh.ex)


############ Hyperdiploid ###############
cyto = read.table("/Users/yellapav//local/resources/cytoband_b37.bed",header = F,sep ="\t") %>% distinct() %>% dplyr::mutate(diff = V3-V2)
sizes = cyto %>% group_by(V1) %>% summarise(size = sum(diff))

hypc = c("3" ,"5","7" ,"9" ,"11" ,"15" , "19" ,"21")
cnv.exome.hyp = cnv.exome.bkup %>% dplyr::filter(nMaj2_A>2)
cnv.wgs.hyp = cnv.wgs.bkup %>% dplyr::filter(nMin1_A>1 | nMin2_A>1 | nMaj1_A > 2 | nMaj2_A>2 | (nMaj1_A > 1 & nMin1_A>0) | (nMaj2_A > 1 & nMin2_A>0))
cnv.hyp = distinct(rbind(cnv.exome.hyp,cnv.wgs.hyp))
cnv.hyp = distinct(cnv.wgs.hyp)
cnv.hyp = cnv.hyp[as.character(cnv.hyp$chr) %in% hypc, ]
hyp.wgs = cnv.hyp %>% dplyr::mutate(diff = endpos-startpos) %>% group_by(V1,chr) %>% 
  summarise(size = sum(diff)) %>% dplyr::left_join( sizes, by=c("chr"="V1")) %>% 
  dplyr::mutate(fraction = size.x/size.y) %>% 
  dplyr::mutate(amp=ifelse(fraction>0.6,1,0)) %>% 
  group_by(V1) %>% summarise(no.amp=sum(amp)) %>% 
  dplyr::filter(no.amp>1) %>% 
  dplyr::mutate(hyper=1) %>% 
  dplyr::select(-no.amp) %>% 
  dplyr::rename(Sample=V1) %>% as.data.frame()

hyp.ex = read.table("/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/weinhold/hyp.txt", header = F, sep = "\t") %>% dplyr::rename(iid = V1,hyper = V2)
hyp.ex = dplyr::left_join(isa.ex,hyp.ex, by = "iid") %>% dplyr::select(-iid) %>% distinct() 
hyp.all = distinct(rbind(hyp.wgs,hyp.ex))
  
##########################################
t = indels.df.all
t$status=1
t = t %>% dplyr::select(TARGET_NAME,VEP_SYMBOL, status) %>% dplyr::rename(Sample=TARGET_NAME)
t.indel = t
indel.t = dcast(t,Sample~VEP_SYMBOL, value.var="status")
##colnames(indel.t)=c("Sample",paste0(grep("Sample",colnames(indel.t),value=T,invert=T),".indel"))


t = snv.df.all
t$status=1
t = t %>% dplyr::select(TARGET_NAME,VEP_SYMBOL, status) %>% dplyr::rename(Sample=TARGET_NAME)
t.snv = t
snv.t = dcast(t,Sample~VEP_SYMBOL, value.var="status")
##colnames(snv.t)=c("Sample",paste0(grep("Sample",colnames(snv.t),value=T,invert=T),".snv"))
t =  distinct(rbind(t.snv,t.indel))
mut.t = dcast(t,Sample~VEP_SYMBOL, value.var="status")


t = cnv.df.all
t$status=1
t = t %>% dplyr::select(Sample,code, status) %>% distinct()
cnv.t = dcast(t,Sample~code, value.var="status")

#colnames(snv.t)=c("Sample",paste0(grep("Sample",colnames(snv.t),value=T,invert=T),".snv"))

supp.final = left_join(sv.df.all,mut.t)
supp.final = left_join(supp.final,cnv.t)
supp.final = left_join(supp.final,hyp.all)
supp.final[is.na(supp.final)]=0

write.table(supp.final,file="/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/autopsy_supplementary_all.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)
write.table(supp.final,file="/Users/yellapav/Desktop/p220_2019/supplementary_ncomm/autopsy_supplementary_all_hyp2.txt", append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE,quote = FALSE)







r = read.table("~/Desktop/p220_2019/weinhold_facets_dp/dp/IID_H159207_T01_01_WE01_DPoutput_10000iters_1000burnin/IID_H159207_T01_01_WE01_10000iters_1000burnin_bestConsensusAssignments.bed",sep="\t",header = T)
r = read.table("~/Desktop/p220_2019/weinhold_facets_dp/dp/IID_H159207_T01_01_WE01_DPoutput_10000iters_1000burnin/IID_H159207_T01_01_WE01_10000iters_1000burnin_bestConsensusAssignments.bed",sep="\t",header = T)
