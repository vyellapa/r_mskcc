

masses = c(A =  331, C =  307, G =  347, T =  322)
probs  = c(A = 0.12, C = 0.38, G = 0.36, T = 0.14)
N  = 100000
sd = 3
nuclt   = sample(length(probs), N, replace = TRUE, prob = probs)
quadwts = rnorm(length(nuclt),
                mean = masses[nuclt],
                sd   = sd)
ggplot(tibble(quadwts = quadwts), aes(x = quadwts)) +
  geom_histogram(bins = 500, fill = "purple")


ggplot(tibble(x = rgamma(10000, shape = 2, rate = 1/3)),
       aes(x = x)) + geom_histogram(bins = 100, fill= "purple")
ggplot(tibble(x = rgamma(10000, shape = 10, rate = 3/2)),
       aes(x = x)) + geom_histogram(bins = 100, fill= "purple")


library(dplyr)
options(scipen=999)
library(GenomicRanges)
setwd("/Users/yellapav/SV_project/data/")
###############################################
###### Create bins ################
cyto.bed=read.table("cytoband_b37.bed",sep="\t")
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

bin.size=1000000


#Create dataframe with bins of bin.size into df called bins.df
bins.df <- data.frame(chrom=character(), 
                      start=numeric(), 
                      stop=numeric(),
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
  bins.df=rbind(bins.df, data.frame(chrom=chroms, start=starts, stop=stops,bin.num=arms))
}

bins.df$bin.num=paste(bins.df$bin.num,bins.df$start,bins.df$stop,sep="_") 



samp="I-H-106917"
nos=8
samp="I-H-130718"
#=nos=16

#samp="I-H-130719"
#nos=8

#samp="I-H-130720"
#nos=10



setwd("/Users/yellapav/Desktop/p220_2019/dirichlet/c2/data/bbg")
t=read.table("I-H-106917-T2-1-D1-2_subclones.txt",sep="\t",header=T) %>% dplyr::mutate(SAMPLE="SAMPLE") %>% dplyr::filter(X=="A")
all = t
for(i in grep(samp,list.files(pattern = "subclones.txt$"), value=T))
{
  print(i)
  r=read.table(i,sep="\t",header=T) 
  r$SAMPLE=i
  all= rbind(all,r)
}

all$SAMPLE=gsub("_subclones.txt","",all$SAMPLE)
all$key = paste(all$SAMPLE,all$chr, all$startpos,all$endpos,sep = ":")
all$length = (all$endpos - all$startpos)/1000000
all = all %>% dplyr::filter(length>=1)

#Overlap with bins.df
asub = all[,c("SAMPLE","chr","startpos","endpos","BAF","LogR","ntot","length","key")] %>% unique()
#asub$key = paste(asub$SAMPLE,asub$chr, asub$startpos,asub$endpos,sep = ":")
#asub$length = (asub$endpos - asub$startpos)/1000000

asub.st = asub %>% dplyr::mutate(st=startpos, sp=startpos+1)
asub.sp = asub %>% dplyr::mutate(st=endpos, sp=endpos+1)
asub.all = rbind(asub.st, asub.sp) %>% unique()




gr.bins = GRanges(seqnames=Rle(bins.df$chrom), IRanges(bins.df$start, bins.df$stop))
gr.sv1 = GRanges(seqnames=Rle(as.character(asub.all$chr)), IRanges(asub.all$st, asub.all$sp))

overlapGenes <- findOverlaps(gr.bins, gr.sv1)
df.sv = data.frame(asub.all[subjectHits(overlapGenes),], bins.df[queryHits(overlapGenes),])

sub.all = df.sv %>% dplyr::mutate(merged.bin=".") %>% dplyr::filter(merged.bin!=".")
cat("Problem if you see rows below..")
sub.all
for(i in unique(df.sv$key)) {
  sub = df.sv %>% dplyr::filter(df.sv$key == i)
  sub$merged.bin = paste(sub$bin.num,collapse=":")
  sub.all=rbind(sub.all,sub)
  
}

sub.all = sub.all %>% dplyr::filter(length>=1)
sub.all %>% group_by(merged.bin,SAMPLE) %>% summarise(n=n()) %>% dplyr::filter(n!=2)

sub.all = sub.all %>% group_by(merged.bin) %>% summarise(n=n()) %>% 
  dplyr::right_join(sub.all,by="merged.bin") %>% 
  dplyr::mutate(status=ifelse(n==nos,"SHARED","UNIQUE"))

sub.all.shared = sub.all %>% dplyr::filter(status=="SHARED")
sub.all.unique = sub.all %>% dplyr::filter(status=="UNIQUE")
shared = all %>% dplyr::filter(key %in% unique(sub.all.shared$key)) %>% dplyr::select(-c("key","length","SAMPLE"))
unique = all %>% dplyr::filter(key %in% unique(sub.all.unique$key)) %>% dplyr::select(-c("key","length","SAMPLE"))

write.table(shared,sprintf("%s_shared.txt",samp),row.names = F,col.names = TRUE,append = F, quote = F,sep = "\t" )
write.table(unique,sprintf("%s_uniqued.txt",samp),row.names = F,col.names = TRUE,append = F, quote = F,sep = "\t" )

