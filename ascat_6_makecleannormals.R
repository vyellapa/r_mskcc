# Retrieve the list of files for the normal samples included
samples = read.table(file='../caveman/samples.txt', header=F, sep='\t')
normals = read.table(file='normals.txt', header=F, sep='\t')[,2]
normfiles = paste("allelecounts/", normals, ".allelecount.txt",sep="")
ac = read.table(normfiles[1],sep="\t")
# Retrieve the list of SNP positions
SNPpos = ac[,c(1,2)]
counts = matrix(nrow = dim(ac)[1],ncol = length(normfiles))
# For each normal samples, read the coverage across all the SNP positions
for (i in 1:length(normfiles)) {
  print(i)
  ac = read.table(normfiles[i],sep="\t")
  counts[,i] = ac[,7]
}
###########################################################
# some samples may have very low coverage, remove them:
# normals[colSums(counts)<median(colSums(counts))/2]
# PD21396b PD21397b
counts = counts[,!grepl('PD21396b', normfiles) & !grepl('PD21397b',normfiles)]
normals = normals[!(normals %in% c('PD21396b','PD21397b'))]
normfiles = paste("allelecounts/", normals, ".allelecount.txt",sep="")

# also may need to remove additional samples because of noise or aberrations, or because they were duplicates
# this is now partially automated!
# for the TMD project, removed 8 samples in the end (because of noise, artefacts or aberrations):
# PD13635, PD13638, PD13642, PD13656, PD136707, PD13711, PD13714 and PD13717
###########################################################
sum(rowSums(counts>=10)>=length(normfiles)*0.9)
# 65429 positions now (with 92 normals) since two were dropped due to low coverage
# Considering only those SNPs with at least 10x coverage,
# select those with total coverage >= 90% of the number of normal samples
hc = which(rowSums(counts>=10)>=length(normfiles)*0.9)
SNPpos = SNPpos[hc,]
countn = counts[hc,]

bafn = matrix(nrow = length(hc),ncol = length(normfiles))
for (i in 1:length(normfiles)) {
  print(i)
  ac = read.table(normfiles[i],sep="\t")
  ac = ac[hc,]
  acgt = ac[,c(3:6)]
  acgts = t(apply(acgt,1,sort))
  bafn[,i] = acgts[,4]/(acgts[,3]+acgts[,4])
}
colnames(bafn) = normals

bafn[is.na(bafn)]=1

# not more than 2 samples misbehaving and at least two heterozygous
#old:
#sum(rowSums(bafn<0.9&bafn>0.6|bafn>0.1&bafn<0.4)<3&rowSums(bafn<0.6&bafn>0.4)>=2)
#new for genomes:
sum(rowSums(countn>=20&(bafn<0.9&bafn>0.68|bafn>0.1&bafn<0.32))<2.5*0.053876*rowSums(countn>=20&bafn<0.9&bafn>0.1))
# For TMD 3168 probes now
# Lucy's data set, 13681 probes for now

hc = which(rowSums(countn>=20&(bafn<0.9&bafn>0.68|bafn>0.1&bafn<0.32))<2.5*0.053876*rowSums(countn>=20&bafn<0.9&bafn>0.1))
SNPpos = SNPpos[hc,]

countn = countn[hc,]

bafn = bafn[hc,]

countn[countn==0] = 1



# throw out SNPs too close (depending on freq in population?? If it's the same SNP, no matter what freq? Maybe if freq is same, remove?
# probably a (near?-)perfect overlap between cases inside [0.1,0.9] would be a good criterion
# needs a way to flag the cases first
d = diff(SNPpos[,2])

#> sum(d<100)
#[1] 861
# there are many of these pairs!

poskes = which(d<=75)
todel = NULL
for (poske in poskes) {
  #print(poske)
  bafnp1 = bafn[poske,]
  bafnp2 = bafn[poske+1,]
  hetp1 = ifelse(bafnp1>0.1&bafnp1<0.9,1,0)
  hetp2 = ifelse(bafnp2>0.1&bafnp2<0.9,1,0)
  # if less than 10% of the calls are different
  if(sum(abs(hetp1-hetp2))<0.05*(sum(hetp1)+sum(hetp2))) {
    prevp = max(poske-1,1)
    nextp = min(poske+1,length(d))
    if(d[prevp]<d[nextp]) {
      todel = c(todel,poske)
    }
    else {
      todel = c(todel,poske+1)
    }
  }
}
todel = unique(todel)

length(todel)

SNPpos = SNPpos[-todel,]
countn = countn[-todel,]
bafn = bafn[-todel,]
# this deletes 126 probes


countn[countn==0] = 1


LogRn = matrix(nrow = dim(countn)[1],ncol=dim(countn)[2])
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = log(countn[,i],2)
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}

# create normal model:
LogRref = apply(LogRn,1,median)
# median should be more stable than mean

# renormalize normals
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = LogRn[,i]-LogRref
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}

colnames(SNPpos) = c("chr","pos")
rownames(SNPpos) = paste("snp",1:dim(SNPpos)[1],sep="")

colnames(LogRn) = normals


# look at standard deviation and remove outliers
sdkes = apply(LogRn,1,sd)

sum(sdkes>2*median(sdkes))
# some debate about whether this should be 1.5x or 2x..
# ~ 424 probes, quite a lot..


# remove and recalculate:
SNPpos = SNPpos[sdkes<2*median(sdkes),]
countn = countn[sdkes<2*median(sdkes),]
bafn = bafn[sdkes<2*median(sdkes),]

LogRn = matrix(nrow = dim(countn)[1],ncol=dim(counts)[2])
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = log(countn[,i],2)
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}

# create normal model:
LogRref = apply(LogRn,1,median)
# median should be more stable than mean

# renormalize normals
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = LogRn[,i]-LogRref
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}

sdkes = apply(LogRn,1,sd)


colnames(LogRn) = normals
rownames(LogRn) = rownames(SNPpos)



# filter samples on high SD?
sdsam = apply(LogRn,2,sd)
# covers both noise and waves

names(which(sdsam>1.5*median(sdsam)))
# note the the output below is with the old samples in!
# [1] "PD13635b" "PD13638b"
# so, these are 2 samples to remove.
# For Lucy's data set, thiis removes 14 normals.
bafn = cbind(SNPpos,bafn)

write.table(bafn,"Normal_BAF.txt",sep="\t",col.names=NA,row.names=T,quote=F)

LogRn = cbind(SNPpos,LogRn)

write.table(LogRn,"Normal_LogR.txt",sep="\t",col.names=NA,row.names=T,quote=F)
