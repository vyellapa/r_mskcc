# Retrieve the list of files for the normal samples included
samples = read.table(file='../caveman/samples.txt', header=F, sep='\t')
normals = read.table(file='normals.txt', header=F, sep='\t')[,2]
normfiles = paste("allelecounts2/", normals, ".allelecount.txt",sep="")
ac = read.table(normfiles[1],sep="\t")
SNPpos = ac[,c(1,2)]
counts = matrix(nrow = dim(ac)[1],ncol = length(normfiles))
for (i in 1:length(normfiles)) {
  print(i)
  ac = read.table(normfiles[i],sep="\t")
  counts[,i] = ac[,7]
}

bafn = matrix(nrow = dim(counts)[1],ncol = length(normfiles))
for (i in 1:length(normfiles)) {
  print(i)
  ac = read.table(normfiles[i],sep="\t")
  acgt = ac[,c(3:6)]
  acgts = t(apply(acgt,1,sort))
  bafn[,i] = acgts[,4]/(acgts[,3]+acgts[,4])
}
colnames(bafn) = normals
bafn[is.na(bafn)]=1

# now, use SNP probes (selected from cases with heterozygous SNPs) to normalize LogR:
# + check for variance of LogR of SNPs and choose comparable CN only probes
hc = rowSums(counts>=20&(bafn<0.9&bafn>0.68|bafn>0.1&bafn<0.32))<2.5*0.053876*rowSums(counts>=20&bafn<0.9&bafn>0.1)
countn = counts
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
  LogRn[,i] = LogRn[,i]-mean(LogRn[hc,i])
}
colnames(LogRn) = normals

sdevSNP = apply(LogRn[hc,],1,sd)
sdevall = apply(LogRn,1,sd)

sum(sdevall>2*median(sdevSNP))
sum(sdevSNP>2*median(sdevSNP))
# these are about in proportion (well, about..).. -> yet keep condition for now..
# this run removes 1182 from all, 157 from SNP - but that's ok, median(sdev) are about equal and SNPs have already been filtered).

hc = sdevall<2*median(sdevSNP)
bafn = bafn[hc,]
LogRn = LogRn[hc,]
SNPpos = SNPpos[hc,]
colnames(SNPpos) = c("chr","pos")
rownames(SNPpos) = paste("snp",1:dim(SNPpos)[1],sep="")
colnames(LogRn) = normals

#renormalize LogR:
LogRref = apply(LogRn,1,median)
# renormalize normals
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = LogRn[,i]-LogRref
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}
bafn = cbind(SNPpos,bafn)
write.table(bafn,"Normal_BAF.txt",sep="\t",col.names=NA,row.names=T,quote=F)
LogRn = cbind(SNPpos,LogRn)
write.table(LogRn,"Normal_LogR.txt",sep="\t",col.names=NA,row.names=T,quote=F)
