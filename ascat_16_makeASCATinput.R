# Retrieve the list of files for the normal samples and exclude the ones in the black list
normals = read.table(file='normals.txt', header=F, sep='\t')[,2]
normfiles = paste("allelecounts_final/", normals, ".allelecount.txt",sep="")

ac = read.table(normfiles[1],sep="\t")
SNPpos = ac[,c(1,2)]
counts = matrix(nrow = dim(ac)[1],ncol = length(normfiles))
for (i in 1:length(normfiles)) {
  print(i)
  ac = read.table(normfiles[i],sep="\t")
  counts[,i] = ac[,7]
}
countn = counts
countn[countn==0]=1

LogRn = matrix(nrow = dim(countn)[1],ncol=dim(countn)[2])
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = log(countn[,i],2)
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}

# create normal model:
LogRref = apply(LogRn,1,median)
# renormalize normals
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = LogRn[,i]-LogRref
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}
colnames(SNPpos) = c("chr","pos")
rownames(SNPpos) = paste("snp",1:dim(SNPpos)[1],sep="")
colnames(LogRn) = normals

# Retrieve the list of files for the tumour samples included
tumors = read.table(file='../caveman/samples.txt', header=F, sep='\t')[,2]
tumorfiles = paste('allelecounts_final/', tumors,".allelecount.txt",sep="")

counts = matrix(nrow = dim(ac)[1],ncol = length(tumorfiles))
for (i in 1:length(tumorfiles)) {
  print(paste(i,": ",tumors[i],sep=""))
  ac = read.table(tumorfiles[i],sep="\t")
  counts[,i] = ac[,7]
}
countn = counts
countn[countn==0]=1

LogRn = matrix(nrow = dim(countn)[1],ncol=dim(countn)[2])
for (i in 1:dim(countn)[2]) {
  LogRn[,i] = log(countn[,i],2)-LogRref
  LogRn[,i] = LogRn[,i]-mean(LogRn[,i])
}
colnames(LogRn) = tumors

bafn = matrix(nrow = dim(countn)[1],ncol = length(tumors))
for (i in 1:length(tumorfiles)) {
  print(paste(i,": ",tumors[i],sep=""))
  ac = read.table(tumorfiles[i],sep="\t")
  acgt = ac[,c(3:6)]
  acgts = t(apply(acgt,1,sort))
  bafn[,i] = acgts[,4]/(acgts[,3]+acgts[,4])
}
colnames(bafn) = tumors
bafn = cbind(SNPpos,bafn)
write.table(bafn,"Tumor_BAF.txt",sep="\t",col.names=NA,row.names=T,quote=F)
LogRn = cbind(SNPpos,LogRn)
write.table(LogRn,"Tumor_LogR.txt",sep="\t",col.names=NA,row.names=T,quote=F)
