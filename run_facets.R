library(facets)

args=commandArgs(TRUE)
fName=args[1]  ## Output file from step 1 (snp-pileup program)
oFile=args[2]  ## Output filename from facets
pngName=args[3] ## PNG output filename


## Change the running directory accordingly
setwd("/Users/yellapav/Desktop/selected/selected_v2_105")

#datafile
rcmat = readSnpMatrix(fName)
xx=preProcSample(rcmat, ndepth=50,gbuild="hg19",ndepthmax=1500, unmatched = TRUE)
oo=procSample(xx,cval=150)
fit=emcncf(oo)
cn=fit$cncf

#Adjust based on overall log ratio
cn$cnlr.median.clust=cn$cnlr.median.clust-oo$dipLogR[[1]]


write.table(cn,file=oFile, append=FALSE, sep="\t", eol="\n", row.names=F, col.names=T)



png(file = pngName, res=300, width=2000, height=1500)
plotSample(x=oo,emfit=fit)

dev.off()

