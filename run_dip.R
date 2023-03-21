# Put this in the directory where the files from step1 are and type Rscript run_dip.R and hit enter

library(facets)
#setwd("/Users/yellapav/Desktop/selected/selected_v2_105")

for(i in list.files(pattern = "*.gz")) { 

oFile=gsub(".gz", "", i);
pngName=gsub(".txt",".png",oFile);
print(c(oFile,pngName,i))

#datafile
rcmat = readSnpMatrix(i)
xx=preProcSample(rcmat, ndepth=50,gbuild="hg19",ndepthmax=1500, unmatched = TRUE)
oo=procSample(xx,cval=150)
fit=emcncf(oo)
cn=fit$cncf
cn$cnlr.median.clust=cn$cnlr.median.clust-oo$dipLogR[[1]]


cat("Purity is: ")
cat(fit$purity)
cat("\n")
cat("Ploidy is: ")
cat(fit$ploidy)
cat("\n")

write.table(cn,file=oFile, append=FALSE, sep="\t", eol="\n", row.names=F, col.names=T)



png(file = pngName, res=300, width=2000, height=1500)
plotSample(x=oo,emfit=fit)

dev.off()

}
