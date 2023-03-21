library(facets)

args=commandArgs(TRUE)
fName=args[1]
oFile=args[2]
pngName=args[3]

setwd("/Users/yellapav/Desktop/selected")
#datafile = system.file("extdata","st.csv.gz" , package="facets")
#datafile = system.file("extdata", "stomach.csv.gz", package="facets")

#datafile
rcmat = readSnpMatrix(fName)
xx=preProcSample(rcmat, ndepth=50,gbuild="hg19",ndepthmax=1500, unmatched = TRUE)
oo=procSample(xx,cval=150)
fit=emcncf(oo)

write.table(fit$cncf,file=oFile, append=FALSE, sep="\t", eol="\n", row.names=F, col.names=T)



png(file = pngName, res=300, width=1800, height=1400)
plotSample(x=oo,emfit=fit)

dev.off()


