library(ASCAT)
setwd("/Users/yellapav/Desktop/p292/uploads/ascat")
ascat.bc = ascat.loadData("f1.LR","f1.BAF")
ascat.plotRawData(ascat.bc)
ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, "AffyCytoScanHD") 
ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg) 
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc) 


