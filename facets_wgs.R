library(facets)
suppressMessages(library(optparse, quietly = TRUE))
suppressMessages(library(getopt, quietly = TRUE))


option_list = list(make_option(c("-i", "--input"), type = "character",
                               help = "input file",
                               metavar = "path", action = "store"),
                   make_option(c("-o", "--output"), type = "character",
                               help = "output file", metavar = "path"),
		   make_option(c("-r", "--rho_psi_out"), type = "character",
                               help = "Intermediate parameters output", metavar = "path"),
                   make_option(c("-p", "--png"), type = "character",
                               help = "output png file", metavar = "path"),
                   make_option(c("-d","--diplogR"), type = "double",
                               help = "alternate diplogR to use",
                               default = 999.9, metavar = "double"));



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);




## Change the running directory accordingly
cur_dir=getwd()
setwd(cur_dir)
cat(getwd())
cat("\n")
#cat(opt$s)
#datafile
rcmat = readSnpMatrix(opt$i)
xx=preProcSample(rcmat, ndepth=20,gbuild="hg19",ndepthmax=950, snp.nbhd=300)

if(opt$d==999.9) {
  oo=procSample(xx,cval=150)
  def = 0.0
  altLogR = oo$alBalLogR[[1]]

} else {
  oo=procSample(xx,cval=150, dipLogR = opt$d)
  def = opt$d
  altLogR = oo$dipLogR[[1]]
}

fit=emcncf(oo)
cn=fit$cncf


#Print purity and ploidy into a file named purity_and_ploidy.txt
#purity = gsub(".txt.gz",".purity.txt",basename(opt$i))

#pur_out=paste(cur_dir,"/",opt$s,sep="")
#cat(opt$r)

fileConn<-file(opt$r)
  writeLines(c("Purity",fit$purity,"\nPloidy",fit$ploidy, "\nDiploid LogR",fit$dipLogR,"\nAlternate Diploid LogR", oo$alBalLogR[[1]]), fileConn)
close(fileConn)


#Adjust based on overall log ratio
cn$log2Fold_default=cn$cnlr.median.clust-oo$dipLogR[[1]]
cn$log2Fold_alternate=cn$cnlr.median.clust-altLogR
cn$log2Fold_user=cn$cnlr.median.clust-def

if(opt$d==999.9) { 
	def = oo$dipLogR[[1]]
} 

cn$log2Fold_toUse=cn$cnlr.median.clust-def

write.table(cn,file=opt$o, append=FALSE, sep="\t", eol="\n", row.names=F, col.names=T)

png(file = opt$p, res=300, width=2000, height=1500)
  plotSample(x=oo,emfit=fit)
dev.off()
