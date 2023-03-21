getMutTables <- function(myFile, onlyPASSED=FALSE, onlyPassedAsrd=FALSE, genome.v="hg19", genomeSeq=Hsapiens, addContext=TRUE) {

                                        # plots mutation-context for all variants in the vcf file
                                        # and separately for the variants that passed

# load the vcf file    
B35.1.2.vcf <- readVcf(myFile, genome.v)

# Filter variants that failed the basic post-processing flags
rd <- rowRanges(B35.1.2.vcf)
fs <- rd$FILTER
fs.passed <- (fs=='PASS')

if (onlyPASSED) {
        
    B35.1.2.vcf <- B35.1.2.vcf[fs.passed,]
    rd <- rowRanges(B35.1.2.vcf)
    fs <- rd$FILTER
    fs.passed <- (fs=='PASS')
    
}
# Use soft filter for ASRD>0.90
info.data <- info(B35.1.2.vcf)
if(onlyPassedAsrd) {
	asrd.passed = info.data$ASRD>0.90
	B35.1.2.vcf <- B35.1.2.vcf[asrd.passed,]
	fs.passed <- (fs=='PASS' & asrd.passed)
}

info.data <- info(B35.1.2.vcf)
rd <- rowRanges(B35.1.2.vcf)
fs <- rd$FILTER
fs.passed <- (fs=='PASS')
rgs <- ranges(rd)
starts <- start(rgs)
ends <-  end(rgs)
chroms <- paste('chr', seqnames(rd), sep='')

fxd <- (fixed(B35.1.2.vcf))
wt <- as.character(rd$REF)
mt <- as.character(unlist(rd$ALT))

barcode <- paste(chroms, '-',starts,'-', mt, sep='')

if (addContext) {
    bb <- as.character(getSeq(genomeSeq, chroms, start=starts-1, end=ends-1))
    ba <- as.character(getSeq(genomeSeq, chroms, start=starts+1, end=ends+1))
    wt.ref <- as.character(getSeq(genomeSeq, chroms, start=starts, end=ends))
    triplets <- as.character(getSeq(genomeSeq, chroms, start=starts-1, end=ends+1))

                                        # check the annotation
    if (sum(!wt.ref==wt)>0) {
        cat('wrong reference genome \n')
        browser()
    }
                                        # 
    mut.table <- data.frame(bbef=as.character(bb), wt=as.character(wt), mt=as.character(mt), baft=as.character(ba), stringsAsFactors=FALSE)
    
    
    mut.table$pyrwt <- as.character(mut.table$wt)
    mut.table$pyrmut <- as.character(mut.table$mt)
    mut.table$pyrbbef <- as.character(mut.table$bbef)
    mut.table$pyrbaft <- as.character(mut.table$baft)
    
    
# the mutations originally not in pyramidine contex
    not.pyr <- ((wt=='G') | (wt=='A'))
    mut.table$pyrwt[not.pyr] <- as.character(toPyr(mut.table$wt[not.pyr]))
    mut.table$pyrmut[not.pyr] <- as.character(toPyr(mut.table$mt[not.pyr]))
    mut.table$pyrbbef[not.pyr] <- as.character(toPyr(mut.table$baft[not.pyr]))
    mut.table$pyrbaft[not.pyr] <- as.character(toPyr(mut.table$bbef[not.pyr]))
    
    all.hist <- generateHist(mut.table, normalise=FALSE)
    passed.hist <- generateHist(mut.table[fs.passed,], normalise=FALSE)
    if (sum(info.data[,'SNP'])>0) {
        snps.hist <-  generateHist(mut.table[info.data[,'SNP'],], normalise=FALSE)
        names(snps.hist) <- mut.order
    } else {
        snps.hist <- NULL
    }
    non.snps.hist <-  generateHist(mut.table[!info.data[,'SNP'],], normalise=FALSE)
    
    names(passed.hist) <- mut.order
    names(non.snps.hist) <- mut.order

    muts <- data.frame(chroms=chroms, starts=starts, ends = ends, wt=wt, mt=mt, pyrwt=mut.table$pyrwt , pyrmut=mut.table$pyrmut, pass=fs.passed, barcode=barcode,
                   context=paste(mut.table$pyrbbef, '[',mut.table$pyrwt, '>',mut.table$pyrmut , ']', mut.table$pyrbaft,sep=''),
                   tumor.freq=rep(NA,length(chroms)), normal.freq=rep(NA,length(chroms)),
                   tumor.reads=rep(NA,length(chroms)), normal.reads=rep(NA,length(chroms)),
                   tumor.depth=rep(NA,length(chroms)), normal.depth=rep(NA,length(chroms)),
                                                           isSnp =rep(NA,length(chroms))
                   )
} else {
    
    mut.table <- NULL
    
    muts <- data.frame(chroms=chroms, starts=starts, ends = ends, wt=wt, mt=mt,
                       barcode=barcode,
                       tumor.freq=rep(NA,length(chroms)), normal.freq=rep(NA,length(chroms)),
                       tumor.reads=rep(NA,length(chroms)), normal.reads=rep(NA,length(chroms)),
                       tumor.depth=rep(NA,length(chroms)), normal.depth=rep(NA,length(chroms)),
                       isSnp =rep(NA,length(chroms)), filters=fs
                       )

    all.hist <- NULL 
    passed.hist <- NULL 
    mut.table <- NULL
    snps.hist <- NULL
    non.snps.hist <- NULL
    
}

# calculate normal frequency
geno.data <- geno(B35.1.2.vcf)

if ('FAZ' %in% names(geno.data)) {
    muts$tumor.depth <- geno.data[['FAZ']][,'TUMOUR'] + geno.data[['RAZ']][,'TUMOUR'] +
        geno.data[['FCZ']][,'TUMOUR'] + geno.data[['RCZ']][,'TUMOUR'] +
            geno.data[['FGZ']][,'TUMOUR'] + geno.data[['RGZ']][,'TUMOUR'] +
                geno.data[['FTZ']][,'TUMOUR'] + geno.data[['RTZ']][,'TUMOUR']
    muts$normal.depth <- geno.data[['FAZ']][,'NORMAL'] + geno.data[['RAZ']][,'NORMAL'] +
        geno.data[['FCZ']][,'NORMAL'] + geno.data[['RCZ']][,'NORMAL'] +
            geno.data[['FGZ']][,'NORMAL'] + geno.data[['RGZ']][,'NORMAL'] +
                geno.data[['FTZ']][,'NORMAL'] + geno.data[['RTZ']][,'NORMAL']

    is.A <- (mt=='A')
    muts$tumor.reads[is.A] <- ((geno.data[['FAZ']][,'TUMOUR'][is.A] +  geno.data[['RAZ']][,'TUMOUR'][is.A]))
    muts$normal.reads[is.A] <- ((geno.data[['FAZ']][,'NORMAL'][is.A] +  geno.data[['RAZ']][,'NORMAL'][is.A]))
    muts$tumor.freq[is.A] <- ((geno.data[['FAZ']][,'TUMOUR'][is.A] +  geno.data[['RAZ']][,'TUMOUR'][is.A]))/ muts$tumor.depth[is.A] 
    muts$normal.freq[is.A] <- ((geno.data[['FAZ']][,'NORMAL'][is.A] +  geno.data[['RAZ']][,'NORMAL'][is.A]))/ muts$normal.depth[is.A]

    is.C <- (mt=='C')
    muts$tumor.reads[is.C] <- ((geno.data[['FCZ']][,'TUMOUR'][is.C] +  geno.data[['RCZ']][,'TUMOUR'][is.C]))
    muts$normal.reads[is.C] <- ((geno.data[['FCZ']][,'NORMAL'][is.C] +  geno.data[['RCZ']][,'NORMAL'][is.C]))
    muts$tumor.freq[is.C] <- ((geno.data[['FCZ']][,'TUMOUR'][is.C] +  geno.data[['RCZ']][,'TUMOUR'][is.C]))/ muts$tumor.depth[is.C] 
    muts$normal.freq[is.C] <- ((geno.data[['FCZ']][,'NORMAL'][is.C] +  geno.data[['RCZ']][,'NORMAL'][is.C]))/ muts$normal.depth[is.C]

    is.G <- (mt=='G')
    muts$tumor.reads[is.G] <- ((geno.data[['FGZ']][,'TUMOUR'][is.G] +  geno.data[['RGZ']][,'TUMOUR'][is.G]))
    muts$normal.reads[is.G] <- ((geno.data[['FGZ']][,'NORMAL'][is.G] +  geno.data[['RGZ']][,'NORMAL'][is.G]))
    muts$tumor.freq[is.G] <- ((geno.data[['FGZ']][,'TUMOUR'][is.G] +  geno.data[['RGZ']][,'TUMOUR'][is.G]))/ muts$tumor.depth[is.G] 
    muts$normal.freq[is.G] <- ((geno.data[['FGZ']][,'NORMAL'][is.G] +  geno.data[['RGZ']][,'NORMAL'][is.G]))/ muts$normal.depth[is.G]
    
    is.T <- (mt=='T')
    muts$tumor.reads[is.T] <- ((geno.data[['FTZ']][,'TUMOUR'][is.T] +  geno.data[['RTZ']][,'TUMOUR'][is.T]))
    muts$normal.reads[is.T] <- ((geno.data[['FTZ']][,'NORMAL'][is.T] +  geno.data[['RTZ']][,'NORMAL'][is.T]))
    muts$tumor.freq[is.T] <- ((geno.data[['FTZ']][,'TUMOUR'][is.T] +  geno.data[['RTZ']][,'TUMOUR'][is.T]))/ muts$tumor.depth[is.T] 
    muts$normal.freq[is.T] <- ((geno.data[['FTZ']][,'NORMAL'][is.T] +  geno.data[['RTZ']][,'NORMAL'][is.T]))/ muts$normal.depth[is.T]
}

# whether the sub was found to be a SNP
muts$isSnp <- info.data[,'SNP']

result<- list()
result$all.hist <- all.hist 
result$passed.hist <- passed.hist 
result$muts <- muts
result$mut.table <- mut.table
result$snps.hist <- snps.hist
result$non.snps.hist <- non.snps.hist
result

}
