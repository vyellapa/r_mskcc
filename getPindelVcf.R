getPindelVcf <- function(myFile, onlyPASSED=FALSE) {

                                        # plots mutation-context for all variants in the vcf file
                                        # and separately for the variants that passed


    
B35.1.2.vcf <- readVcf(myFile, "hg19")

# matrix with the info to each column
info(B35.1.2.vcf)

#filters failed for each variant
rd <- rowData(B35.1.2.vcf)


fs <- fixed(B35.1.2.vcf )$FILTER
fs.passed <- (fs=='PASS')

if (onlyPASSED) {
        
    B35.1.2.vcf <- B35.1.2.vcf[fs.passed,]
    rd <- rowData(B35.1.2.vcf)
    fs <- fixed(B35.1.2.vcf )$FILTER
    fs.passed <- (fs=='PASS')
    
}

inf <- info(B35.1.2.vcf)

rgs <- ranges(rd)
starts <- start(rgs)
ends <-  end(rgs)
chroms <- paste('chr', seqnames(rd), sep='')

fxd <- (fixed(B35.1.2.vcf))
wt <- as.character(fxd$REF)
mt <- as.character(unlist(fxd$ALT))

muts <- data.frame(Chromosome=chroms, chromStart=starts, chromEnd = ends, Type=inf$VT, wt=wt, mt=mt, pass=fs.passed)



result<- list()
result$muts <- muts
result

}
