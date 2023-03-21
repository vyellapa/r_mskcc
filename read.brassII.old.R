read.brassII <- function(FILE.REARR) {


    rearrs <- try(suppressWarnings(read.table(FILE.REARR, header=FALSE, sep='\t', stringsAsFactors=F)))
    if (!inherits(rearrs, "try-error")) {
        rearrs$V1 <- as.character(rearrs$V1 )
        rearrs$V4 <- as.character(rearrs$V4 )
        rearrs$V1 [rearrs$V1=='23'] <- 'X'
        rearrs$V4[rearrs$V4=='24'] <- 'Y'
        cat(paste( dim(rearrs)[1], ' rearrs \n'))
        
        pf <- rep(4,nrow(rearrs))
        f.32 <- rearrs$V1!=rearrs$V4; pf[f.32] <- 32 # translocations
        f.1 <- rearrs$V1==rearrs$V4 & rearrs$V9=='+' & rearrs$V10=='-'; pf[f.1] <- 1 # inversion +- genomic strand
        f.8 <- rearrs$V1==rearrs$V4 & rearrs$V9=='-' & rearrs$V10=='+'; pf[f.8] <- 8 # inversion -+ genomic strand
        f.2 <- rearrs$V1==rearrs$V4 & rearrs$V9=='+' & rearrs$V10=='+'; pf[f.2] <- 2 # deletion ++ genomic strand
        rearrs$pf <- pf
        
        rearrs.formatted <- data.frame(Chromosome=rearrs$V1, chromStart=rearrs$V2, chromEnd=rearrs$V2, Chromosome.1=rearrs$V4, chromStart.1=rearrs$V5, chromEnd.1=rearrs$V5, pf=rearrs$pf, flag1=rearrs$V9, flag2=rearrs$V10)
        return(rearrs.formatted)
    } else {
        return(data.frame())
    }
    
    
}
