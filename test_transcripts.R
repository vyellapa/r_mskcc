require("Homo.sapiens")

## TODO: to speed up unit tests I need some small pieces of annotation
## for testing a small subset of granges and the matching metadata
## Ideally I want a Homo.sapiens package that uses the small subset DB
## from GenomicFeatures.  Perhaps a "mini-me" package for testing? -
## but making this is a bit of a project.

x <- Homo.sapiens
txdb <- OrganismDbi:::.getTxDb(x)


## some internal testing (make sure helpers work as expected)

test_compressMetadata <- function(){
    cols <- c("SYMBOL","GENENAME", "TXCHROM", "PMID")
    txs <- transcripts(txdb, columns="tx_id")[1:100]  ## shortened
    meta <- select(x, keys=as.character(mcols(txs)$tx_id), cols, "TXID") 
    f <- factor(meta[["TXID"]],levels=mcols(txs)[["tx_id"]])
    res <- OrganismDbi:::.compressMetadata(f, meta, "TXID")
    checkTrue(is(res, "DataFrame"))
    checkTrue(dim(res)[2] ==4)
    checkTrue(dim(res)[1] ==100)
    checkTrue(all(colnames(res) %in% cols))
}

## .combineMetadata is an important helper function.
test_combineMetadata <- function(){
    cols <- c("SYMBOL","GENENAME", "TXCHROM", "PMID")
    txs <- transcripts(txdb, columns="tx_id")[1:100]  ## shortened
    meta <- select(x, keys=as.character(mcols(txs)$tx_id), cols, "TXID") 
    res <- OrganismDbi:::.combineMetadata(txs,meta,avoidID="TXID",
                                          joinID="tx_id", columns=cols)
    checkTrue(is(res, "DataFrame"))
    checkTrue(dim(res)[2] ==4)
    checkTrue(dim(res)[1] ==100)
    checkTrue(all(colnames(res) %in% c(cols))) 
}


## These tests are slow so I will need a smaller thing to test with...
test_transcripts <- function(){
    library(Homo.sapiens); h = Homo.sapiens; cols = c("TXNAME","SYMBOL")
    res <- transcripts(h, columns=cols)
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) > 80000)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("TXNAME","SYMBOL")))
}

test_exons <- function(){
    library(Homo.sapiens); h = Homo.sapiens; cols = c("TXCHROM","REFSEQ")
    res <- exons(h, columns=cols)
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) > 200000)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("TXCHROM","REFSEQ")))
}

test_cds <- function(){
    library(Homo.sapiens); h = Homo.sapiens; cols = c("GENENAME","SYMBOL")
    res <- cds(h, columns=cols)
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) > 200000)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))
}





test_transcriptsBy <- function(){
    library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
    res <- transcriptsBy(h, by="gene", cols)    
    checkTrue(class(res) == "CompressedGRangesList")
    checkTrue(length(res) > 20000)
    ## check inner mcols
    checkTrue(all(colnames(mcols(res[[1]])) %in%
                  c("tx_name","GENENAME","SYMBOL")))
    ## check outer mcols
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENEID","GENENAME","SYMBOL")))

    ## extra check for case where we only have one field.
    cols = c("SYMBOL")
    res2 <- transcriptsBy(h, by="gene", cols)
    checkTrue(class(res) == "CompressedGRangesList")
    checkTrue(length(res) > 20000)
    ## check inner mcols
    checkTrue(all(colnames(mcols(res2[[1]])) %in%
                  c("tx_name","SYMBOL")))
    ## outer mcols
    checkTrue(all(colnames(mcols(res2)) %in%
                  c("GENEID","SYMBOL")))
    ## 
}

test_exonsBy <- function(){
    library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
    res <- exonsBy(h, by="gene", cols)
    ## TODO: look more closely at this one.  The metadata looks off...
    checkTrue(class(res) == "CompressedGRangesList")
    checkTrue(length(res) > 20000) 
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))
}

test_cdsBy <- function(){
    library(Homo.sapiens);h=Homo.sapiens;by="gene";cols = c("GENENAME","SYMBOL")
    res <- cdsBy(h, by="gene", cols)
    checkTrue(class(res) == "CompressedGRangesList")
    checkTrue(length(res) > 19000)
    checkTrue(all(colnames(mcols(res)) %in%
                  c("GENENAME","SYMBOL")))
}




#################
## This test adresses a bug that relates to the access of TXID and
## other column values that were used for "joining" data together.

test_rangeMethods_for_JoinFailures <- function(){
    library(Homo.sapiens);h<-Homo.sapiens; cols<-"TXID"    
    res <- transcripts(h, columns=cols)
    checkTrue(names(mcols(res))=="TXID")
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) > 10000)  ## large
    
    res <- transcriptsBy(h, by="gene", columns=cols)
    checkTrue("TXID" %in% names(mcols(res[[1]])))
    checkTrue(class(res) == "CompressedGRangesList")
    checkTrue(length(res) > 10000)  ## large
}



## Other Bug: for exons() 'filter' argument does not work on Homo.sapiens objects...

test_filterArg <- function(){
    ## this works:
    res <- exons(Homo.sapiens, filter=list(gene_id="1"))
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) < 20)  ## small

    ## so does this! 
    res <- exons(Homo.sapiens, columns="SYMBOL", filter=list(gene_id="1"))
    checkTrue("SYMBOL" %in% names(mcols(res)))
    checkTrue(class(res) == "GRanges")
    checkTrue(length(res) < 20)  ## small
}


## Fast testing: BiocGenerics:::testPackage(pattern="^test_transcripts.*\\.R$")
