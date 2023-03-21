require("Homo.sapiens")
x <- Homo.sapiens

################################################################################
## VERY simple minded tests for VERY simple wrapper funs
## Basically are they returning something?  If so then they are proabbly OK.
.testPopulatedGranges <- function(res){
    all(checkTrue(length(res) > 1),
        checkTrue(class(res) == "GRanges"))
}

test_microRNAs <- function(){
    res <- microRNAs(x)
    .testPopulatedGranges(res)
}

require("FDb.UCSC.tRNAs")
test_tRNAs <- function(){
    res <- tRNAs(x)
    .testPopulatedGranges(res)
}

test_promoters <- function(){
    res <- promoters(x)
    .testPopulatedGranges(res)
}

test_threeUTRsByTranscript <- function(){
    res <- threeUTRsByTranscript(x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "CompressedGRangesList")
}


test_asBED <- function(){
    res <- asBED(x)
    .testPopulatedGranges(res)
}

test_asGFF <- function(){
    res <- asGFF(x)
    .testPopulatedGranges(res)
}

test_disjointExons <- function(){
    res <- disjointExons(x)
    .testPopulatedGranges(res)
}

test_distance <- function(){
    gr <- GRanges(c("chr19", "chr20"),
                  IRanges(c(60000000, 40000000),  width=200))
    res <- distance(gr, x,id=c('1','100'),type="gene")
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "integer")
}


## This adds new suggests
require(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

test_extractTranscriptSeqs <- function(){
    res <- extractTranscriptSeqs(genome, x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "DNAStringSet")
}

test_extractUpstreamSeqs <- function(){
    res <- extractUpstreamSeqs(genome, x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "DNAStringSet")
}

test_intronsByTranscript <- function(){
    res <- intronsByTranscript(x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "CompressedGRangesList")
}

test_fiveUTRsByTranscript <- function(){
    res <- fiveUTRsByTranscript(x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "CompressedGRangesList")
}

test_isActiveSeq <- function(){
    ## getter works
    res <- isActiveSeq(x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "logical")
    checkTrue(all(res) == "TRUE")
    ## and check setter works
    isActiveSeq(x)[seqlevels(x)] <- FALSE
    res <- isActiveSeq(x)
    checkTrue(length(res) > 1)
    checkTrue(class(res) == "logical")
    checkTrue(all(res) == "FALSE")
    ## Because the tests run in some 'other' order, I need to set this
    ## back to default!
    isActiveSeq(x)[seqlevels(x)] <- TRUE
}

## TODO: OK I can't win. :/ This unit test proves that I really do
## require the storage of an actual TxDb in a slot (IOW just
## preferentially pulling these object from storage is not enough to
## solve all my problems since that will only resolve the read only
## issues).
## The reason is because if I just always get things from the source
## files then the resource effectively becomes a 'read only' object.
## So in order to have functioning activeSeq values I need to store
## the object dujour somewhere (and swap/change it when the TxDb<- is
## used)
## Changes to fix from here:
## 1) add slot to OrganismDb (maybe MultiDb instead)?  - done
## 2) .getTxDb needs to pull from that slot - done
## 3) .updateTxDb needs to put an actual TxDb into that slot - automatic
## 4) constructor for OrganismDb needs to also populate that slot (maybe tricky since it needs to know which thing is a TxDb - instantiate and check class ?) - done

## 5) All 'other' access cases (like select etc. are fine to use the loadDb 1st approach of getting the TxDb (And in fact this is preferred) - nothing to do but test?

## PROBLEM: the above only solves this for OrganismDb (not MultiDb) UNLESS I add a slot to MultiDb called what? TxDbIfPresent?

## IF I do that (add a slot for TxDb to MultiDb), then my instructions remain basically the same (except for 1&4) but it will work better overall.  I think I want to do it this way for TxDbs (because they should work correctly for MultiDbs too (without having leaky scope issues)


## BUT adding it to MultiDb means adding an Herve style Union class
## examples in: S4Vectors, grep for 'orNull'
