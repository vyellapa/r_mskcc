setGeneric("getBsGenome", function(x, ...) standardGeneric("getBsGenome"))

##'  Returns the name of the aligner used
##' @param x A \code{\linkS4class{SeqStore}} object
##' @return The name of the aligner used, eg. \code{bwa}, \code{bowtie}, etc.
setGeneric("aligner", function(x, ...) standardGeneric("aligner"))


##' Query \code{\linkS4class{SeqStore}} object to see if data is from a paired end run
##'
##' @exportMethod isPaired
##'
##' @param x A \code{\linkS4class{SeqStore}} object
##' @return A logical indicating if data is paired.
setGeneric("isPaired", function(x, ...) standardGeneric("isPaired"))

setGeneric("header", function(x, ...) standardGeneric("header"))

setGeneric("swapStrand", function(x, ...) standardGeneric("swapStrand"))



##' @exportClass TagSeqExperiment
setClass("TagSeq", contains="SummarizedExperiment")

##' @exportClass SAGEseqExperiment
setClass("SAGEseq", contains="TagSeq")

##' @exportClass CompressedReads
setClass("CompressedReads", contains="GRanges")

## An (annotated) compressed container for the tags (per chromosome).
## The values() DataFrame will store:
##  * strand (factor)
##  * chr (Rle<character>)
##  * is.unique (Rle<logical>)
##  * tag.type (tagType factor)
##  * symbol (character)
##  * tag.type.anti (character) -- the annotation on the opposite strand
##  * symbol.anti (character) -- the annotated gene on the opposite strand
##  * tag.type.extra -- add 1,2,3 for splitting of utrs?
## setClass("AnnotatedTags", contains="TagSequences",
##          representation(.tags='GRanges'),
##          prototype(.tags=GRanges()))

## -----------------------------------------------------------------------------
## Methods


