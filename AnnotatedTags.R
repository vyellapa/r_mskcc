## setMethod("compressReads", c(x="GRanges"),
## function(x, ...) {
##   pieces <- lapply(split(x, seqnames(x)), function(piece) {
##     strands <- lapply(split(piece, strand(piece)), function(sp) {
##       sp <- sp[order(start(sp))]
##       .granges <- reduce(sp)
##       .rle <- Rle(start(sp))

##       if (!all(start(.granges) == runValue(.rle))) {
##         stop("Compression broken: ", seqnames(sp[1]), ", strand ",
##              strand(sp[1]))
##       }
##       values(.granges) <- DataFrame(count=runLength(.rle))
##       .granges
##     })
##     names(strands) <- NULL
##     do.call(c, strands)
##   })
##   names(pieces) <- NULL
##   compressed <- do.call(c, pieces)
##   class(compressed) <- 'CompressedReads'
##   compressed
## })
 


## .required.AnnotatedTags.meta <- c(
##   count='integer',
##   strand='factor',
##   chr='Rle',
##   is.unique='logical',
##   tag.type='Rle',
##   symbol='character',
##   tag.type.anti='Rle',
##   symbol.anti='Rle',
##   tag.type.extra='integer')

## .emptyAnnotatedTagsDataFrame <- function() {
##   df <- DataFrame()
##   for (what in names(.required.AnnotatedTags.meta)) {
##     df[[what]] <- new(.required.AnnotatedTags.meta[[what]])
##   }
##   df
## }

## setMethod("initialize", "AnnotatedTags", 
## function(.Object, ...,
##          .tags=GRanges()) {
##   # if (ncol(values(.tags)) == 0L) {
##   #   values(.tags) <- 
##   # }
##   callNextMethod(.Object, .tags=.tags, ...)
## })

## ##' Build an annotated tags object
## AnnotatedTags <- function(ranges, meta=NULL, ...) {
##   if (!is.null(meta)) {
##     values(ranges) <- meta
##   }
  
##   new('AnnotatedTags', .tags=ranges)
## }

## ##' Ensures that there is an DataFrame attached to the GRanges object containg
## ##' the appropriate columns
## setValidity("AnnotatedTags", function(object) {
##   ## names() are the column names required
##   ## values are the type of object
##   meta.required <- .required.AnnotatedTags.meta  
##   meta <- values(object)
##   meta.missing <- !(names(meta.required) %in% names(meta))
##   if (any(meta.missing)) {
##     return(paste(paste(names(meta.required[meta.missing]), collapse=","),
##                  "are missing from values()", sep=" "))
##   }
  
##   for (what in names(meta.required)) {
##     if (!inherits(meta[[what]], meta.required[[what]])) {
##       return(paste("Illegal type for ", what, ". ", meta.required[[what]],
##                    "was expected, but got ", is(meta[[what]])[1], sep=""))
##     }
##   }
  
##   TRUE
## })

## setAs('AnnotatedTags', 'GRanges', function(from) {
##   ## TODO: Include more meta information in the expanded AnnotatedTags object
##   meta <- values(from)
##   snames <- rep(as.vector(meta$chr), meta$count)
##   starts <- rep(start(from), meta$count)
##   ends <- rep(end(from), meta$count)
##   strands <- rep(as.factor(meta$strand), meta$count)
##   GRanges(seqnames=snames, ranges=IRanges(starts, ends), strand=strands)
## })

## setMethod("values", c(x="AnnotatedTags"),
## function(x, ...) {
##   elementMetadata(x@.tags, ...)
## })

## setMethod("elementMetadata", c(x="AnnotatedTags"),
## function(x, ...) {
##   elementMetadata(x@.tags, ...)
## })

## setMethod("ranges", c(x="AnnotatedTags"),
## function(x, ...) {
##   ranges(x@.tags)
## })

## setMethod("start", c(x="AnnotatedTags"),
## function(x, ...) {
##   start(x@.tags)
## })

## setMethod("end", c(x="AnnotatedTags"),
## function(x, ...) {
##   end(x@.tags)
## })

## setMethod("width", c(x="AnnotatedTags"),
## function(x) {
##   width(x@.tags)
## })
