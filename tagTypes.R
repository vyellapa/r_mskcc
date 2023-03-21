##' Returns the factors used in the tag annotation table.
##'
##' name      
##' ----------
##' utr       
##' utr5      
##' utr5*     (inferred utr5 extension)
##' cds       
##' utr3      
##' utr3*     (inffered utr3 extension)
##' intron    
##' overlap   
##' intergenic

setGeneric("tagType", function(x, ...) standardGeneric("tagType"))
setMethod("tagType", c(x="missing"),
function(x, ...) {
  factor(levels=c('utr', 'utr5', 'utr5*', 'cds', 'utr3', 'utr3*', 'intron',
                  'overlap', 'intergenic'))
})

setMethod('tagType', c(x='character'),
function(x, ...) {
  factor(x, levels=levels(tagType()))
})
 
## setMethod("tagType")
## tagType <- function()
## ## Default method for generic
## getTagTypeTable <- function(x, ...) {
##   ## The order these are generated as is important!
##   tags <- c('utr', 'utr5', 'cds', 'utr3', 'intron', 'overlap',
##             'out-of-bounds-3', 'out-of-bounds-5', 'intergenic')
##   data.frame(id=seq(tags), name=tags, stringsAsFactors=FALSE)
## }
## setGeneric("getTagTypeTable")
## 
## setMethod("getTagTypeTable", c(x="DGEseq"),
## function(x, ...) {
##   warning("We should, perhaps, have experiment-specific behavior here, ",
##       "but we don't")
##   callNextMethod()
## })
