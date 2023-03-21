##' This is the representation of the DGE protocol developed by Markus Hafner.
##' The Anchor Enzyme refers to the enzyme (NlaIII, DpnII) used to make the
##' initial cut from the transcript tag.
##' Constructor for a DGEseq experiment
##'
##' @param directory The directory the alignments and meta data are in.
##' @param temp.store Sets the temp store for the internal database
##' @param with.db A logical indicating whether there is an associated
##' db that stores metadata
##' @param ... Variables to override any of the slot values for the
##' object that's created.
##'
##' @return An object containing information about the experiment.
SAGEseq <- function(directory, with.db=FALSE, ...) {
  ## Check for bam files and ignore any from restriction-info (lumiDGE)
  ss <- SeqStore(directory, 'SAGEseq', is.stranded=TRUE,
                 is.paired=FALSE, ...)
  args <- list(...)
  meta <- metadata(ss)
  restriction.site <- meta$restriction.site
  if (is.null(restriction.site)) {
    if (is.null(args$restriction.site)) {
      stop("Can't find restriction site in metadata.",
           "You can pass it in to function")
    }
    restriction.site <- args$restriction.site
  }
  restriction.site <- unlist(strsplit(meta$restriction.site, ","))
  restriction.site <- gsub("\\s+", "", restriction.site)
  
  ## Check that valid restriction sites are provided.
  ## TODO: The value of the restriction sites should be stored in some
  ##       meta file in the experiment SeqStore directory itself.
  rs.keep <- sapply(restriction.site, function(rs) nchar(rs) > 0)
  restriction.site <- restriction.site[rs.keep]
  if (length(restriction.site) == 0) {
    stop("Unknown restriction sites for experiment")
  }
  
  ss@restriction.site <- restriction.site
  ss
}

setGeneric("restrictionSite", function(x, ...) {
  standardGeneric("restrictionSite")
})

setMethod("restrictionSite", c(x="SAGEseq"),
function(x, ...) {
  x@restriction.site
})

setGeneric("anchorSite", function(x, ...) standardGeneric("anchorSite"))
setMethod("anchorSite", c(x="SAGEseq"),
function(x, ...) {
  x@restriction.site
})
