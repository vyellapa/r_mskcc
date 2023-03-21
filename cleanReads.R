## 
## These should be moved to seqtools, if kept at all
##

##' The sequences have a run of 5' Ts we need to strip
##'
##' @param x The sequences to clean
##' @param min.nonT The minimum number of non-T characters in the sequence
##' to even consider for cleanning -- reads that don't pass this test are
##' removed right away.
setGeneric("cleanReads", function(x, min.nonT=4, min.length=10, ...) {
  standardGeneric("cleanReads")
})

## The ShortRead version should do something with the quality scores?
setMethod("cleanReads", c(x="ShortReadQ"),
function(x, min.nonT, min.length, ...) {
  ## bounds <- cleanReads(sread(x), min.nonT, min.length, bounds.only=TRUE, ...)
  ## narrow(x, bounds$start, bounds$end)
  fn <- selectMethod("cleanReads", c(x="DNAStringSet"))
  fn(x, min.nonT=min.nonT, min.length=min.length, ...)
})

## This function can actually work on a DNAStringSet or a ShortReadQ object
setMethod("cleanReads", c(x="DNAStringSet"),
function(x, min.nonT, min.length, bounds.only=FALSE, ...) {
  ## The initial reads had several T-{mostly|only} reads, so we only keep
  ## reads that have more than 4 non-T characters
  verbose <- checkVerbose(...)
  reads <- if (inherits(x, 'ShortReadQ')) sread(x) else x
  
  freqs <- alphabetFrequency(reads)
  not.enough.info <- rowSums(freqs[, c('A', 'C', 'G')]) < min.nonT
  
  if (verbose) {
    n.axe <- sum(not.enough.info)
    cat("Removing ", n.axe, " TTTT reads [", n.axe / length(x), "%]\n",
        sep="")
  }
  
  ## We use a simple regular expression to nuke a run of TTTNTT*'s
  ## (probably not the best, but effective nonetheless)
  strip.ts <- regexpr("^[TN]+", reads, perl=TRUE)
  t.length <- attr(strip.ts, 'match.length')
  
  starts <- ifelse(t.length == -1, 1, t.length + 1)
  starts <- ifelse(starts > width(reads), width(reads), starts)
  ends <- width(reads)

  too.short <- rep(FALSE, length(x))
  if (!is.null(min.length)) {
    too.short <- (ends - starts  + 1) < min.length
    n.axe <- sum(too.short)
    if (verbose) {
      cat(n.axe, " reads are too short\n")
    }
  }

  ## Many reads have lots of N's -- axe any read with N
  still.n <- freqs[, 'N'] > 0
  
  axe <- not.enough.info | too.short | still.n
  if (verbose) {
    n.axe <- sum(axe)
    cat(n.axe, " reads do not meet criteria ",
        "[", n.axe / length(axe), "%]\n", sep="")
  }
  
  if (bounds.only) {
    ret.val <- list(start=starts, ends=ends, axe=axe)
  } else {
    ## Cut short sequences
    keep <- !axe
    k.starts <- starts[keep]
    k.ends <- ends[keep]
    ret.val <- narrow(x[keep], k.starts, k.ends)
  }
  
  ret.val
})
