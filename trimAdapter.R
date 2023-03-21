##' Removes 3' adapter contamination from reads
##' 
##' If the sequence does not have a sufficient (w.r.t score.threshold) match to
##' the adapter, then it is returned unchanged.
##' 
##' (This function was adapted from girafe version 1.1.5)
##' 
##' TODO: Check effect of tweaking score.threshold
##' 
##' @param x A FASTQ file, DNAStringSet,  or ShortRead object of reads
##' @param adapter The character/DNAString representing the adapter to strip
##' @param match.score Score to add when characters match in pairwiseAlignment
##' @param mismatch.score Score to substract when mismatch
##' @param score.threshold Sensitivity of matching against adapter
setGeneric("trimAdapter", signature=c("x"),
function(x, adapter='TCGTATGCCGTCTTCTGCTTG', match.score=1, mismatch.score=2,
         score.threshold=2, gap.penalty=-5, gap.extension.penalty=-2, ...) {
  standardGeneric("trimAdapter")
})

##' @param .bounds.only For internal use only. If set to TRUE, then the function
##' is being used as "an engine" from a calling function to get the start/end
##' bounds to use as bounds for the "real" reads at eachposition, which is
##' what is returned.
setMethod("trimAdapter", c(x="XStringSet"),
function(x, adapter, match.score, mismatch.score, score.threshold, gap.penalty,
         gap.extension.penalty, .bounds.only=FALSE, ...) {
  verbose <- checkVerbose(...)
  if (is.character(adapter)) {
    adapter <- DNAString(adapter)
  }
  if (!inherits(adapter, "XString")) {
    stop("Unknown adapter type")
  }

  read.length <- unique(width(x))
  
  if (length(read.length) > 1) {
    stop(paste("Expected all reads in object",
               parse(substitute(x)),
               "to be of the same lengths! Found lengths:",
               paste(read.length, collapse=", "),"\n"))
  }
  
  mat <- nucleotideSubstitutionMatrix(match=match.score,
                                      mismatch=mismatch.score)
  
  ## giraffe uses gapOpening=0, gapExtension=-Inf
  pa <- pairwiseAlignment(x, adapter, type="overlap",
                          substitutionMatrix=mat,
                          gapOpening=gap.penalty,
                          gapExtension=gap.extension.penalty)
  
  areCompleteOverlap <- (score(pa) > score.threshold) &
                        (start(pattern(pa)) == 1) &
                        (end(pattern(pa)) == read.length)

  kstarts <- integer(length(x)) + 1L

  kends <- ifelse(score(pa) < score.threshold, read.length,
                  ifelse(end(pattern(pa)) == read.length,
                         end(pattern(pa)) - width(pattern(pa)), read.length))
  
  if (sum(areCompleteOverlap) > 0) {
    kstarts[areCompleteOverlap] <- 1L
    kends[areCompleteOverlap] <- read.length
  }

  if (length(kstarts) != length(x)) {
    stop("The length of the resulting boundaries are not same length as ",
         "original object")
  }
  
  if (.bounds.only) {
    x2 <- list(start=kstarts, end=kends)
  } else {
    x2 <- narrow(x, start=kstarts, end=kends)
  }

  if (verbose) {
    found <- kends < read.length
    n.found <- sum(found)
    cat("Adapter found in ", n.found, "reads ",
        "(", n.found / length(kends), "%)\n", sep="")
  }
  
  x2
})

setMethod("trimAdapter", c(x="ShortReadQ"),
function(x, adapter, match.score, mismatch.score, score.threshold, gap.penalty,
         gap.extension.penalty, ...) {
  bounds <- trimAdapter(sread(x), adapter=adapter, match.score=match.score,
                        mismatch.score=mismatch.score,
                        score.threshold=score.threshold,
                        gap.penalty=gap.penalty,
                        gap.extension.penalty=gap.extension.penalty,
                        .bounds.only=TRUE, ...)
  narrow(x, start=bounds$start, end=bounds$end)
})


## Loads a ShortReadQ file from x, or attempts to turn a character vector into
## a DNAStringSet to nuke adapters from
##
## Optionally writes result to a file (use for writing trimmed FASTQ files).
setMethod("trimAdapter", c(x="character"),
function(x, adapter, match.score, mismatch.score, score.threshold, gap.penalty,
         gap.extension.penalty,  out.file=NULL,
         qualityType=c("SFastqQuality", "FastqQuality", "Auto" ), ...) {         
  verbose <- checkVerbose(...)
  if (length(x) == 1 && file.exists(x)) {
    ## Read as ShortReadQ and fire away
    if (verbose) {
      cat("Reading FASTQ file ... ")
      ct <- proc.time()['elapsed']
    }
    qualityType <- match.arg(qualityType)
    x <- readFastq(dirname(x), basename(x), qualityType=qualityType)
    if (verbose) {
      cat("[", proc.time()['elapsed'] - ct, "]\n", sep="")
    }
  } else {
    x <- tryCatch(DNAStringSet(x), error=function(e) NA)
    if (!is(x, 'DNAStringSet')) {
      stop("Do not know how to handle the character argument. ",
           "If it is a file, it can't be read:\n",
           x[1])
    }
  }
  
  xt <- trimAdapter(x, adpater=adapter, match.score=match.score,
                    mismatch.score=mismatch.score,
                    score.threshold=score.threshold,
                    gap.penalty=gap.penalty,
                    gap.extension.penalty=gap.extension.penalty, ...)

  ## write back to file?
  if (!is.null(out.file)) {
    if (verbose) {
      cat("Writing to FASTQ file `", out.file, "` ... ", sep="")
      ct <- proc.time()['elapsed']
    }
    writeFastq(xt, out.file)
    if (verbose) {
      cat("[", proc.time()['elapsed'] - ct, "]\n", sep="")
    }
  }
  
  invisible(xt)
})
