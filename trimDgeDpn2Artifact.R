setGeneric("trimDgeDpn2Artifact",
function(x, match.pattern='CGATC', trim=1L, ...) {
  standardGeneric("trimDgeDpn2Artifact")
})

setMethod("trimDgeDpn2Artifact", c(x="XStringSet"),
function(x, match.pattern, trim, .bounds.only=FALSE, ...) {
  verbose <- checkVerbose(...)
  m <- vmatchPattern(match.pattern,
                     narrow(x, start=1, width=nchar(match.pattern)))
  found <- countIndex(m) == 1L
  starts <- ifelse(found, 1L + trim, 1L)
  
  if (.bounds.only) {
    ret.val <- starts
  } else {
    ret.val <- narrow(x, start=starts)
  }

  if (verbose) {
    n.found <- sum(found)
    cat("Number of restriction sites found: ", n.found, "(",
        (n.found / length(found)) * 100, "% of reads)\n", sep="")
  }
  
  ret.val
})

setMethod("trimDgeDpn2Artifact", c(x="ShortReadQ"),
function(x, match.pattern, trim, ...) {
  starts <- trimDgeDpn2Artifact(sread(x), match.pattern=match.pattern, trim=trim,
                                .bounds.only=TRUE, ...)
  narrow(x, start=starts)
})

setMethod("trimDgeDpn2Artifact", c(x="character"),
function(x, match.pattern, trim, ...) {
  ## TODO: Trim a FASTQ file of DpnII artifacts?
  xt <- trimDgeDpn2Artifact(DNAStringSet(x), match.pattern=match.pattern,
                            trim=trim, ...)
  invisible(xt)
})
