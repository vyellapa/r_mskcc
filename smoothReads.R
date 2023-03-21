################################################################################
## Smoothers to clean out tag-sequence-like data
##
## Test Data
## gr <- GRanges(seqnames='chr1', strand='+',
##               ranges=IRanges(start=c(sample(18:22, 10, replace=TRUE),
##                                      sample(100:108, 10, replace=TRUE)),
##                              width=sample(18:20, 20, replace=TRUE)))

##' Locates islands/clusters of tags and removes any such island whose mean
##' coverage does not pass the \code{min.cover} threshold
smoother.islandCoverage <- function(reads, min.cover=5, window=9,
                                    honor.strand=TRUE, ...) {
  if (!is.numeric(min.cover)) {
    warning("Illegal value for smoother.idlandCoverage::min.cover, set to 1")
    min.cover <- 1L
  }

  if (honor.strand) {
    ## A little recursion. Each iteration will drop down into the code block
    ## below.
    reads <- split(reads, strand(reads))
    smoothed <- seqapply(reads, smoother.islandCoverage, min.cover, window,
                         honor.strand=FALSE)
    smoothed <- unlist(smoothed, use.names=FALSE)
    if (length(smoothed) > 1) {
      smoothed <- smoothed[order(start(smoothed))]
    }
    return(smoothed)
  }


  if (length(reads) == 0) return(reads)

  covr <- coverage(reads)[[1]]
  ## covr.mean <- runmean(covr, window, 'constant')
  ## bounds <- slice(covr.mean, lower=min.cover)
  bounds <- slice(covr, lower=min.cover)

  if (length(bounds) != 0) {
    seq.names <- seqnames(reads[1])
    b.strand <- '*'
  } else {
    seq.names <- character()
    b.strand <- strand()
  }

  take <- GRanges(seqnames=seq.names, ranges=bounds, strand=b.strand)
  subsetByOverlaps(reads, take)
}


##' Converts tags that are presumably the same, but have "jagged ends" to all
##' have uniform start/end positions.
##'
##' This assumes that all tags that *even slightly* overlap eachother are meant
##' to be the same.
smoother.repackTags <- function(reads, honor.strand=TRUE, ...) {
  .reads <- reads
  ## reduce()-ing the reads nukes the metadata. If we are repacking a
  ## CompressedReads object, the validObject function will complain
  class(.reads) <- 'GRanges'
  if (!honor.strand) {
    strand(.reads) <- '*'
  }
  fence.posts <- reduce(.reads)
  ranges(reads) <- ranges(fence.posts)[match(.reads, fence.posts)]
  reads
}

##' Remove reads that don't have correct anchor/restrition site on 5' ends.
##'
##' This can happen when alignments can have some mismatches
smoother.anchorSite <- function(reads, bs.chr, anchor.site,
                                anchor.max.mismatch=0L, as.logical=FALSE, ...) {
  if (is.null(anchor.site) || nchar(anchor.site) == 0L) {
    return(if (as.logical) rep(TRUE, length(reads)) else reads)
  }

  seqs <- Views(bs.chr, ranges(reads))
  is.rev <- strand(reads) == '-'
  keep <- rep(FALSE, length(seqs))
  strands <- as.character(unique(strand(reads)))
  anchor.site <- lapply(anchor.site, DNAString)

  for (.strand in strands) {
    strand.take <- as.logical(strand(reads) == .strand)
    for (rs in anchor.site) {
      if (.strand == '-') {
        rs <- reverseComplement(rs)
      }
      fix <- if (.strand == '-') 'end' else 'start'
      .seqs <- resize(seqs[strand.take], width=nchar(rs), fix=fix)
      found <- vcountPattern(rs, .seqs, max.mismatch=anchor.max.mismatch)
      keep[strand.take] <- keep[strand.take] | as.logical(found)
    }
  }

  if (as.logical) {
    keep
  } else {
    reads[keep]
  }
}
