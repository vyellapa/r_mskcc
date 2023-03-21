## 
## Not really using these so anymore
## 
setValidity("CompressedReads", function(object) {
  meta <- values(object)
  if (!'count' %in% colnames(meta)) {
    return("`count` missing from elementMetadata")
  }
  TRUE
})

setGeneric("compressReads",
function(x, ...) {
  standardGeneric("compressReads")
})

setAs("CompressedReads", "GRanges",function(from) {
  if (length(from) == 0L) {
    return(GRanges())
  }
  ## times <- values(from)$count
  ## GRanges(seqnames=rep(seqnames(from), times),
  ##         ranges=IRanges(start=rep(start(from), times),
  ##                        end=rep(end(from), times)),
  ##         strand=rep(strand(from), times))
  class(from) <- "GRanges"
  from
})

setAs("CompressedReads", "IRanges", function(from) {
  if (length(from) == 0L) {
    return(IRanges())
  }
  times <- values(from)$count
  IRanges(start=rep(start(from), times), end=rep(end(from), times))
})

setAs("GRanges", "CompressedReads", function(from) {
  compressReads(from)
})

setAs("IRanges", "CompressedReads", function(from) {
  compresRanges(from)
})

##' Created a \code{\link{CompressedReads}} object from \code{GRanges} reads
##' by keeping only the unique ranges, and adding a $count value in
##' the \code{elementMetadata} for the number of times that unique range
##' appares in \code{x}
##'
##' @param x The reads (\code{GRanges}) object
setMethod("compressReads", c(x="GRanges"),
function(x, with.values=TRUE, ...) {
  chr.reads <- lapply(split(x, seqnames(x)), function(sreads) {
    if (length(sreads) == 0) {
      return(NULL)
    }
    seqname <- seqnames(sreads)[1]
    cstrand <- lapply(split(sreads, strand(sreads)), function(.reads) {
      iranges <- ranges(.reads)
      values(iranges) <- values(.reads)

      i <- compressRanges(iranges, with.values=with.values)
      meta <- elementMetadata(i)
      elementMetadata(i) <- NULL

      if (length(i) > 0L) {
        gr <- GRanges(seqnames=seqname, ranges=i, strand=strand(.reads)[1])
        values(gr) <- meta
        gr
      } else {
        NULL
      }
    })
    do.call(c, unname(cstrand[!sapply(cstrand, is.null)]))
  })
  cr <- do.call(c, unname(chr.reads[!sapply(chr.reads, is.null)]))
  class(cr) <- 'CompressedReads'
  cr
})

setMethod("compressReads", c(x="IRanges"),
function(x, ...) {
  compressRanges(x, ...)
})

compressRanges <- function(x, with.values=TRUE, ...) {
  stopifnot(inherits(x, 'IRanges'))
  if (length(x) == 0L) {
    i <- IRanges()
    values(i) <- DataFrame(count=rep(0L, 0))
  } else {
    x <- x[order(x)]
    i <- unique(x)
    orig.counts <- values(x)$count
    has.counts <- !is.null(orig.counts)

    sh <- subjectHits(findOverlaps(x, i, type='equal'))

    if (!(any(duplicated(sh)))) {
      ## Nothing has changed
      if (!has.counts) {
        values(x)$count <- 1L
      }
      return(x)
    }

    if (has.counts) {
      orig.count <- dlply(idata.frame(data.frame(idx=sh, count=orig.counts)),
                          .(idx), summarise, count=sum(count))
      counts <- unlist(unname(orig.count), use.names=FALSE)
    } else {
      counts <- runLength(Rle(sh))
    }

    if (with.values) {
      values(i)$count <- counts
    } else {
      values(i) <- DataFrame(count=counts)
    }
  }
  i
}

setGeneric("uncompress", function(x, ...) standardGeneric("uncompress"))
setMethod("uncompress", c(x="CompressedReads"),
function(x, with.values=FALSE, ...) {
  f <- selectMethod('uncompress', 'GRanges')
  x <- f(x, with.values=FALSE, ...)
  if (!is.null(values(x)$tpm)) {
    values(x)$tpm <- values(x)$tpm / values(x)$count
  }
  values(x)$count <- 1L
  x
})

setMethod("uncompress", c(x="GRanges"),
function(x, with.values=FALSE, ...) {
  count <- values(x)$count
  if (is.numeric(count)) {
    x <- rep(x, count)
  }
  x
})


