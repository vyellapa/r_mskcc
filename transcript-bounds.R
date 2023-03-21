
##' Emperically determine the transcribed length of transcription boundaries.
##'
##' The transcription unit (keyed by `tx.key`) are summarized by their
##' min(start):max(end) combos.
##'
##' @param include NULL or a character string of columns to include from
##' \code{annotated.granges} and attach to the aggregated output, you might
##' want this to be c('entrez.id', 'symbol'), for instance
aggregateTxBounds <- function(annotated.granges, tx.key='entrez.id',
                              anno.keep=c('utr3', 'utr3*'),
                              include=NULL) {
  stopifnot(tx.key == 'entrez.id')
  dt <- as(annotated.granges, 'data.table')

  ## Doing this before we get going so a "valid" empty data.table can be
  ## returned when nrow(dt) == 0L
  keep.cols <- c('seqnames', 'strand', 'start', 'end', 'width', tx.key)
  if ('symbol' %in% colnames(dt)) {
    keep.cols <- c('symbol', keep.cols)
  }
  keep.cols <- unique(keep.cols)

  dumb <- lapply(keep.cols, function(x) NA)
  names(dumb) <- keep.cols
  dumb <- data.table(dumb)
  dumb[[tx.key]] <- 'dumb'

  dt <- subset(dt, !(is.null(dt[[tx.key]]) | is.na(dt[[tx.key]])))
  if (nrow(dt) == 0L) {
    return(dumb)
  }

  if (!is.null(anno.keep)) {
    dt <- subset(dt, exon.anno %in% anno.keep)
  }
  if (nrow(dt) == 0L) {
    return(dumb)
  }

  setkeyv(dt, unique(c('seqnames', 'strand', tx.key)))
  bounds <- dt[, by=key(dt), {
    bound.start <- min(start)
    bound.end <- max(end)
    list(start=bound.start, end=bound.end, width=bound.end - bound.start + 1L)
  }]

  ## bounds[, unique(keep.cols), with=FALSE]
  if (is.character(include)) {
    if (!inherits(annotated.granges, 'GRanges')) {
      annotated.granges <- tryCatch({
        as(annotated.granges, 'GRanges')
      }, error=function(e) NULL)
    }

    if (is.null(annotated.granges)) {
      warning("The annotated.granges object could not be xformed into a ",
              "GRanges object to retrieve the extra annos. This is being ",
              "ignored -- not extra annos are returned")
    } else {
      meta <- values(annotated.granges)
      include <- setdiff(include, names(bounds))
      include <- include[include %in% names(meta)]
      if (length(include) == 0) {
        warning("Columns to include were not found in values()")
      } else {
        meta <- as.data.frame(meta[, include, drop=FALSE])
        gbounds <- as(bounds, 'GRanges')
        xref <- findOverlaps(gbounds, annotated.granges, select="first")
        ## as.data.frame(as.list()) because a 1 column DataFrame doesn't
        ## turn into a normal data.frame well
        bounds <- cbind(as.data.frame(bounds), meta[xref,,drop=FALSE])
      }
    }
  }

  bounds
}

##' Combines annotated transcription bounds with observed transcription.
##'
##' The start of the transcription unit is taken from the annotations
##' The end is inferred from the observation. If extend only is TRUE,
##' then the end of the tx bounds can only be extended, otherwise the bounds
##' will be shrunk to the farthest observed transcription in `observed`, even
##' if it ends before the annotatd tx unit.
##'
##' NOTE: If you are trying to build extended 3'utrs, the txBounds for
##' non-coding transcripts (ones that have just 'utr' exon.anno's) will not
##' be present in the returned data.frame
##'
##' @param annotation An annotated GRanges-like object built from the annotated
##' genome.
##' @param observed An annotated GRanges-like object of the observed events
##' (reaods) that was built from annotateReads (or annotateExperiment).
##' @param extend.only Logical indicating whether to only change the annotated
##' transcription by extending it (\code{=TRUE})
##' @param tx.key The pimary-key like column(s) to use when aggregating
##' tx bounds.
##' @param anno.keep The annotated regions to keep and calculate bounds over.
##' @param fix If \code{start}: only the (strand-aware) ends of the tx units
##' are changed. \code{end}: Only the starts are changed. \code{none}: The
##' tx unit is set to the bounds of observed transcription.
##' @return A data.frame with the start/end boundaries.
inferTxBounds <- function(annotation, observed, extend.only=TRUE,
                          tx.key='entrez.id', anno.keep=c('utr3'),
                          observed.keep=c('utr3', 'utr3*'),
                          fix=c('start', 'end', 'none'),
                          include='symbol') {
  fix <- match.arg(fix)
  stopifnot(tx.key == 'entrez.id')

  o.bounds <- aggregateTxBounds(observed, tx.key, observed.keep,
                                include=include)
  if (fix == 'none') {
    return(o.bounds)
  }
  o.bounds <- as.data.table(o.bounds)
  setkeyv(o.bounds, c('seqnames', 'strand', tx.key))

  ## if (aggregate.annotated) {
  ##   a.bounds <- aggregateTxBounds(annotation, tx.key, anno.keep)
  ## }
  a.bounds <- aggregateTxBounds(annotation, tx.key, anno.keep)
  a.bounds <- as(a.bounds, 'data.table')
  setkeyv(a.bounds, key(o.bounds))

  tx.bounds <- merge(o.bounds, a.bounds, suffixes=c("", ".1"))
  tx.bounds <- as.data.frame(tx.bounds)
  is.rev <- tx.bounds$strand == '-'
  is.fwd <- !is.rev

  ## Start with starts and ends from the annotated data
  starts <- tx.bounds$start.1
  ends <- tx.bounds$end.1

  if (fix == 'start') {
    if (extend.only) {
      ## Extend to observed transcription, minimum length of tx.unit is the
      ## annotated one
      ends[is.fwd] <- pmax(tx.bounds$end[is.fwd], tx.bounds$end.1[is.fwd])
      starts[is.rev] <- pmin(tx.bounds$start[is.rev], tx.bounds$start.1[is.rev])
    } else {
      ends[is.fwd] <- tx.bounds$end[is.fwd]
      starts[is.rev] <- tx.bounds$start[is.rev]
    }
  } else if (fix == 'end') {
    ## fix == 'end'
    if (extend.only) {
      starts[is.fwd] <- pmin(tx.bounds$start[is.fwd], tx.bounds$start.1[is.fwd])
      ends[is.rev] <- pmax(tx.bounds$end[is.rev], tx.bounds$end.1[is.rev])
    } else {
      starts[is.fwd] <- tx.bounds$start[is.fwd]
      ends[is.rev] <- tx.bounds$end[is.rev]
    }
  } else {
    ## fix == 'none'
    stop("This should have already been taken care of")
  }

  tx.bounds$start <- starts
  tx.bounds$end <- ends

  tx.bounds$extended.by <- ifelse(is.fwd, tx.bounds$end - tx.bounds$end.1,
                                  tx.bounds$start.1 - tx.bounds$start)

  drop.cols <- grep("\\.1$", colnames(tx.bounds))
  if (length(drop.cols) > 0L) {
    tx.bounds <- tx.bounds[, -drop.cols, drop=FALSE]
  }

  ## Reset the width to the newly calculated starts/ends
  transform(tx.bounds, width=end - start + 1L)
}

test.inferTxBounds <- function() {
  annotated <- GRanges(c('chr1', 'chr1', 'chr2'),
                       IRanges(c(10, 50, 100), width=20),
                       c('+', '+', '-'),
                       entrez.id=c('100', '101', '102'))

  observed <- GRanges(c('chr1', 'chr1', 'chr2'),
                      IRanges(c(45, 60, 95), width=1L),
                      c('+', '+', '-'),
                      entrez.id=c('100', '101', '102'))
  itx <- inferTxBounds(annotated, observed, anno.keep=NULL, observed.keep=NULL)
  stopifnot(all(itx$start == c(10L, 50L, 95L)))
  stopifnot(all(itx$end == c(45L, 69L, 119L)))

  itx.shrink <- inferTxBounds(annotated, observed, extend.only=FALSE,
                              anno.keep=NULL, observed.keep=NULL)
  stopifnot(all(itx.shrink$start == c(10L, 50L, 95L)))
  stopifnot(all(itx.shrink$end == c(45L, 60L, 119L)))
}

