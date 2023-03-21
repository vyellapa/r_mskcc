## TODO: Improve internal priming score
##
## Using hump-finding-axe, we flag internal priming events when
## 8/10 in a sliding window are A's starting from peak of event to its
## 3'end.
##
## Examine:
##   1. Where do polya signal lie w.r.t to the peak of the event
##      - the 3' end of the event?
##      (use only events that !is.primed)
##   2. Look at A-enrichment (avg. of slidigin window A count) in
##      a. event
##      b. peak of event to 3'
##
## Looks like internal priming, but saved by polyA rescue:
##  * chr21:33,640,671-33,640,784
##    The 3' end is a very rich-A tail
##    The polyA signal (AGTAAA) is in a low-coverage area under the 3' tail
##    And its antisense!
##  * The peak immediately 5' (chr21:33,640,882-33,640,992) also looks like
##    internal priming. Also very A rich region.
##    Wide event (97bp) and low coverage.
##
## Older::
## This looks like internal priming but has a max of 7/10 A's in a row
##   chr1:93,730,436-93,730,504
## Think about creating a longer range score -- maybe it incorporates the
## avg # of A's compared to background for these longer events?
## Engage the brain and do something smarter, dang it!
##
## It really seems that "long/wide" events are due to internal priming!
## Look at antisenes reads here: chr1:182,351,914-182,352,666

## Some (older) areas that are tricky:
##  chr1:59,041,067-59,041,309
##    looks like 5' tail of event is due to internal priming, but 3' end
##    is real
##  chr1:45,244,272-45,244,514
##    5' end of event is internally primed?

## TODO: events +/- sliding window distance can overrun (or under run) the
##       start and end of chromosomes -- fix!
internalPrimingSummary <- function(events, bsg, ip.functions=NULL, chrs=NULL,
                                   .parallel=TRUE, ...) {
  if (is.null(chrs)) {
    chrs <- seqlevels(events)
  }
  if (is.null(ip.functions)) {
    ip.functions <- list(
      externally.primed=flagInternallyPrimedEventsByBoundedQuantiles,
      internally.primed=flagInternalInternallyPrimedEventsByBoundedQuantiles
    )
  }

  if (.parallel && length(chrs) > 1L) {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }

  ip.args <- list(...)
  processed <- foreach(chr=chrs, .packages=c('Biostrings'),
                       .options.multicore=list(preschedule=FALSE)) %loop% {
    idxs <- which(as.logical(seqnames(events) == chr))
    if (length(idxs) == 0L) {
      return(NULL)
    }
    these <- events[idxs]
    dna <- unmasked(bsg[[chr]])
    ips <- lapply(ip.functions, function(ip.function) {
      do.call(ip.function, c(list(x=these, dna=dna), ip.args))
    })
    names(ips) <- names(ip.functions)
    c(list(idxs=idxs), ips)
  }
  names(processed) <- chrs
  processed <- processed[!sapply(processed, is.null)]

  ## Generate a named-result-list
  ## Figure out how many columns we need, and their names
  col.names <- unlist(sapply(names(ip.functions), function(fname) {
    paste(fname, names(processed[[1]][[fname]]), sep=".")
  }))

  result <- lapply(col.names, function(cn) rep(0L, length(events)))
  names(result) <- col.names
  result <- as(result, 'DataFrame')

  for (chr.result in processed) {
    idxs <- chr.result$idxs
    for (ipf.name in names(ip.functions)) {
      for (rname in names(chr.result[[ipf.name]])) {
        ans <- chr.result[[ipf.name]][[rname]]
        col.name <- paste(ipf.name, rname, sep=".")
        result[idxs, col.name] <- ans
      }
    }
  }

  result
}

##' @param quantile.offset If 'range', then the numbers in quantile.{start|end}
##' are from WITHIN the range of x
.reconstructQuantileBounds <- function(x, quantile.start, quantile.end,
                                       quantile.offset=c('range', 'none'),
                                       ...) {
  quantile.offset <- match.arg(quantile.offset)

  ## Get/check quantile info
  if (is.character(quantile.start)) {
    quantile.start <- values(x)[[quantile.start]]
    if (is.null(quantile.start)) {
      stop("Can't find quantile.starts using `quantile.start` param")
    }
  }
  if (is.character(quantile.end)) {
    quantile.end <- values(x)[[quantile.end]]
    if (is.null(quantile.end)) {
     stop("Can't find quantile.end using `quantile.end` param")
    }
  }
  if (!is.numeric(quantile.start) || length(quantile.start) != length(x)) {
    stop("Illegal input for quantile.start")
  }
  if (!is.numeric(quantile.end) || length(quantile.end) != length(x)) {
    stop("Illegal input for quantile.end")
  }

  if (quantile.offset == 'range') {
    quantile.start <- start(x) + quantile.start - 1L
    quantile.end <- start(x) + quantile.end - 1L
  }

  return(IRanges(quantile.start, quantile.end))
}


##' Flags internal priming events using start/end positions found
##' in \code{values(x)[,c(quantile.start, quantile.end)]}. These
##' are defined from a call to \code{boundEventsWithSignalByQuantile}
##'
##' This function runs one-chromosome-at-a-time.
##'
##' @param x A \code{GRanges}-like object
##' @param dna A \code{DNAString} with the genomic sequence for the
##' given chromosome
##' @param a.count The number of A's in the sliding window to flag as primed
##' @param ip.window The size of the sliding window
##' @param ip.distal.downstream How far to go downstream from the annotated
##' event in x to look for priming in "singlar" peak events, or peaks that are
##' the most distal in a series of consecutive peaks of a transcribed event.
##' @param ip.internal.downstream How far to go downstream for peaks that are
##' internal in transcribed events.
##' @param flag.by If any: return TRUE/FALSE vector indicating flag, if
##' 'position', returns the 5'-most position in the peak where the flag was TRUE
##' @param quantile.start The column name to fetch from \code{values(x)}, or
##' an integer vector defining the new starts for events in \code{x}
##' @param quantile.end Same as \code{quantile.start}, but for the ends.
flagInternallyPrimedEventsByBoundedQuantiles <-
function(x, dna, a.window.count=5, abs.a.count=6, ip.window=6,
         ip.distal.downstream=12, ip.internal.downstream=12,
         flag.by=c('position', 'any'), quantile.start='quantile.start',
         quantile.end='quantile.end', min.upstream=10L, ...) {
  x <- as(x, 'GRanges')
  flag.by <- match.arg(flag.by)
  seqname <- unique(as.character(seqnames(x)))
  stopifnot(length(seqname) == 1L)
  if (inherits(dna, 'BSgenome')) {
    dna <- unmasked(dna[[seqname]])
  }
  stopifnot(inherits(dna, "DNAString"))

  if (!is.null(quantile.start) && !is.null(quantile.end)) {
    qbounds <- .reconstructQuantileBounds(x, quantile.start, quantile.end,
                                          ...)
  } else {
    qbounds <- ranges(x)
  }

  ## make sure we scan at least `min.upstream` from end of event
  if (is.numeric(min.upstream)) {
    is.fwd <- as.logical(strand(x) == '+')
    qends <- end(qbounds)
    qends[is.fwd] <- pmin(end(x[is.fwd]) - min.upstream + 1L, qends[is.fwd])
    end(qbounds[is.fwd]) <- pmax(start(qbounds[is.fwd]), qends[is.fwd])

    is.rev <- !is.fwd
    qstarts <- start(qbounds)
    qstarts[is.rev] <- pmax(start(x[is.rev]) + min.upstream - 1L, qstarts[is.rev])
    start(qbounds[is.rev]) <- pmin(end(qbounds[is.rev]), qstarts[is.rev])
  }

  if (flag.by == 'position') {
    is.primed.window <- integer(length(qbounds))
  } else {
    is.primed.window <- logical(length(qbounds))
  }
  is.primed.abs <- is.primed.window

  for (.strand in c('+', '-')) {
    idx <- which(as.logical(strand(x) == .strand))
    meta <- values(x[idx])
    txbounds <- ranges(x[idx])
    sbounds <- qbounds[idx]
    if (length(sbounds) > 0) {
      if (.strand == '-') {
        primed.by <- 'T'
        take.pos <- min
        end(sbounds) <- start(sbounds)
        start(sbounds) <- start(txbounds) - ip.window - ip.distal.downstream
        ## start(sbounds[!meta$is.distal]) <- start(sbounds[!meta$is.distal]) -
        ##   (ip.internal.downstream - ip.distal.downstream)
      } else {
        primed.by <- 'A'
        take.pos <- max
        start(sbounds) <- end(sbounds)
        end(sbounds) <- end(txbounds) + ip.window + ip.distal.downstream
        ## end(sbounds[!meta$is.distal]) <- end(sbounds[!meta$is.distal]) +
        ##   (ip.internal.downstream - ip.distal.downstream)
      }

      views <- Views(dna, sbounds)
      lfreq <- letterFrequencyInSlidingView(views, ip.window, primed.by)
      abs.pattern <- paste(rep(primed.by, abs.a.count), collapse="")
      abs.pattern <- sprintf('%s+', abs.pattern)
      abs.matches <- gregexpr(abs.pattern, as(views, 'character'))

      if (flag.by == 'any') {
        ## Sliding window
        ipw <- sapply(lfreq, function(freq) {
          max(freq) >= a.window.count
        })

        ## absolute count
        ipa <- sapply(abs.matches, function(xx) xx[1] != -1L)
      } else {
        ## sliding window
        axe.pos <- sapply(lfreq, function(freq) {
          ## return the index of the most 5' position >= a.count
          ## returns 0 if threshold is not broken
          idxs <- which(freq >= a.window.count)
          if (length(idxs) > 0) take.pos(idxs) else 0L
        })

        ipw <- start(sbounds) + axe.pos - 1L
        if (.strand == '-') {
          ipw <- ipw + ip.window
        }
        ipw[axe.pos == 0L] <- 0L

        ## absolute count
        if (.strand == '+') {
          pos <- sapply(abs.matches, function(xx) {
            xx[length(xx)] - 1L
          })
        } else {
          pos <- sapply(abs.matches, function(xx) {
            if (xx[1L] == -1L) {
              return(-1L)
            }
            xx[1L] + attr(xx, 'match.length')[1L] - 1L
          })
        }
        ipa <- start(sbounds) + pos
        ipa[pos < 0] <- 0L
      }
      is.primed.window[idx] <- ipw
      is.primed.abs[idx] <- ipa
    }
  }

  data.frame(window=is.primed.window, abs=is.primed.abs)
}

## Slide a window *inside* the event and see how many windows are A-rich
## kill these. This function ignores the quantile stuff and scans from
## the peak position to the end of the event, using more permissive filter.
flagInternalInternallyPrimedEventsByBoundedQuantiles <-
function(x, dna, internal.a.window.count=4, internal.ip.window=5,
         internal.threshold=.6, internal.flag.by=c('position', 'score', 'any'),
         quantile.start='quantile.start', quantile.end='quantile.end',
         return.val=c('score', 'call'), ...) {
  x <- as(x, 'GRanges')
  internal.flag.by <- match.arg(internal.flag.by)
  seqname <- unique(as.character(seqnames(x)))
  stopifnot(length(seqname) == 1L)
  if (inherits(dna, 'BSgenome')) {
    dna <- unmasked(dna[[seqname]])
  }
  stopifnot(inherits(dna, "DNAString"))
  return.val <- match.arg(return.val)

  if (!is.null(quantile.start) && !is.null(quantile.end)) {
    qbounds <- .reconstructQuantileBounds(x, quantile.start, quantile.end,
                                          ...)
  } else {
    qbounds <- ranges(x)
  }
  ## ## Set the 5' bound to be the peak position
  ## pbounds.fwd <- .reconstructQuantileBounds(x, 'peak.pos', quantile.start)
  ## pbounds.rev <- .reconstructQuantileBounds

  qbounds <- ranges(x)

  ## No peak.pos with edge-peak caller
  ## peak.pos <- values(x)$peak.pos
  ## start(qbounds) <- ifelse(as.logical(strand(x) == '+'),
  ##                          peak.pos - floor((peak.pos - start(qbounds)) / 2),
  ##                          start(qbounds))

  ## end(qbounds) <- ifelse(as.logical(strand(x) == '-'),
  ##                        peak.pos + floor((end(qbounds) - peak.pos) / 2),
  ##                        end(qbounds))

  start(qbounds) <- ifelse(as.logical(strand(x) == '+'),
                           start(qbounds) + floor(0.25 * width(qbounds)),
                           start(qbounds))
  end(qbounds) <- ifelse(as.logical(strand(x) == '-'),
                         end(qbounds) - floor(0.25 * width(qbounds)),
                         end(qbounds))

  if (internal.flag.by == 'position') {
    internal <- integer(length(qbounds))
  } else if (internal.flag.by == 'score') {
    internal <- numeric(length(qbounds))
  } else {
    internal <- logical(length(qbounds))
  }

  for (.strand in c('+', '-')) {
    idx <- which(as.logical(strand(x) == .strand))
    meta <- values(x[idx])
    txbounds <- ranges(x[idx])
    sbounds <- qbounds[idx]

    if (length(sbounds) > 0) {
      if (.strand == '-') {
        primed.by <- 'T'
        take.pos <- min
      } else {
        primed.by <- 'A'
        take.pos <- max
      }

      views <- Views(dna, sbounds)
      lfreq <- letterFrequencyInSlidingView(views, internal.ip.window, primed.by)

      if (internal.flag.by == 'any') {
        ## Sliding window
        fail <- sapply(lfreq, function(lf) {
          (sum(lf >= internal.a.window.count) / length(lf)) >= internal.threshold
        })
      } else if (internal.flag.by == 'position') {
        ## sliding window
        fail <- sapply(lfreq, function(lf) {
          if (is.null(lf)) {
            return(0)
          }
          wfail <- which(lf >= internal.a.window.count)
          if (length(wfail) / length(lf) >= internal.threshold) {
            take.pos(wfail)
          } else {
            0L
          }
        })
      } else {
        fail <- sapply(lfreq, function(lf) {
          sum(lf >= internal.a.window.count) / length(lf)
        })
      }

      internal[idx] <- fail

    }
  }

  data.frame(primed=internal)
}


################################################################################
## Old
##' This function works on GRanges that are produced from
##' as(`summarizeSignal()`, 'GRanges')
##' which is to say that the ranges need to have a values() column
##' named peak.pos, which identifies where the SINGULAR peak is inside
##' each range
##'
##' @param x A GRanges object with annotated peaks
##' @param bsg The genome to use to get sequences from
##' @param a.count The number of A's in the sliding window to flag as primed
##' @param ip.window The size of the sliding window
##' @param ip.distal.downstream How far to go downstream from the annotated
##' event in x to look for priming in "singlar" peak events, or peaks that are
##' the most distal in a series of consecutive peaks of a transcribed event.
##' @param ip.internal.downstream How far to go downstream for peaks that are
##' internal in transcribed events.
##' @param flag.by If any: return TRUE/FALSE vector indicating flag, if
##' 'position', returns the 5'-most position in the peak where the flag was TRUE
flagInternallyPrimedPeakAnnotations <- function(x, bsg, a.count=8, ip.window=10,
                                                ip.distal.downstream=0,
                                                ip.internal.downstream=5,
                                                flag.by=c('position', 'any'),
                                                .parallel=FALSE, ...) {
  if (!inherits(x, 'GRanges') && !"peak.pos" %in% colnames(values(x))) {
    stop("Invalid object x")
  }

  if (.parallel && length(chrs) > 1L) {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }

  if (!'is.distal' %in% colnames(values(x))) {
    ## Take a guess at what's "internal" peaks
    ## if (length(x) == 1L) {
    ##   values(x)$is.distal <- TRUE
    ## } else {
    ##   values(x)$is.distal <- end(x[-length(x)]) - start(x[-1L]) < 5
    ## }
    if (ip.distal.downstream != ip.internal.downstream) {
      message("distal peaks haven't been annotated -- did you not use ",
              "getAnnotatedSignal(..., do.split=TRUE)?")
    }
    values(x)$is.distal <- TRUE
  }
  flag.by <- match.arg(flag.by)
  chrs <- unique(as.character(seqnames(x)))
  result <- foreach(chr=chrs, .packages=c('BSgenome')) %loop% {
    cat("===", chr, "===\n")
    bchr <- unmasked(bsg[[chr]])
    sr <- lapply(c('+', '-'), function(.strand) {
      idx <- which(seqnames(x) == chr & strand(x) == .strand)
      xx <- x[idx]
      primed.by <- if (.strand == '+') 'A' else 'T'
      ip.downstream <- ifelse(values(xx)$is.distal, ip.distal.downstream,
                              ip.internal.downstream)
      if (.strand == '+') {
        ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## TODO: start from 1/2 from peak to end to look for priming
        ##       a better idea is to go to the position where x% (25%?) of
        ##       reads in the total peak are
        ## look.start <- values(xx)$peak.pos
        peak.half2end <- values(xx)$peak.pos + ((end(xx) - values(xx)$peak.pos) / 2)
        peak.half2end <- as.integer(floor(peak.half2end))
        look.start <- pmin(peak.half2end, end(xx) - 10L)
        look.end <- end(xx) + ip.window + ip.downstream - 1L
        take.pos <- min
      } else {
        ## Is look.end right for (-) strand?
        look.start <- start(xx) - ip.window - ip.downstream + 1L
        ## look.end <- values(xx)$peak.pos
        peak.half2start <- values(xx)$peak.pos - ((values(xx)$peak.pos - start(xx)) / 2)
        peak.half2start <- as.integer(ceiling(peak.half2start))
        look.end <- pmax(peak.half2start, start(xx) + 10L)
        take.pos <- max
      }
      look.it <- IRanges(look.start, look.end)
      views <- Views(bchr, look.it)
      lfreq <- letterFrequencyInSlidingView(views, ip.window, primed.by)
      if (flag.by == 'any') {
        to.axe <- sapply(lfreq, function(freq) max(freq) >= a.count)
      } else {
        to.axe <- sapply(lfreq, function(freq) {
          ## return the index of the most 5' position >= a.count
          ## returns 0 if threshold is not broken
          idxs <- which(freq >= a.count)
          if (length(idxs) == 0) 0L else take.pos(idxs)
        })
        shifted.axe <- look.start + to.axe - 1L
        if (.strand == '-') {
          shifted.axe <- shifted.axe + ip.window - 1L
        }
        shifted.axe[to.axe == 0L] <- 0L
      }
      list(idx=idx, axe=shifted.axe)
    })
    unlist(sr, recursive=FALSE)
  }
  result <- unlist(result, recursive=FALSE)

  if (flag.by == 'any') {
    axe <- logical(length(x))
  } else {
    axe <- integer(length(x))
  }

  for (idx in seq(1, length(result), by=2)) {
    idxs <- result[[idx]]
    to.axe <- result[[idx + 1L]]
    axe[idxs] <- to.axe
  }

  axe
}

if (FALSE) {
  ## Figure out wtf is going on with read patterns
  source('mayr/setup.R')
  expts <- loadMayr(nickd)
  expt <- expts$nick.mcf10a
  r <- getReadsOnChromosome(expt, 'chr1', smooth.by=NULL)
  c1 <- list(fwd=coverage(r[strand(r) == '+'])[[1]],
             rev=coverage(r[strand(r) == '-'])[[1]])
  sig.bounds <- lapply(c1, slice, iTPM(3, expt))
  shifted.widths <- lapply(sig.bounds, function(x) {
    data.frame(start=start(x), end=end(x), width=width(x), count=viewMaxs(x))
  })
  shifted.widths[[1]]$strand <- '+'
  shifted.widths[[2]]$strand <- '-'
  shifted.widths <- do.call(rbind, shifted.widths)
  rownames(shifted.widths) <- NULL

  trimmed <- subset(shifted.widths, count < quantile(count, .95) &
                    count > quantile(count, .05))
  trimmed.avg <- with(trimmed, median(rep(width, count)))
  trimmed.sd <- with(trimmed, sd(rep(width, count)))

  c1.smooth <- lapply(c1, function(x) Rle(kdeSmooth(x, bandwidth=8, kernel='triweight')))
  ##############################################################################
  ## For ranges of "normal" width
  some <- subset(trimmed, width > (trimmed.avg - 5) & width < (trimmed.avg + 5))
  some <- some[order(some$count, decreasing=TRUE),]

  coverages <- lapply(1:10, function(idx) {
    info <- some[idx,]
    test.reads <- subsetByOverlaps(r, GRanges('chr1', IRanges(info$start, info$end), strand='*'))
    strand.dist <- table(strand(test.reads))
    print(strand.dist)
    test.reads <- ranges(test.reads[strand(test.reads) == info$strand])
    test.reads <- shift(test.reads, -(start(test.reads[1]) - 1))
    this.cov <- as.integer(coverage(test.reads))
    this.cov
  })

  plot.signals <- function(signal.list) {
    cols <- rainbow(length(signal.list))
    xmax <- max(sapply(signal.list, length))
    for (idx in 1:10) {
      this.cov <- signal.list[[idx]]
      if (idx == 1) {
        plot(this.cov / max(this.cov), col=cols[idx], xlim=c(0, xmax), type='l')
      } else {
        lines(this.cov / max(this.cov), col=cols[idx])
      }
    }
  }

  cov.tri.10 <- lapply(coverages, kdeSmooth, bandwidth=10, kernel='triweight')
  quartz()
  plot.signals(cov.tri.10)

  cov.tri.5 <- lapply(coverages, kdeSmooth, bandwidth=5, kernel='triweight')
  quartz()
  plot.signals(cov.tri.5)

  cov.padd <- lapply(coverages, function(x) c(0, 0, 0, 0, 0, x, 0, 0, 0, 0, 0))
  plot.signals(cov.padd)
  quartz()
  plot.signals(lapply(cov.padd, kdeSmooth, bandwidth=10, kernel='normal'))
  ##############################################################################
  ## For "wider" ranges
  ## Use reads over S100A6 as a test case -- this gene is super expressed
  s10.reads <- subsetByOverlaps(r, GRanges('chr1', IRanges(153507051, 153507230), strand="*"))
  neg.cov <- as.integer(coverage(s20.reads[strand(s10.reads) == '-'])$chr1[153507051:153507230])
  pos.cov <- as.integer(coverage(s10.reads[strand(s10.reads) == '+'])$chr1[153507051:153507230])

  plot(as.integer(neg.cov) / max(neg.cov), type='l')

  ## beta
  rb <- rbeta(1000, 2, 10)
  rbi <- round(rb * 100)
  rbi.ranges <- IRanges(rbi, width=36)
  cov.rbi <- as.integer(coverage(rbi.ranges))
  cov.rbi <- cov.rbi / max(cov.rbi)
  lines(cov.rbi, col='red')
  points(as.integer(pos.cov), col='red')
  points(as.integer(neg.cov) + as.integer(pos.cov), col='green')

  ## poisson?
  rp <- rpois(1000, 10)
  ir.pois <- IRanges(rp, width=36)
  cov.pois <- as.integer(coverage(ir.pois))
  conv.pois <- convolve(neg.cov, cov.pois, type='open')
  plot(conv.pois, xlim=c(0, 180))

  cov.rbi <- (cov.rbi / max(cov.rbi)) / 2
  conv.rbi <- convolve(neg.cov, rev(cov.rbi), type='open')
  plot(conv.rbi, xlim=c(0, 180))
  ## smooth with gausian?
  smoothing.gauss <- rnorm(1000, 10)
  ig <- as.integer(smoothing.gauss)
  points(smoothing.gauss + 50, col='blue')
}

