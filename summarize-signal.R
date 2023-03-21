## WIN: For doing hump finding (in SAPAS/mcf10a experiment).
##
##   reads <- getReadsOnChromosome(sapas$mcf10a.sapas, 'chr1',
##                                 unique.only=FALSE, max.mismatch=NULL,
##                                 smooth.by=NULL, .dataset='realign.only')
##   region <- GRanges('chr1', IRanges(25170692, 25170822))
##   visualizeAntisense(reads, region, threshold=8, kernel='normal')
##
##   Why:
##   The middle region is bogus, but forms an artificial "antenna" where the
##   end of first event and start of third overlap. This overlap region is
##   correctly assigned a (read) count of 0.
##
##   seqnames             ranges strand  | count  covr event.start event.end
##   chr1   [25170692, 25170731]      +  |   333 458.7    25170692  25170817
##   chr1   [25170732, 25170777]      +  |     0 598.0    25170692  25170817
##   chr1   [25170778, 25170817]      +  |   286 394.0    25170692  25170817


## TODO: Move signals into an internal environment, so we can keep signals
##       over multiple chromosomes
setClass("AnnotatedSignal", contains="GRanges",
         representation=representation(
           signal.fwd='Rle',
           signal.rev='Rle'),
         prototype=prototype(
           signal.fwd=Rle(rep(0, 0)),
           signal.rev=Rle(rep(0, 0))))

setMethod("show", "AnnotatedSignal", function(object) {
  output <- sprintf("AnnotatedSignal over %s with %d islands.\n",
                    as.character(object@seqnames[1]),
                    length(object@ranges))
  cat(output)
})

setMethod("c", c(x="AnnotatedSignal"),
function(x, ..., .ignore.signals=TRUE, recursive=TRUE) {
  if (!.ignore.signals) {
    warning("We are ignoring signals -- you have no choice.")
    .ignore.signals <- TRUE
  }
  cc <- selectMethod("c", "GenomicRanges")
  args <- unname(list(x, ...))
  if (.ignore.signals) {
    args <- lapply(args, function(x) {
      x@signal.fwd <- Rle(rep(0, 0))
      x@signal.rev <- Rle(rep(0, 0))
      x
    })
  } else {
    ## Implement later
  }
  suppressWarnings(do.call(cc, args))
})

##' Drops the signal and converts it to more like a CompressedReads
##' object. The transcribed events are broken into their maxima-pieces
##'
##' A 'is.distal' column is added to the values() of the returned object
##' indicating whether or not (TRUE/FALSE) the peak is the most distal
##' one in its event. If the peak is solitary in the event, then it is
##' defacto TRUE, otherwise it is FALSE if there is a peak more downstream
##' from this one in the same transcribed event.
setAs("AnnotatedSignal", "GRanges", function(from) {
  chrs <- as.character(unique(seqnames(from)))

  events <- lapply(chrs, function(chr) {
    x <- from[seqnames(from) == chr]
    xmeta <- values(x)
    peaks <- unlist(xmeta$peak.info)
    meta <- values(peaks)
    parts <- xmeta$peak.info@partitioning

    ## DEBUG: Some events after the following shifts are off by one
    ## (on [-] strand ony?)
    ps <- shift(peaks, meta$shift)

    ## remap the meta$event.idx to point to the appropriate range in the
    ## x GRanges object.

    ## the signal was built by first looking at coverage on (+) strand, then
    ## doing the (-) after, so all regions are in ascending order on (+),
    ## then resets to "start" on (-). The switch.idx is the position of the
    ## last (+) range in the unwound/unlisted `peaks`
    switch.after <- which(diff(meta$event.idx) < 0)[1L]
    single.stranded <- length(x) == 1 || is.na(switch.after)
    if (single.stranded) {
      strands <- rep(strand(x)[1L], nrow(meta))
      these.idx <- meta$event.idx
    } else {
      pos.idx <- meta$event.idx[1:switch.after]
      neg.idx <- meta$event.idx[-(1:switch.after)]
      neg.idx <- neg.idx + which(strand(x) == '-')[1L] - 1L
      these.idx <- c(pos.idx, neg.idx)
      strands <- c(rep('+', length(pos.idx)), rep('-', length(neg.idx)))
    }

    iranges <- IRanges(start(ps), end(ps))
    event.start <- start(x)[these.idx]
    event.end <- end(x)[these.idx]

    ## figure out 'distal' ranges
    .x <- data.table(strand=strands, group=togroup(parts), key='group')
    is.distal <- .x[, {
      len <- length(group)
      if (len == 1L) {
        distal <- TRUE
      } else{
        distal <- logical(len)
        if (strand[1] == '+') {
          distal[length(distal)] <- TRUE
        } else {
          distal[1L] <- TRUE
        }
      }
      distal
    }, by='group']$V1
    gr <- GRanges(chr, iranges, strands,
                  count=meta$event.read.count, covr=meta$peak.covr,
                  event.start=event.start, event.end=event.end,
                  is.distal=is.distal)
    values(gr)$peak.pos <- values(gr)$event.start + meta$peak.pos -1L
    gr
  })

  suppressWarnings(do.call(c, events))
})

## Ensure that processors don't sit around waiting when wrapping up
## one experiment and going on to another
doSignalFirehose <- function(expts, chrs, out.dir=NULL, ...) {
  mcopts <- list(preschedule=FALSE)
  foreach(expt=expts, .options.multicore=mcopts, .inorder=FALSE) %:%
    foreach(chr=chrs, .options.multicore=mcopts, .inorder=FALSE) %dopar% {
      cat("===", experimentName(expt), ":", chr, "===\n")
      if (!is.null(out.dir)) {
        save.path <- file.path(metaDirectory(expt, ...), 'smooth.signal')
      } else {
        save.path <- NULL
      }
      summarizeSignal(expt, chrs=chr, do.save=TRUE, save.path=save.path,
                      .parallel=FALSE, .trace=FALSE, ...)
  }
  invisible(NULL)
}

##' Processes the read events, converts them into smoothed signals and annotates
##' signals with peak information.
summarizeSignal <- function(expt, kernel='normal', bandwidth=18,
                            chrs=seqnames(expt), unique.only=FALSE,
                            max.mismatch=NULL, with.readcount=FALSE,
                            cached.reads=NULL, do.save=TRUE,
                            save.path=NULL, .parallel=TRUE, .trace=TRUE, ...) {
  if (.parallel && length(chrs) > 1L) {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }

  if (do.save && is.null(save.path)) {
    save.path <- file.path(metaDirectory(expt, ...), 'annotated.signal')
  }
  result <- foreach(chr=chrs, .packages=c('GenomicRanges', 'IRanges'),
                    .options.multicore=list(preschedule=FALSE),
                    .inorder=FALSE) %loop% {
    if (.trace) {
      cat("===", chr, "===\n")
    }
    if (is.null(cached.reads)) {
      reads <- getReadsOnChromosome(expt, chr, unique.only=unique.only,
                                    max.mismatch=max.mismatch, smooth.by=NULL,
                                    ...)
    } else {
      reads <- cached.reads[seqnames(cached.reads) == chr]
    }

    if (length(reads) == 0L) {
      return(NULL)
    }

    ## values(reads) <- values(reads)[, c('tag.X0', 'tag.X1', 'tag.NM')]
    values(reads) <- NULL

    info <- lapply(c('+', '-'), function(.strand) {
      r <- reads[strand(reads) == .strand]
      if (length(r) == 0L) {
        return(list(islands=GRanges(), signal=Rle(rep(0, 0))))
      }

      read.tree <- IntervalTree(ranges(r))
      covr <- coverage(ranges(r))
      islands <- slice(covr, lower=0L, includeLower=FALSE, rangesOnly=TRUE)
      smooth <- kdeSmooth(covr, kernel, bandwidth)
      views <- Views(smooth, islands)

      ## NOTE: Using viewMaxs of the coverage within each event produces read
      ##       counts for "events" that are <= the counts we get from
      ##       compressReads(smoother.repackTags(reads)). In sufficiently wide
      ##       events, the reads won't have coverage (first case) as high as
      ##       the number of reads within the event bounds (second case).
      ##       countOverlaps should produce a number closer to the 2nd case.
      covr.maxs <- viewMaxs(views)
      island.counts <- countOverlaps(islands, ranges(r))

      ## if (FALSE) {
      ##   debug <- which(covr.maxs > 30)
      ##   debug <- sample(debug, 50)
      ##   views <- views[debug]
      ##   covr.maxs <- covr.maxs[debug]
      ##   island.counts <- island.counts[debug]
      ## }

      peak.info <- lapply(seq_len(length(views)), function(idx) {
        bounds <- views[idx]
        signal <- views[[idx]]
        info <- tryCatch(findPeakBoundsInSignal(signal, bounds),
                         error=function(e) NULL)
        if (is.null(info)) {
          return(IRanges())
        }
        events <- shift(info, values(info)$shift)
        values(info)$peak.covr <- viewMaxs(Views(signal, info))
        if (with.readcount) {
          values(info)$event.read.count <- countOverlaps(events, read.tree)
        } else {
          values(info)$event.read.count <- rep(NA_integer_, length(info))
        }
        values(info)$event.idx <- idx
        info
      })

      n.peaks <- sapply(peak.info, length)
      ## Find and remove "error islands"
      is.good <- n.peaks > 0L
      peak.info <- peak.info[is.good]
      views <- views[is.good]
      n.peaks <- n.peaks[is.good]
      covr.maxs <- covr.maxs[is.good]
      event.read.counts <- island.counts[is.good]

      ret.value <- GRanges(chr, as(views, 'IRanges'), strand=.strand,
                           peak.info=IRangesList(peak.info),
                           n.peaks=n.peaks,
                           event.covr=covr.maxs,
                           event.read.counts=island.counts)

      list(islands=ret.value, signal=smooth)
    })
    asignal <- c(info[[1]]$islands, info[[2]]$islands)
    class(asignal) <- 'AnnotatedSignal'
    slot(asignal, 'signal.fwd') <- info[[1]]$signal
    slot(asignal, 'signal.rev') <- info[[2]]$signal

    if (do.save) {
      fn <- signalFN(save.path, chr, kernel, bandwidth, unique.only, max.mismatch, ...)
      cat("... saving", fn, "\n")
      save(asignal, file=fn)
      chr
    } else {
      asignal
    }
  }

  invisible(result)
}

##' Finds boundaries of events as defined by the quantiles given for
##' distribution of reads across the event.
##'
##' This function only works over GRanges (and signals) defined over
##' one chromosome at a time.
##'
##' @param x A \code{GRanges} or \code{AnnotatedSignal} object. If
##' it's a \code{GRanges} object, pass in signals to use through the
##' \code{signal.fwd, signal.ref} parameters
##' @param start.quantile
##' @param end.quantile with respect to the left-most position of event,
##' (ie. strand is ignored)
##' @param covr \code{character} or an \code{integer} vector of
##' \code{length(start.quantile) - 1L} which indicates breaks in covr where
##' each value in the quantiles applies to. You should really stick with using
##' the \code{count} column, which is added after signal post processing, since
##' wide events can have a lower coverage but larger count since they are spread
##' thin. If \code{length(covr) > 1} then the pmax (or pmin) of these vectors
##' is taken for each 'island'
##' @param signal.fwd The signal on the fwd strand of the chromosome.
##' @param signal.rev The signal on the rev strand of the chromosome.
boundEventsWithSignalByQuantile <-
function(x, start.quantile=c(0.2, 0.1, 0.05), end.quantile=1-start.quantile,
         covr.breaks=c(0, 6, 20), covr=c('count', 'covr'), ptake=pmax,
         signal.fwd=NULL, signal.rev=NULL, ...) {
  stopifnot(length(start.quantile) == length(end.quantile))

  if (is.null(signal.fwd) || is.null(signal.rev)) {
    if (!inherits(x, "AnnotatedSignal")) {
      stop("Provide an AnnotatedSignal object, or send signals in explicitly")
    }
    if (is.null(signal.fwd)) {
      signal.fwd <- x@signal.fwd
    }
    if (is.null(signal.rev)) {
      signal.rev <- x@signal.rev
    }
    if (any(!values(x)$is.distal)) {
      warning("The boundedEvents will be longer than the AnnotatedSignal",
              immediate.=TRUE)
    }
  }

  ## Splits the signal into "pieces"
  x <- as(x, 'GRanges')
  stopifnot(is(x, 'GRanges'))
  bounds <- ranges(x)

  if (length(start.quantile) > 1L) {
    if (is.null(covr)) {
      covr <- 'covr'
    }
    if (is.character(covr)) {
      if (!all(covr %in% colnames(values(x))) || length(covr) == 0) {
        stop("Columns used for covr/count calculations not found")
      }
      if (length(covr) > 1) {
        covr <- do.call(ptake, lapply(covr, function(name) values(x)[[name]]))
      } else {
        covr <- values(x)[[covr]]
      }
    }

    if (!is.numeric(covr)) {
      stop("Need a coverage vector here")
    }
    if (length(covr) != length(x)) {
      stop("Coverage vector isn't same length as `x`")
    }
    which.quantile <- as.integer(cut(covr, covr.breaks))
    too.high <- is.na(which.quantile)
    if (any(too.high)) {
      which.quantile[too.high] <- length(covr.breaks)
    }
  } else {
    which.quantile <- rep(1L, length(x))
  }

  is.fwd <- as.logical(strand(x) == '+')

  squant <- start.quantile[which.quantile]
  equant <- end.quantile[which.quantile]
  ## TODO: Try to make {start|end}.quantile strand specific
  ## start.quantile <- ifelse(is.fwd, squant, equant)
  ## end.quantile <- ifelse(is.fwd, equant, squant)
  start.quantile <- squant
  end.quantile <- equant

  fwd <- TRUE
  for (idxs in list(which(is.fwd), which(!is.fwd))) {
    if (fwd) {
      sigs <- Views(signal.fwd, bounds[idxs])
      fwd <- FALSE
    } else {
      sigs <- Views(signal.rev, bounds[idxs])
    }

    inner.idx <- seq_along(idxs)
    boundses <- sapply(inner.idx, function(idx) {
      outter.i <- idxs[idx]
      signal <- sigs[[idx]]
      percent <- cumsum(signal) / sum(signal)
      start.idx <- which(percent >= start.quantile[outter.i])[1L]
      end.idx <- which(percent >= end.quantile[outter.i])[1L]
      c(start.idx, end.idx)
    })
    end(bounds[idxs]) <- start(bounds[idxs]) + boundses[2, ] - 1L
    start(bounds[idxs]) <- start(bounds[idxs]) + boundses[1,] - 1L
  }

  bounds
}

getAnnotatedSignal <- function(expt, seqname, kernel='normal', bandwidth=18,
                               unique.only=FALSE, max.mismatch=NULL,
                               do.split=FALSE) {
  if (missing(seqname)) {
    seqname <- seqnames(expt)
  }
  signals <- lapply(seqname, function(sn) {
    fn <- signalFN(expt, sn, kernel, bandwidth, unique.only, max.mismatch)
    x <- load.it(fn)
    if (do.split) {
      x <- as(x, 'GRanges')
    }
    x
  })
  if (length(signals) > 1) {
    signals <- suppressWarnings(do.call(c, signals))
  } else {
    signals <- signals[[1]]
  }
  signals
}

## The files that this function relies on are generated by
## these were generated from remapping
getProcessedSignal <- function(expt, chr=NULL, do.split=TRUE,
                               ...) {
  if (!do.split) {
    cat("Sorry -- it's split")
  }
  x <- load.it(file.path(metaDirectory(expt, ...), 'annotated.signal',
                         'all.signals.rda'))
  if (!is.null(chr)) {
    if (!chr %in% seqlevels(expt)) {
      return(GRanges())
    }
    x <- x[seqnames(x) == chr]
  }
  if (!'tpm' %in% colnames(values(x))) {
    values(x)$tpm <- TPM(values(x)$count, expt)
  }
  x
}

getSmoothSignal <- function(x, chr, kernel='normal',
                            bandwidth=18L, ...) {
  stopifnot(inherits(x, 'TailSeq'))
  if (missing(chr)) {
    stop("chromosome name and strand required")
  }
  sig.dir <- file.path(metaDirectory(x, ...), 'smooth.signal')
  if (!dir.exists(sig.dir)) {
    stop("No smooth.signal directory")
  }
  file.pattern <- sprintf('.*\\.%s.bandwidth-%d.*%s\\.rda',
                          kernel, bandwidth, chr)
  the.file <- list.files(sig.dir, file.pattern, full.names=TRUE)

  if (length(the.file) != 1L) {
    if (length(the.file) > 1) {
      msg <- "Ambiguous pattern [%s] matches many files in %s"
    } else {
      msg <- "Pattern [%s] does not match many files in %s"
    }
    msg <- sprintf(msg, file.pattern, sig.dir)
    stop(msg)
  }

  load.it(the.file)
}

## The GRanges defined in \code{events} signify the regions that have some
## read coverage. The \code{is.real} vector is as long as \code{events} and
## indicates wether each region is "real" ie. one that isn't internally primed,
## or ones that are but have a polya signal.
.calculateEventCoverage <- function(expt, chr, strand, events, is.real,
                                    unique.only=FALSE, max.mismatch=NULL,
                                    kernel=NULL, bandwidth=18, ...) {
  zero <- Rle(0, seqlengths(expt)[chr])
  take.idx <- as.logical(seqnames(events) == chr & strand(events) == strand)
  is.real <- is.real[take.idx]
  if (sum(is.real) == 0L) {
    return(zero)
  }
  events <- events[take.idx]
  reads <- getReadsOnChromosome(expt, chr, strand=strand,
                                unique.only=unique.only,
                                max.mismatch=max.mismatch, smooth.by=NULL, ...)
  if (length(reads) == 0) {
    return(zero)
  }
  event.idx <- assignUniqueOverlaps(reads, events, ...)
  if (any(is.na(event.idx))) {
    warning("NA events in ", chr, " -- there are reads in regions we didn't ",
            "mark in `events`\n")
  }
  take.real <- is.real[event.idx]
  take.real[is.na(take.real)] <- FALSE
  real.reads <- reads[take.real]
  covr <- coverage(ranges(real.reads))
}

##' Creates signal files for genome-browser-viewing from the original reads.
##'
##' Use this function after internally-primed/artifactual regions are
##' removed and pass the "clean" ranges in through \code{events}.
##'
##' 2011-04-27:
##'
##' An example workflow for a virgin experiment would be::
##'
##'   > registerDoMC(8) ## use ika
##'   > bsg <- getBsGenome(expt)
##'   > annotated.genome <- getAnnotatedGenome(gcr)
##'   > summarizeSignal(expt)
##'   > ## wait for (4) hours
##'   > peaks <- getAnnotatedSignals(expt)
##'   > peaks <- peaks[seqnames(peaks) != 'chrM']
##'   > anno <- annotateReads(peaks, annotated.genome, nuke.antisense=FALSE)
##'   > values(anno)$is.primed <- flagInernallyPrimedPeakAnnotations(anno, bsg)
##'   > real.events <- with(values(anno), hasPolyA | !is.primed)
##'   > serializeSignal(expt, events, 'some.name', 'bigWig')
serializeSignal <- function(expt, events, real.events=values(events)$is.real,
                            unique.only=FALSE, max.mismatch=NULL,
                            kernel=NULL, bandwidth=18, chrs=NULL,
                            prefix=experimentName(expt), save.path='.',
                            format=c('bigWig', 'wig', 'bedGraph'), ...) {
  format <- match.arg(format)
  if (is.null(chrs)) {
    chrs <- as.character(unique(seqnames(events)))
  }

  ## filename handling
  fn.base <- prefix
  if (!is.null(kernel)) {
    fn.base <- paste(fn.base, ".", kernel, "-bw", bandwidth, sep="")
  }
  fn.base <- file.path(save.path, fn.base)

  covr <- foreach(chr=chrs, .packages=c("GenomicRanges")) %dopar% {
    cat("===", chr, "(+) ===\n")
    .calculateEventCoverage(expt, chr, '+', events, real.events,
                            unique.only=unique.only, max.mismatch=max.mismatch,
                            kernel=kernel, bandwidth=bandwidth)
  }
  names(covr) <- chrs

  ## serialize bigWig
  rd <- as(RleList(covr), 'RangedData')
  genome(rd) <- genome(expt)
  fn <- paste(fn.base, 'fwd', format, sep=".")
  cat("saving to:", fn, "\n")
  export.bw(rd, fn)
  rm(covr, rd)
  gc()

  covr <- foreach(chr=chrs, .packages=c("GenomicRanges")) %dopar% {
    cat("===", chr, "(-) ===\n")
    .calculateEventCoverage(expt, chr, '-', events, real.events,
                            unique.only=unique.only, max.mismatch=max.mismatch,
                            kernel=kernel, bandwidth=bandwidth)
  }
  names(covr) <- chrs

  ## serialize bigWig
  rd <- as(RleList(covr), 'RangedData')
  genome(rd) <- genome(expt)
  fn <- paste(fn.base, 'rev', format, sep=".")
  cat("saving to:", fn, "\n")
  export.bw(rd, fn)
  rm(rd, covr)
  gc()

  invisible(NULL)
}


eventDiff <- function(granges.1, granges.2) {
  o <- findOverlaps(granges.1, granges.2)
  take <- setdiff(1:length(granges.1), unique(queryHits(o)))
  granges.1[take]
}


