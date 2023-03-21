###############################################################################
## This is defunkt -- using biosignals package now
###############################################################################

##' Finds peaks/troughs across an entire coverage vector at once.
##'
##' This function only works on one coverage vector at a time, so each strand
##' should be done separately if your coverages are strand specific.
##'
##' @param cvr.smooth A \emph{smoothed} coverage vector. It is important that it
##' is smoothed -- this function won't check that for you.
##' @param which.strand '+' or '-' indicating which strand this signal is from
##' @param eps The "tolerance". Differnces between two numbers that are less
##' than abs(eps) are considered to be 0.
##' @param lower The number of reads required to define an event boundary. Make
##' this number proportional to the number of experiments.
##' @param min.width The minimum width of events -- using a number that is the
##' same size as the minimum length of reads in the experiment isn't a bad idea.
specializedSignalHumps <- function(cvr.smooth, which.strand, eps=1e-6, lower=5,
                                   min.width=21L, quantile.breaks=NULL,
                                   with.quantiles=!is.null(quantile.breaks)) {
  if (missing(which.strand)) {
    stop("strand information required.")
  }

  bounds <- slice(cvr.smooth, lower, includeLower=lower > 0, rangesOnly=TRUE)
  bounds <- bounds[width(bounds) >= min.width]

  ## Cut "long tails" off of events to avoid finding peaks in them -- most
  ## are due to noise. See the "bumpy-signal-need-quantiles" ACE thing.
  quants <- defineQuantilePositions(cvr.smooth, bounds, quantile.breaks)
  qbounds <- .reconstructQuantileBounds(bounds, quants$quantile.02,
                                        quants$quantile.98)

  ## Fast peak calling requires converting the Rle to a numeric vector.
  monster <- as.numeric(cvr.smooth)
  monster.info <- specialized.signal.humps.rcpp(monster, eps)

  ## Associate the min/max fenceposts to regions defined in `bounds`
  mins <- monster.info$minima
  ## If the coverage vector ended at the last non-zero position, then there
  ## is no "final" minimum after that event, so we add that here.
  if (as.logical(cvr.smooth[length(cvr.smooth)] > 0)) {
    mins <- c(mins, length(cvr.smooth))
  }

  maxs <- monster.info$maxima
  ## There are many spurios minima, so we filter them out by only taking
  ## the ones that fencepost the peaks
  mins.left <- findInterval(maxs, mins)
  peaks <- IRanges(maxs, width=1L)
  values(peaks) <- DataFrame(left=mins[mins.left], right=mins[mins.left + 1L])

  ## Score peaks based on how much they clear their proximal minima
  is.noise <- peakIsNoise(start(peaks), values(peaks)$left,
                          values(peaks)$right, monster)
  values(peaks)$is.noise <- is.noise
  ## Removal of noise here is new (2011-10-17)
  peaks <- peaks[!is.noise]

  ## sanity check
  is.bad <- values(peaks)$left >= start(peaks) | start(peaks) >= values(peaks)$right

  if (any(is.bad)) {
    stop(sum(is.bad), "bad boundary calls")
  }

  ##############################################################################
  ## Only take peaks that are within the 'trimmed' events
  o <- findOverlaps(peaks, qbounds)
  pidx <- queryHits(o)
  bidx <- subjectHits(o)

  ## (1) Figure out which peaks belong to "multi-modal" events
  dt <- data.table(distal=FALSE, bidx=bidx, key='bidx')
  if (which.strand == '-') {
    is.distal <- dt[, {
      .distal <- distal
      .distal[1L] <- TRUE
      list(distal=.distal)
    }, by=bidx]$distal
  } else {
    is.distal <- dt[, {
      .distal <- distal
      .distal[length(.distal)] <- TRUE
      list(distal=.distal)
    }, by=bidx]$distal
  }

  ## (2) Tighten the boudnaries of the events.
  ## The Rcpp peak caller make events wider under some sceanrios:
  ##
  ##  (i) primarily fix the tails that were cut off by the slice() operation
  ##      where there values are non-zero, but less than the parameter set in
  ##      `lower`.
  ##
  ## (ii) prior to doing the `qbounds` mojo above, we also had to accomdate
  ##      noise in wide/shallow tails
  ##
  ## Boundaries either the ones found by the peak finder, or the ones specified
  ## as the start()/end() coords in `bounds`, whichever is "tigther"
  meta <- values(peaks[pidx])
  p.starts <- meta$left
  p.ends <- meta$right
  is.noise <- meta$is.noise
  pp <- start(peaks[pidx])

  b.starts <- start(bounds)[bidx]
  b.ends <- end(bounds)[bidx]

  ## Tighten the boundaries of the events --
  re.start <- ifelse(p.starts < b.starts, b.starts, p.starts)
  re.end <- ifelse(p.ends > b.ends, b.ends, p.ends)

  rebounds <- IRanges(re.start, re.end)
  ## Check noise calls
  if (FALSE) {
    ## TODO: You'll see that sometimes peaks that are noise should be
    ## merged with a proximal peak.
    redf <- data.frame(left=start(rebounds), right=end(rebounds),
                       width=end(rebounds) - start(rebounds) + 1L, peak=pp,
                       noise=is.noise,
                       cvr.left=monster[start(rebounds)], cvr.peak=monster[pp],
                       cvr.right=monster[end(rebounds)], distal=is.distal)
    redf.n <- subset(redf, noise)
  }

  is.bad <- pp < start(rebounds) | pp > end(rebounds)
  if (any(is.bad)) {
    stop(sum(is.bad), "rebound boundaries")
  }

  ## The start/end of consecutive peaks in a "singular" transcribed event
  ## are the same, let's separate them.
  if (which.strand == '-') {
    start(rebounds) <- start(rebounds) + ifelse(is.distal, 0, 1L)
  } else {
    end(rebounds) <- end(rebounds) - ifelse(is.distal, 0, 1L)
  }

  ## (3) Package up for shipping out
  values(rebounds) <- DataFrame(peak.pos=pp, is.distal=is.distal)

  if (with.quantiles) {
    qdf <- defineQuantilePositions(monster, rebounds, quantile.breaks)
    values(rebounds) <- cbind(values(rebounds), qdf)
  }

  ## Remove 'noise' if it's not the distal peak
  rebounds <- rebounds[!is.noise | values(rebounds)$is.distal]

  ## If there are any overlapping regions at this point, they should be merged
  mm <- as.matrix(findOverlaps(rebounds, ignoreSelf=TRUE, ignoreRedundant=TRUE))
  if (nrow(mm) > 0) {
    cat("... some overlapping ranges left -- merging\n")
    rebounds <- mergeOverlaps(rebounds, mm)
    if (length(findOverlaps(rebounds, ignoreSelf=TRUE)) > 0) {
      stop("You need to do a recursive merge of overlapping peaks!")
    }
  }
  rebounds
}

## This is not recursive -- assumes only 1 overlap to merge
mergeOverlaps <- function(x, match.matrix=NULL) {
  if (is.null(match.matrix)) {
    match.matrix <- findOverlaps(x, ignoreSelf=TRUE, ignoreRedundant=TRUE)
    match.matrix <- as.matrix(match.matrix)
  }
  stopifnot(is.matrix(match.matrix))
  if (nrow(match.matrix) == 0) {
    cat("... no merging necessary\n")
    return(x)
  }

  new.ones <- lapply(1:nrow(match.matrix), function(i) {
    i1 <- match.matrix[i,1]
    i2 <- match.matrix[i,2]
    template <- x[c(i1, i2)]
    meta <- values(template)
    meta$peak.pos <- as.integer(mean(meta$peak.pos))
    meta$is.distal <- any(meta$peak.pos)
    nranges <- range(template)
    values(nranges) <- meta[1,]
    nranges
  })
  new.ones <- do.call(c, unname(new.ones))
  x <- x[-unique(as.integer(match.matrix))]
  x <- c(x, new.ones)
  x[order(x)]
}

peakIsNoise <- function(peak.pos, left, right, values, percent.clear=0.05,
                        min.count=10, max.count=50) {
  peak.cvr <- as.numeric(values[peak.pos])
  min.left <- as.numeric(values[left])
  min.right <- as.numeric(values[right])

  pclear.count <- pmin(peak.cvr - min.left, peak.cvr - min.right)
  pclear.percent <- pclear.count / peak.cvr

  pclear.count < min.count | pclear.percent < percent.clear
}

## This is in biosignals now
if (FALSE && !exists('specialized.signal.humps.rcpp')) {
  ## Specialized Rcpp version to define "event characteristics".
  ##
  ## This assumes that x is one "signal/transcribed-event", so all values must
  ## be non-negative.
  ##
  ##   (i) A minima will be added at the last zero element (or first nonzero)
  ##       element on the left side to "start" the event.
  ##  (ii) A minima will be added at the end of `x` if there is a maxima
  ##       reported without a trailing minima
  ##  (ii) A "flat" event will have a peak called at its center.
  ##
  ## All indices reported by this function are 1-based.
  specialized.signal.humps.rcpp <- cxxfunction(signature(x="numeric", eps="numeric"), '
  Rcpp::NumericVector xx(x);
  Rcpp::NumericVector epsilon(eps);
  double prev_slope = 0;
  double prev_slope_at_lag = 0;
  double lag_start = -1;
  double prev_value = 0;
  double slope = 0;
  double value = 0;
  double prevmax = -1;
  double precision = epsilon[0];
  int nxx = xx.size();
  int i;
  bool check_for_lost_maxima;
  int last_zero_pos = -1;
  int has_signal;
  int previous_turn = 0; // 1 = peak, -1 = minima
  int is_climbing = 0;

  std::vector<int> maxima;
  std::vector<int> minima;

  for (i = 0; i < nxx; i++) {
    check_for_lost_maxima = false;

    value = xx[i];
    if (value <= precision) {
      value = 0;
    }

    slope = value - prev_value;
    if ((-precision <= slope) && (slope <= precision)) {
      slope = 0;
    }

    //Rprintf("idx:%d, value:%.3f, slope:%.3f, prev.slope:%.3f\\n", i, value, slope, prev_slope);

    if (value == 0) {
      last_zero_pos = i;
      if (prev_value > 0) {
        // We just "finished" a transcribed region, add a minimum at this index.
        //Rprintf("  Finished an event\\n");
        minima.push_back(i + 1); // C -> R indexing
        check_for_lost_maxima = true;
      }
    }

    if (value > 0 || check_for_lost_maxima) {
      // We are in a transcribed region.
      has_signal = 1;

      if (slope == 0) {
        // Slope ~ 0 -- this must be some type of plateu
        //Rprintf("  0 slope -- plateu\\n");
        if ((prev_slope == 0) && (lag_start == -1)) {
          lag_start = i;
        }
      } else {
        // not inside a plateau
        if (prev_value == 0) {
          // This index is the first non-zero element of a new event
          //Rprintf("  Starting new event\\n");
          if (i == 0) {
            // mitochondria screwed me w/ non zero value at pos 1!
            minima.push_back(1);
          } else {
            minima.push_back(i);
          }
        } else {
          // We are in the thick of some signal landscape, check if we are
          // climbing, reached a maxima, on the decline, hit a minima, etc.
          if (slope < 0) {
            if (prev_slope > 0) {
              // this is a maximum
              //Rprintf("  easy maximum...\\n");
              maxima.push_back(i);
            } else if ((prev_slope == 0) && (lag_start > 0) && (previous_turn != 1)) {
              // was in a plateu, but now we dropped off.
              // the check against previous_turn ensures no maxima is called
              // at the edge of a "cliff" that occurs after a maximum
              // set the maxima to be halfway between start of plateu and here
              //Rprintf(" negotiating maximum position from plateu. i:%d, lag_start;%d\\n", i, lag_start);
              maxima.push_back((int) (lag_start + ((i - lag_start - 1) / 2)));
            }
            previous_turn = 1;
          } else {
            // minimum? -- slope > 0
            if (prev_slope < 0) {
              if (previous_turn != -1) {
                minima.push_back(i);
                previous_turn = -1;
              }
            }
          }
        }
        lag_start = -1;
      } // if-else in a plateu
    } // if-else in a 0-valued element
    prev_value = value;
    prev_slope = slope;
  } // for loop

  return Rcpp::List::create(
    Rcpp::Named("maxima") = maxima,
    Rcpp::Named("minima") = minima);

  ', plugin="Rcpp")
}


##' Find quantile information for transcribed events in the signal.
##'
##' @param signal A coverage vector
##' @param bounds An IRanges instance defining the boundaries in signal
##' to look for quantile information from, or a number indicating the minimum
##' count to use for \code{slice}-ing the signal into events
defineQuantilePositions <- function(signal, bounds, quantile.breaks=NULL) {
  if (is.numeric(bounds) && length(bounds) == 1) {
    bounds <- slice(signal, bounds, includeLower=bounds > 0)
  }
  stopifnot(is(bounds, 'IRanges'))
  if (missing(quantile.breaks) || is.null(quantile.breaks)) {
    quantile.breaks <- c(0.05, 0.10, 0.15, 0.20, 0.80, 0.85, 0.90, 0.95)
  }
  if (any(quantile.breaks <= 0) || any(quantile.breaks >= 1)) {
    stop("Quantils must be within (0,1) range.")
  }
  if (is.unsorted(quantile.breaks)) {
    stop("Quantiles need to be in ascending order")
  }
  if (any(start(bounds) < 1L) || any(end(bounds) > length(signal))) {
    stop("Bounding info is out of bounds.")
  }
  if (!is.numeric(signal)) {
    signal <- as.numeric(signal)
  }

  bquants <- quantiles.indices.from.events.rcpp(signal, start(bounds),
                                                end(bounds), quantile.breaks)

  qdf <- do.call(rbind, bquants)
  colnames(qdf) <- paste('quantile', gsub("0\\.", "", quantile.breaks), sep=".")
  as(qdf, 'DataFrame')
}

if (FALSE && !exists('quantiles.indices.from.events.rcpp')) {

  ##' Calculates the indices of the quantile marks for the coverage events in `x`
  ##' which are start/end fenceposted by the start/end params.
  ##'
  ##'
  quantiles.indices.from.events.rcpp <- cxxfunction(signature(x="numeric", start="integer", end="integer", breaks="numeric"), '
    Rcpp::NumericVector xx(x);
    Rcpp::IntegerVector bstarts(start);
    Rcpp::IntegerVector bends(end);
    Rcpp::NumericVector percentiles(breaks);
    int np = percentiles.size();
    int nranges = bstarts.size();
    int istart;
    int iend;
    int iwidth;
    double value;
    std::vector<double> tmpvals;
    double total;
    double current_break;
    double cumsum;
    int i, j, nfound, tmp;
    Rcpp::List output(nranges);

    for (i = 0; i < nranges; i++) {
      tmpvals.clear();
      istart = bstarts[i] - 1;  // R -> C indexing
      iend = bends[i] - 1;      // R -> C indexing
      iwidth = iend - istart + 1;
      total = 0;
      for (j = 0; j < iwidth; j++) {
        value = xx[istart + j];
        total += value;
        tmpvals.push_back(value);
      }

      current_break = percentiles[0];
      j = 0;
      nfound = 0;
      std::vector<int> this_result;
      cumsum = 0;

      while ((j < iwidth) && (nfound < np)) {
        cumsum += tmpvals[j];
        if (current_break < (cumsum / total)) {
          this_result.push_back(j + 1); // C -> R indexing
          if (++nfound < np) {
            current_break = percentiles[nfound];
          }
        }
        j++;
      }

      if (this_result.size() < np) {
        if (this_result.size() == 0) {
          tmp = 1;
        } else {
          tmp = this_result[this_result.size() - 1];
        }
        for (j = this_result.size(); j < np; j++) {
          this_result.push_back(tmp);
        }
      }

      //this_result.push_back(total);
      output[i] = Rcpp::wrap(this_result);
    }

    return output;', plugin="Rcpp")
}

##' Finds the peaks/troughs in a (smoothed)siganl/coverage vector.
##'
##' @param signal A numeric or Rle coverage/signal vector.
##' @param bounds A length-1 IRanges-like object that gives the bounds that
##' this coverage vector lies within, ie: the position of of \code{signal[1]}
##' is \code{start(bounds)}
##'
##' @return An \code{IRanges} object with as many ranges as maxima found in
##' \code{signal}. start/end are the minima of the events. The positions are
##' reported with respect to their position in \code{signal}.
##'
##' \code{values(result)$peak.pos} has the position of the peak in the event
##'
##' \code{values(result)$shift} is the value to move the coordinates listed
##' in the result back to the "original" coordinates of the signal as indicated
##' in \code{bounds}. For example shift(result, values(result)$shift)
findPeakBoundsInSignal <- function(signal, bounds=NULL) {
  if (!is.null(bounds)) {
    bounds <- as(bounds, 'IRanges')
  }

  signal <- as.numeric(signal)
  len <- length(signal)
  turns <- turnfun(signal)
  peak.idx <- turns$maxima
  trough.idx <- turns$minima
  nonzero <- which(signal > 0)

  if (length(peak.idx) == 0) {
    ## uniform coverage over the entire event(?)
    ## or this weird one:
    ## ::::::::.......:::::::::
    peak.idx <- as.integer(floor(len / 2L))
    trough.idx <- integer()
  }

  ## Ensure all peaks have "surrounding" troughs
  if (length(trough.idx) == 0L) {
    ## There is only 1 peak
    trough.idx <- c(nonzero[1], tail(nonzero, 1))
  } else {
    if (trough.idx[1] > peak.idx[1]) {
      trough.idx <- c(nonzero[1], trough.idx)
    }
    if (tail(trough.idx, 1) < tail(peak.idx, 1)) {
      trough.idx <- c(trough.idx, tail(nonzero, 1))
    }
  }

  if (!all(findInterval(peak.idx, trough.idx), seq(peak.idx))) {
    stop("Peak/Trough ordering is hosney'd")
  }

  starts <- head(trough.idx, -1L) + 1L
  starts[1] <- starts[1] - 1L
  ends <- trough.idx[-1L] ## don't want to start at 2
  info <- IRanges(starts, ends)
  values(info) <- DataFrame(peak.pos=peak.idx)

  ## Add bounding info
  if (!is.null(bounds)) {
    values(info)$shift <- start(bounds[1]) - 1L
  }

  ## TODO: Marks maxima as significant.
  ##       We need  a good way to ignore "insignifant peaks" -- where we
  ##       find a local maximum, but it's too close to its surrounding troughs,
  ##       either in distance or height, look at::
  ##
  ## > expt <- nexpts$nick.mcf10a.s6
  ## > visualizeAntisense(expt, 'chr1', 197124317, 197124494, bsg,
  ##                      unique.only=FALSE, max.mismatch=2, bandwidth=18)
  ## and
  ## > visualizeAntisense(expt, 'chr1', 212989585, 212989680, bsg, unique.only=F,
  ##                      max.mismatch=2, bandwidth=18)

  info
}

findPeakBoundsInSignal.turnpoints <- function(signal, bounds=NULL) {
  if (!is.null(bounds)) {
    bounds <- as(bounds, 'IRanges')
  }

  signal <- as.numeric(signal)
  len <- length(signal)
  turns <- turnpoints(signal)
  peak.idx <- maxima(turns, as='index')
  trough.idx <- minima(turns, as='index')
  nonzero <- which(signal > 0)

  if (length(peak.idx) == 0) {
    ## uniform coverage over the entire event(?)
    ## or this weird one:
    ## ::::::::.......:::::::::
    peak.idx <- as.integer(floor(len / 2L))
    trough.idx <- integer()
  }

  ## Ensure all peaks have "surrounding" troughs
  if (length(trough.idx) == 0L) {
    ## There is only 1 peak
    trough.idx <- c(nonzero[1], tail(nonzero, 1))
  } else {
    if (trough.idx[1] > peak.idx[1]) {
      trough.idx <- c(nonzero[1], trough.idx)
    }
    if (tail(trough.idx, 1) < tail(peak.idx, 1)) {
      trough.idx <- c(trough.idx, tail(nonzero, 1))
    }
  }

  if (!all(findInterval(peak.idx, trough.idx), seq(peak.idx))) {
    stop("Peak/Trough ordering is hosney'd")
  }

  starts <- head(trough.idx, -1L) + 1L
  starts[1] <- starts[1] - 1L
  ends <- trough.idx[-1L] ## don't want to start at 2
  info <- IRanges(starts, ends)
  values(info) <- DataFrame(peak.pos=peak.idx)

  ## Add bounding info
  if (!is.null(bounds)) {
    values(info)$shift <- start(bounds[1]) - 1L
  }

  ## TODO: Marks maxima as significant.
  ##       We need  a good way to ignore "insignifant peaks" -- where we
  ##       find a local maximum, but it's too close to its surrounding troughs,
  ##       either in distance or height, look at::
  ##
  ## > expt <- nexpts$nick.mcf10a.s6
  ## > visualizeAntisense(expt, 'chr1', 197124317, 197124494, bsg,
  ##                      unique.only=FALSE, max.mismatch=2, bandwidth=18)
  ## and
  ## > visualizeAntisense(expt, 'chr1', 212989585, 212989680, bsg, unique.only=F,
  ##                      max.mismatch=2, bandwidth=18)

  info
}

################################################################################
## `peakedEvents` code below is a dumb way to find peaks. It assumes there is
## at most one dominant peak and 'slices' the event near the top of the peak
## to find other peaks close by in magnitude.
##
## Using the turnpoints method is a much better way to look at peak information
## since we can now easily find the peaks and troughs in an event
## (after smoothing)

## This misses a number of events whose width is larger than normal,
## but has lower expression -- the coverage of any given bp is < tpm3,
## but the tpm of the entire event is > 3tpm
##
## NOTES: "Peaking" the events might be able to help remove internal-priming
##        events, too: look at mcf10a region chr6:116,759,009-116,759,095
## DEBUG: look at `not.in.shared` in debug/peakedEvents.mcf10a.rda
##        note that peaks close together after the "peak calling" are a pain:
##        look at narrow diffs betweens successive peaks.
setGeneric("peakedEvents",
function(x, min.count, min.width=15L, .parallel=TRUE, ...) {
  standardGeneric("peakedEvents")
})

setMethod("peakedEvents", c(x="TailSeq"),
function(x, min.count, min.width, .parallel,
         chrs=levels(seqnames(x)), ...) {
  if (missing(min.count)) {
    min.count <- iTPM(3, x)
  }

  stopifnot(is.numeric(min.count))
  stopifnot(is.numeric(min.width))
  stopifnot(is.logical(.parallel))

  events <- llply(chrs, .parallel=.parallel, .inform=TRUE, function(chr) {
    cat("...", chr, "\n")
    reads <- getReadsOnChromosome(x, chr, smooth.by=NULL, ...)
    if (length(reads) == 0L) {
      return(GRanges())
    }
    peakedEvents(reads, min.count, min.width, .parallel=FALSE)
  })
  suppressWarnings(do.call(c, events))
})

setMethod("peakedEvents", c(x="GRanges"),
function(x, min.count, min.width, .parallel, ...) {
  if (length(x) == 0L) {
    return(x)
  }
  if (missing(min.count)) {
    if ('count' %in% colnames(values(x))) {
      min.count <- iTPM(3, sum(values(x)$count))
    } else {
      min.count <- iTPM(3, length(x))
    }
  }

  stopifnot(is.numeric(min.count))
  stopifnot(is.numeric(min.width))
  stopifnot(is.logical(.parallel))

  chrs <- names(seqlengths(x))
  events <- llply(chrs, .parallel=.parallel, function(chr) {
    reads <- x[seqnames(x) == chr]
    if (length(reads) == 0L) {
      return(GRanges())
    }
    stranded <- lapply(c('+', '-'), function(.strand) {
      these <- reads[strand(reads) == .strand]
      if (length(these) == 0L) {
        return(GRanges())
      }
      covr <- coverage(these)[[1]]
      ## `slice` is the primary means by which peaks are identified, but this
      ## can cause one "peak" to be split into two++ due to some
      ## vagaries in coverage so I'm trying to trim these out -- there's
      ## obviously a smarter way to do this.
      events <- slice(covr, lower=min.count)
      events <- events[width(events) >= min.width]
      peaked <- lapply(seq_along(events), function(idx) {
        view <- events[idx]
        this.covr <- view[[1L]]
        thresh <- floor(max(this.covr) - 1.5 * sd(this.covr))
        over <- as(this.covr >= thresh, 'IRanges')
        if (length(over) > 1L) {
          over <- over[which.max(width(over))]
        }
        bounds <- shift(over, start(view) - 1L)
      })
      peaked <- as(do.call(c, peaked), "IRanges")
      if (length(peaked) != length(events)) {
        stop("Number of peaks != number of events")
      }

      ## Assign a reasonable `count` to each of these peaks.
      GRanges(chr, peaked, .strand, count=viewMaxs(events))
    })
    do.call(c, stranded)
  })
  suppressWarnings(do.call(c, events))
})


## peakedEvents <- function(expt, min.count=iTPM(3, expt), upper.quantile=0.9,
##                          min.width=15L, .parallel=TRUE, ...) {
##   min.count <- as.integer(min.count)
##   events <- llply(chr=as.character(seqnames(expt))) %dopar% {
##     cat("...", chr, "\n")
##     reads <- getReadsOnChromosome(expt, chr, smooth.by=NULL, ...)
##     if (length(reads) == 0L) {
##       return(GRanges())
##     }
##     stranded <- lapply(c('+', '-'), function(.strand) {
##       these <- reads[strand(reads) == .strand]
##       if (length(these) == 0L) {
##         return(GRanges())
##       }
##       covr <- coverage(these)[[1]]
##       ## `slice` is the primary means by which peaks are identified, but this
##       ## can cause one "peak" to be split into two++ due to some
##       ## vagaries in coverage so I'm trying to trim these out -- there's
##       ## obviously a smarter way to do this.
##       events <- slice(covr, lower=min.count)
##       events <- events[width(events) >= min.width]
##       peaked <- lapply(seq_along(events), function(idx) {
##         view <- events[idx]
##         this.covr <- view[[1L]]
##         thresh <- floor(max(this.covr) - 1.5 * sd(this.covr))
##         over <- as(this.covr >= thresh, 'IRanges')
##         if (length(over) > 1L) {
##           over <- over[which.max(width(over))]
##         }
##         bounds <- shift(over, start(view) - 1L)
##       })
##       peaked <- as(do.call(c, peaked), "IRanges")
##       if (length(peaked) != length(events)) {
##         stop("Number of peaks != number of events")
##       }

##       ## Assign a reasonable `count` to each of these peaks.
##       GRanges(chr, peaked, .strand, count=viewMaxs(events))
##     })
##     do.call(c, stranded)
##   }
##   suppressWarnings(do.call(c, events))
## }
