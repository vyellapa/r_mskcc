## TODO: Remove poly-A stuff from the TagSeq package
##
## Information about poly(a) signals taken from::
##   Tian B, et al. A large-scale analysis of mRNA polyadenylation of human and
##   mouse genes. NAR, 2005. (Table 2)
kPolyA <- data.frame(
  signal=c("AAUAAA","AUUAAA","UAUAAA","AGUAAA","AAGAAA","AAUAUA","AAUACA",
           "CAUAAA","GAUAAA","AAUGAA","UUUAAA","ACUAAA","AAUAGA"),
  hs.frequency=c(53.18, 16.78, 4.37, 3.72, 2.99, 2.13, 2.03, 1.92,
                 1.75, 1.56, 1.20, 0.93, 0.60),
  mm.frequency=c(59.16, 16.11, 3.79, 3.28, 2.15, 1.71, 1.65, 1.80,
                 1.16,0.90,1.08,0.64,0.36),
  stringsAsFactors=FALSE)
kPolyA$dna.signal <- gsub("U", "T", kPolyA$signal)

##' Collect information on polyA signals for each event
##'
##' Every time a signal matches within a range in \code{event}, there
##' will be a corresponding row produced with the matching information.
##' For instance, if the same signal matches in two different places within
##' the same event, there will be two rows generated for that signal with
##' the matching information.
##'
##' NOTE: This function only looks WITHIN the events -- the combinePrimingInfo
##'       function measures from the cleavage site to validated the pAS, but
##'       it seems that we should be doing that here(!)
##'
##' If you just want one pA signal per event, follow this call with
##' \code{summarizePolyAPositionStatistics}
##'
##' @param events A GRanges object. The ranges of which will be used
##' to define the regions in which to look for poly(A) signals
##' @param bsg A BSgenome object to get sequence information from.
##' @param polya A vector of poly(A) signals (or any, really) to look
##' for inside \code{events}.
##' @param adjust.bounds A function used on \code{events} to readjust
##' its boundaries when looking for the poly(A) signal.
##' @param ... Additional parameters to pass to the \code{adjust.bounds}
##' function.
##'
##' @return A \code{data.frame} object with mathcing information.
##'
##' \describe{
##'   \item{}{
##'     \code{{start|end}.event} match the start/ends of \code{events}
##'   }
##'   \item{} {
##'     \code{{start|end}.signal} are the positions inside the event
##'     where the signal is. For instance:
##'       bsg.chr[start.event:end.event][start.signal:end.signal]
##'     is where the signal is
##'   }
##'   \item{}{
##'     \code{signal} The signal that is found in this result
##'   }
##'   \item{}{ ... there's more but you don't care ...}
##' }
collectPolyAPositionStatistics <-
  function(events, bsg, polya=kPolyA$dna.signal, cleavage.site='CA',
           cleavage.window=3L, .parallel=FALSE) {
  if (!inherits(events, 'GRanges')) {
    events <- as(events, 'GRanges')
  }
  chrs <- unique(as.character(seqnames(events)))

  pa.fwd <- DNAStringSet(polya)
  pa.rev <- reverseComplement(pa.fwd)
  cleave.fwd <- DNAString(cleavage.site)
  cleave.rev <- reverseComplement(cleave.fwd)

  if (.parallel && length(chrs) > 1L)  {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }

  all.pos <- foreach(chr=chrs, .packages="Biostrings",
                     .options.multicore=list(preschedule=FALSE)) %loop% {
    cat("---", chr, "---\n")
    bsg.chr <- unmasked(bsg[[chr]])
    e.chr <- events[seqnames(events) == chr]
    if (length(e.chr) == 0L) {
      return(NULL)
    }
    str.pos <- lapply(c('+', '-'), function(s) {
      ## Finding PolyA signals in events on chromosome one strand at a time
      e.strand <- e.chr[strand(e.chr) == s]
      if (length(e.strand) == 0L) {
        return(NULL)
      }

      if (s == '+') {
        pa.pattern <- pa.fwd
        cleave.at <- cleave.fwd
      } else {
        pa.pattern <- pa.rev
        cleave.at <- cleave.rev
      }

      dna <- as(Views(bsg.chr, ranges(e.strand)), 'XStringSet')

      ## Iterate over the polya signals, and try to match them against the
      ## selected genomic sequence
      matches <- lapply(1:length(pa.pattern), function(idx) {
        pattern <- pa.pattern[[idx]]
        the.signal <- as.character(pa.fwd[[idx]])
        m <- vmatchPattern(pattern, dna)
        mdf <- as.data.frame(m)
        if (nrow(mdf) == 0L) {
          return(NULL)
        }
        mdf <- transform(mdf, signal=the.signal, seqnames=chr, strand=s)
        colnames(mdf) <- gsub('start', 'start.signal', colnames(mdf))
        colnames(mdf) <- gsub('end', 'end.signal', colnames(mdf))
        mdf
      })
      matches <- do.call(rbind, matches)

      if (!is.null(matches) && nrow(matches) > 0) {
        meta <- values(e.strand)
        matches$start.event <- start(e.strand)[matches$group]
        matches$end.event <- end(e.strand)[matches$group]
        if ('count' %in% colnames(meta)) {
          matches$count <- meta$count[matches$group]
        }
        if ('entrez.id' %in% colnames(meta)) {
          matches$entrez.id <- meta$entrez.id[matches$group]
        }
        if ('exon.anno' %in% colnames(meta)) {
          matches$exon.anno <- meta$exon.anno[matches$group]
        }
        if (s == '+') {
          matches$dist.from.tx.start <- matches$start.signal - 1L
        } else {
          event.width <- matches$end.event - matches$start.event + 1L
          matches$dist.from.tx.start <- event.width - matches$end.signal
        }
      } else {
        matches <- NULL
      }

      axe <- which(colnames(matches) == 'space')
      if (length(axe) > 0) {
        matches <- matches[, -axe, drop=FALSE]
      }
      rownames(matches) <- NULL


      #########################################################################
      ## How many cleavage dimers are in the allotted window in the nieghborhood
      ## around the cleavage site?
      ## if (!is.null(matches) && nrow(matches) > 0) {
      ##   cleavage.range <- IRanges(matches$start.event, matches$end.event)
      ##   if (s == '+') {
      ##     cleavage.range <- resize(cleavage.range, width=1, fix='end')
      ##   } else {
      ##     cleavage.range <- resize(cleavage.range, width=1, fix='start')
      ##   }
      ##   cleavage.range <- cleavage.range + cleavage.window
      ##   dna <- as(Views(bsg.chr, cleavage.range), 'XStringSet')
      ##   matches$cleavage.signals <- vcountPattern(cleave.at, dna)
      ## }
      ##
      matches
    })
    ## all.pos[[chr]] <- do.call(rbind, str.pos)
    do.call(rbind, str.pos)
  }

  do.call(rbind, all.pos)
}

## pa.stats has a `dist.from.tx.start` column we can juse use to check
## that hte poly(A) signal is "in bounds"
isRealPolyA <- function(pa.stats, min.from.end=10L, max.from.end=45L, ...) {
  with(pa.stats,
       dist.from.tx.start >= min.from.end &
       dist.from.tx.start <= max.from.end)
}

##' Poly(A) signals found in "suspect" regions of events are
##' removed.
##'
##' As of 2011-05-05, suspect regions are defined as:
##'
##'  (i) Signals found < 10bp from the end of the transcribed region
##' (ii) Signals found at or 3' to regions flagged as "internally primed"
##'      for the event. We now find the 5'-most region that can be
##'      flagged as internally primed and report that position in
##'      values(events)$is.primed
##'
##' @param pa.stats A data.frame-like object that has the poly(A) stats
##' as calculated in summarizePolyAPositionStatistics
##' @param events A GRanges/data.frame object that has the information
##' for the "humps" in the data (event start/end and peak.pos)
##' @param min.from.end The minimum distance from the end (3') of the
##' event that the poly(A) signal must be in order to be considered "valid".
##' @param max.from.end The maximum distance from the end of the event
##' that the poly(A) signal can be in in order to "count" as usable.
##' @param mask.primed If TRUE, then the poly(A) signal must be further
##' upstream than the most 5' basepair that is above our internal-priming
##' threshold. This requires a "primed.pos" column for \code{events}.
##'
##' NOTE: The defualt value for max.from.end was 35bp when the hump-rescue-1
##'       datasets were generated (2011-05-11)
##'
##' @return A logical vector indicating whether each poly(A) signal listed
##' in \code{pa.stats} lands in a region that we believe a poly(A) signal
##' should be in.
isRealPolyA.overlaps <- function(pa.stats, events, min.from.end=10L,
                                 max.from.end=45L, mask.primed=FALSE, ...) {
  ## A "mask" GRanges object will be built and then used to find
  ## the remove
  ## any poly(A) signals that overlap with it.
  ## mask <- resize(events, width=min.from.end, fix='end')

  bad <- width(events) < max.from.end
  if (any(bad)) {
    bad.events <- events[bad]
    events <- events[!bad]
  }

  good <- resize(events, width=max.from.end, fix='end')
  good <- resize(good, width=max.from.end - min.from.end + 1L, fix='start')

  if (mask.primed) {
    ## Make masks wider (up stream) if there is a flagged internal priming
    ## event upstream from its current tx.start
    starts <- start(good)
    ends <- end(good)
    primed.pos <- values(good)$is.primed

    is.fwd <- as.logical(strand(events) == '+')
    is.rev <- !is.fwd
    is.primed <- as.logical(primed.pos)
    primed.fwd <- is.primed & is.fwd
    primed.rev <- is.primed & is.rev

    ends[primed.fwd] <- pmin(primed.pos[primed.fwd], ends[primed.fwd])
    starts[primed.rev] <- pmax(primed.pos[primed.rev], starts[primed.rev])

    good.keep <- ends >= starts
    good <- good[good.keep]
    start(good) <- starts[good.keep]
    end(good) <- ends[good.keep]
  }

  ## Turn the position of the poly(A) signals into GRanges object
  pa.starts <- with(pa.stats, start.event + start.signal -1L)
  pa.ends <- pa.starts + pa.stats$width -1L
  pa <- GRanges(pa.stats$seqnames, IRanges(pa.starts, pa.ends),
                pa.stats$strand)

  pa <- resize(as(pa, "GRanges"), width=1L, fix='center')
  yes <- countOverlaps(pa, good)
  yes > 0L
}

locatePolyA <- function(pa.stats, events, min.from.end=10L,
                        max.from.end=45L, mask.primed=FALSE,
                        use.quantiles=FALSE, qstart='quantile.start',
                        qend='quantile.end', ...) {
  ## A "mask" GRanges object will be built and then used to find
  ## the remove
  ## any poly(A) signals that overlap with it.
  ## mask <- resize(events, width=min.from.end, fix='end')

  ## It seems that most all of these are internally primed!
  ## but let's ignore that for now
  ## bad <- width(events) < max.from.end
  ## if (any(bad)) {
  ##   bad.events <- events[bad]
  ##   events <- events[!bad]
  ## }
  max.from.end <- pmin(width(events), max.from.end)
  final.width <- max.from.end - min.from.end + 1L
  too.narrow <- final.width < 1
  if (any(too.narrow)) {
    cat(sum(too.narrow), "events are too narrow\n")
    final.width[too.narrow] <- 1L
    ## events <- events[!too.narrow]
    ## max.from.end <- max.from.end[!too.narrow]
    ## max.from.end[too.narrow] <-
  }

  good <- resize(events, width=max.from.end, fix='end')
  good <- resize(good, width=final.width, fix='start')

  if (use.quantiles) {
    if (!all(c(qstart, qend)) %in% colnames(values(events))) {
      stop("Can't find quantile info")
    }
    qbounds <- IRanges(values(events)[[qstart]], values(events)[[qend]])
    is.fwd <- as.logical(strand(good) != '-')
    is.rev <- !is.fwd
    ends <- end(good)
    starts <- stat(good)

    ends[is.fwd] <- pmin(ends[is.fwd], end(qbounds[is.fwd]))
    starts[is.rev] <- pmax(starts[is.rev], starts(qbounds[is.rev]))
    good.keep <- ends >= starts
    good <- good[good.keep]
    start(good) <- starts[good.keep]
    end(good) <- ends[good.keep]
  }

  ## Turn the position of the poly(A) signals into GRanges object
  pa.starts <- with(pa.stats, start.event + start.signal -1L)
  pa.ends <- pa.starts + pa.stats$width -1L
  pa <- GRanges(pa.stats$seqnames, IRanges(pa.starts, pa.ends),
                pa.stats$strand)

  pa.slim <- resize(as(pa, "GRanges"), width=1L, fix='center')

  uo <- assignUniqueOverlaps(pa.slim, events)
  take.which <- list(which.max, which.min)
  dt <- data.table(pa.idx=seq_along(uo), event.idx=uo,
                   take=as.vector(ifelse(strand(pa.slim) == '+', 1L, 2L)),
                   signal=pa.stats$signal, start=start(pa),
                   key='event.idx')
  dt <- subset(dt, !is.na(event.idx))

  event.summary <- dt[, {
    idx <- take.which[[take[1]]](start)
    list(polya.count=length(start), polya.start=start[idx],
         polya.signal=signal[idx])
  }, by=event.idx]

  event.summary <- as.data.frame(event.summary)
  event.summary$polya.signal <- as.character(event.summary$polya.signal)
  result <- lapply(event.summary[,-1], function(column) {
    fn <- getFunction(class(column))
    wut <- fn(length(events))
    wut[event.summary$event.idx] <- column
    wut[too.narrow] <- fn(sum(too.narrow))
    wut
  })

  as.data.frame(result)
}

##' Trim event bounds to look for poly(A) signals.
##'
##'  (i) If the event isn't primed, the poly(A) signal must be
##'      \code{min.from.end=10}bp upstream from end of event
##' (ii) If it is primed, the poly(A) signal must be 5' of the position
##'      indicated in values(events)$is.primed
trimEventsForPolyASearch <- function(events, min.from.end=10L, ...) {
  stop("Don't adjust the events, but rather use removeSuspectPolyA ",
       "to remove poly(A) signals from consideration")
  if (min.from.end > 0L) {
    new.widths <- ifelse(width(events) < (min.from.end + 1L), 1L,
                         width(events) - min.from.end)
    events <- resize(events, width=new.widths)
  }

  starts <- start(events)
  ends <- end(events)
  primed.pos <- values(events)$is.primed

  is.fwd <- as.logical(strand(events) == '+')
  is.rev <- !is.fwd
  is.primed <- as.logical(primed.pos)
  primed.fwd <- is.primed & is.fwd
  primed.rev <- is.primed & is.rev

  ends[is.fwd] <- ifelse(is.primed[is.fwd], primed.pos[is.fwd],
                         ends[is.fwd])
  starts[is.rev] <- ifelse(is.primed[is.rev], primed.pos[is.rev],
                           starts[is.rev])


  replace.fwd <- ends[primed.fwd] > primed.pos[primed.fwd]
  replace.rev <- starts[primed.rev] < primed.pos[primed.rev]
  if (any(replace.fwd)) {
    ends[is.fwd][replace.fwd] <- primed.pos[is.fwd][replace.fwd]
  }
  if (any(replace.rev)) {
    starts[is.rev][replace.rev] <- primed.pos[is.rev][replace.rev]
  }

  start(events) <- starts
  end(events) <- ends

  events
}

##' Picks just one poly(A) signal from each "event".
##'
##' This function works on the data.frame that is returned from
##' \code{collectPolyAPositionStatistics}.
##'
##' "events" are uniquefly defined by the seqnames/strant/start.event/end.event
##' columns.
##'
##' @param pa.stats \code{data.frame} with matching info, as returned
##' @param take Which signal to take, the (event) "proximal", "distal",
##' or one in the center(?)
##' @param polya A character vector of signals to use (if you don't want)
##' to use all of them listed in \code{pa.stats}.
##' @param start.event The name of the column in pa.stats that indicates
##' the start bounds of the event
##' @param end.event The name of the column in pa.stats that indicates
##' the end bounds of the event
summarizePolyAPositionStatistics <-
function(pa.stats, take='proximal', polya=NULL, start.event='start.event',
         end.event='end.event') {
  take <- match.arg(take, c('proximal', 'distal', 'center', 'all'))
  if (!is.null(polya)) {
    pa.stats <- subset(pa.stats, signal %in% polya)
  }
  pa.stats <- as.data.table(pa.stats)
  ## Adding dist.from.tx.start to key orders the data from distal
  ## to proximal.
  setkeyv(pa.stats, c('seqnames', 'strand', start.event, end.event,
                      'dist.from.tx.start'))
  iter.by <- head(key(pa.stats), -1L)
  if (take == 'center') {
    ## this is uber slow -- is it .SD's fault?
    pa.summary <- pa.stats[, {
      idx <- max(1L, floor(length(signal) / 2L))
      lapply(.SD, '[', idx)
    }, by=iter.by]
  } else {
    ## fancyer data.table indexing
    idx.dt <- pa.stats[, iter.by, with=FALSE]
    setkeyv(idx.dt, head(key(pa.stats), -1L))
    idx.dt <- unique(idx.dt)
    .mult <- if (take == 'proximal') 'first' else 'last'
    pa.summary <- pa.stats[idx.dt, mult=.mult]
  }

  pa.summary
}

plotPolyAFrequency <- function(pa.stats, take='all', polya=NULL,
                               polya.order=kPolyA$dna.signal,
                               expand.count=TRUE,
                               title.append=NULL) {
  if (take != 'all' || !is.null(polya)) {
    pa.stats <- summarizePolyAPositionStatistics(pa.stats, take, polya)
  }
  pa.stats <- as.data.table(pa.stats)
  setkeyv(pa.stats, 'signal')

  if (!expand.count) {
    pa.stats$count <- 1L
  }

  counts <- pa.stats[, list(count=sum(count)), by='signal']

  title <- 'PolyA signals found in transcribed events'
  if (take != 'all') {
    title <- paste('Most', take, title)
  }
  if (!is.null(title.append)) {
    title <- paste(title, title.append)
  }

  if (!is.null(polya.order)) {
    xref <- match(counts$signal, polya.order)
    if (any(is.na(xref))) {
      warning("polya reordoring won't work")
    } else {
      sxref <- sort(xref)
      counts$signal <- factor(counts$signal, levels=polya.order[sxref])
    }
  }

  g <- ggplot(as.data.frame(counts), aes(signal, count)) +
    geom_bar(aes(fill=signal), stat='identity') +
      ylab("Count") + xlab("PolyA signal") + theme_bw() +
        opts(title=title, axis.text.x=theme_text(angle=-45, hjust=0, vjust=1),
             legend.text=theme_text(size=16)) + opts(axis.text.x=theme_blank())
  print(g)
  invisible(pa.stats)
}

## See event-chars vignette for another take
plotPolyADistanceToCleavage <- function(pa.stats, take='all', polya=NULL,
                                        polya.order=kPolyA$dna.signal,
                                        expand.count=TRUE, title.append) {
  if (take != 'all' || !is.null(polya)) {
    pa.stats <- summarizePolyAPositionStatistics(pa.stats, take, polya)
  }

  setkeyv(pa.stats, 'signal')

  if (!expand.count) {
    pa.stats$count <- 1L
  }

  x <- pa.stats[, list(signal=rep(signal, count),
                       dist=rep(dist.from.tx.start, count),
                       exon.anno=rep(exon.anno, count))]

  g <- ggplot(as.data.frame(x), aes(x=dist)) + geom_density(fill=NA) +
    facet_grid(exon.anno ~ signal, scale="free")
}



