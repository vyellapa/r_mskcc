combinePrimingInfo <- function(evts, ip.stats, pa.stats,
                               event.a.rich=0.65, max.a.rich=0.8,
                               rescue.polya=kPolyA$dna.signal[1:6],
                               pa.max.from.end=45L, pa.min.from.end=10L,
                               anno.rescue=NULL, rescue.from.end=40,
                               as.mask=FALSE) {
  ##    The default tresholds to flag events for priming is:
  ##      event.a.rich (iArich) >= 0.65
  ##      a hit from externally.primed.window (ewindow.primed)
  if (is(evts, "SummarizedExperiment")) {
    se <- evts
    evts <- rowData(se)
  } else {
    se <- NULL
  }
  stopifnot(inherits(evts, 'GenomicRanges'))
  stopifnot(inherits(ip.stats, 'DataFrame') && nrow(ip.stats) == length(evts))
  if (is.character(rescue.polya) || isTRUE(rescue.polya)) {
    stopifnot(inherits(pa.stats, 'data.frame'))
  }

  meta <- values(evts)

  rename <- c(externally.primed.window='external',
              externally.primed.abs='eA6',
              internally.primed.primed='internal')

  for (name in names(rename)) {
    change.idx <- which(names(ip.stats) == name)
    if (length(change.idx) == 1L) {
      names(ip.stats)[change.idx] <- rename[name]
    }
  }

  if (nrow(ip.stats) != length(evts)) {
    stop("events do not match ip stats")
  }

  ip.stats$external[is.na(ip.stats$external)] <- 0L
  ip.stats$internal[is.na(ip.stats$internal)] <- 0L

  window.ip.flag <- ip.stats$external > 0
  internal.ip.flag <- ip.stats$internal > max.a.rich

  ## Add pA signal stats
  pas <- identifyValidPaSignals(evts, pa.stats,
                                pa.min.from.end=pa.min.from.end,
                                pa.max.from.end=pa.max.from.end)

  ip.stats$n.signal <- rep(0L, nrow(ip.stats))
  ip.stats$n.signal[pas$query] <- pas$N

  ip.stats$signal <- rep(NA_character_, nrow(ip.stats))
  ip.stats$signal[pas$query] <- as.character(pas$first)

  priming.status <- rep('none', length(evts))

  a.rich <- ip.stats$internal >= event.a.rich
  priming.status[a.rich] <- 'a.rich'

  external.p <- ip.stats$external > 0
  priming.status[external.p] <- 'external'

  ## Finally, if any event has an internal/internal A score >= max.a.rich
  ## nuke it w/o regards for poly(A) rescue
  a.rich.max <- ip.stats$internal >= max.a.rich
  priming.status[a.rich.max] <- 'a.rich.max'

  ip.stats$ip.status <- priming.status

  axe <- ip.stats$ip.status != 'none' & ip.stats$n.signal == 0
  axe <- axe | a.rich.max

  if (!is.null(anno.rescue)) {
    on.end <- flagEventsOnAnnotedGeneEnd(evts, anno.rescue, rescue.from.end)
    ip.stats$annotated.end <- on.end > 0
    axe <- axe & on.end == 0
  }

  ip.stats$ip.axe <- axe
  ip.stats
}

countValidPaSignals <- function(x, pa.stats, pa.min.from.end=10L,
                                pa.max.from.end=45L,
                                pa.rank=kPolyA$dna.signal,
                                ...) {
  stopifnot(inherits(x, 'GenomicRanges'))
  stopifnot(inherits(pa.stats, 'data.frame'))
  stopifnot(pa.max.from.end > pa.min.from.end)

  identifyValidPaSignals(x, pa.stats,
                         pa.min.from.end=pa.min.from.end,
                         pa.max.from.end=pa.max.from.end,
                         count.only=TRUE, )
}

## TODO: Change this function to return "unwound" pA info
##       Currently we get a result like so:
##
##       queryHits N  first second third fourth
##               1 1 ATTAAA     NA    NA     NA
##               4 1 CATAAA     NA    NA     NA
##               5 1 AATAGA     NA    NA     NA
##               6 1 GATAAA     NA    NA     NA
##               9 1 ATTAAA     NA    NA     NA
##              10 1 TTTAAA     NA    NA     NA
##              11 2 AATAAA ATTAAA    NA     NA
##              19 1 AATGAA     NA    NA     NA
##              20 2 AATAAA TTTAAA    NA     NA
##              22 2 AATAAA AATACA    NA     NA
##
##       But we probably want:
##       querHits  N  AATAAA ATTAAA TATAAA AGTAAA ...
##              1  5       3      1      0      1 ...
##
##       Can you dig it?
identifyValidPaSignals <- function(x, pa.stats, pa.min.from.end=10L,
                                   pa.max.from.end=45L,
                                   count.only=FALSE,
                                   pa.rank=kPolyA$dna.signal, ...) {
  stopifnot(inherits(x, 'GenomicRanges'))
  stopifnot(inherits(pa.stats, 'data.frame'))
  stopifnot(pa.max.from.end > pa.min.from.end)

  pas <- transform(pa.stats, start=start.signal - 1L, end=end.signal-1L, width=6)
  pa.gr <- as(pas, "GRanges")
  pa.gr <- GenomicRanges::shift(pa.gr, values(pa.gr)$start.event)

  x.gr <- resize(x, width=pa.max.from.end - pa.min.from.end, fix='end')
  is.fwd <- as.logical(strand(x.gr) == '+')
  x.gr <- GenomicRanges::shift(x.gr, ifelse(is.fwd, -pa.min.from.end, pa.min.from.end))

  if (count.only) {
    ans <- countOverlaps(x.gr, pa.gr)
  } else {
    mm <- as.matrix(findOverlaps(x.gr, pa.gr))
    mm <- transform(data.table(mm, key="queryHits"),
                    signal=values(pa.gr)$signal[subjectHits])

    mm$signal.rank <- match(mm$signal, pa.rank)
    default <- NA_character_

    ans <- mm[, {
      N <- length(signal)
      rank <- order(signal.rank)
      first <- as.character(signal[rank[1L]])
      second <- default
      third <- default
      fourth <- default
      if (N > 1L) {
        second <- as.character(signal[rank[2L]])
      }
      if (N > 2L) {
        third <- as.character(signal[rank[3L]])
      }
      if (N > 3L) {
        fourth <- as.character(signal[rank[4L]])
      }
      list(N=length(signal), first=first, second=second, third=third,
           fourth=fourth)
    }, by='queryHits']
  }

  ans
}

identifyCleavagePoint <- function(x, bsg, cleavage.site='CA',
                                  cleavage.window=5L) {
  stopifnot(inherits(x, 'GenomicRanges'))
  stopifnot(is(bsg, "BSgenome"))
  cleavage.window <- as.integer(cleavage.window)

  chrs <- unique(as.character(seqnames(x)))
  cs.counts <- integer(length(x))

  cleave.fwd <- DNAString(cleavage.site)
  cleave.rev <- reverseComplement(cleave.fwd)

  is.fwd <- as.logical(strand(x) == '+')
  is.rev <- !is.fwd

  cc <- mclapply(chrs, function(chr) {
    bchr <- unmasked(bsg[[chr]])
    this.chr <- as.logical(seqnames(x) == chr)

    idxs.fwd <- NULL
    idxs.rev <- NULL
    counts.fwd <- NULL
    counts.rev <- NULL

    for (.strand in c('+', '-')) {
      fix <- if (.strand == '+') 'end' else 'start'
      cleave.at <- if (.strand == '+') cleave.fwd else cleave.rev
      this.strand <- if (.strand == '+') is.fwd else is.rev
      idxs <- which(this.chr & this.strand)

      if (length(idxs) > 0) {
        bounds <- resize(ranges(x[idxs]), width=1L, fix=fix) + cleavage.window
        dna <- Views(bchr, bounds)
        counts <- vcountPattern(cleave.at, dna)

        if (.strand == '+') {
          idxs.fwd <- idxs
          counts.fwd <- counts
        } else {
          idxs.rev <- idxs
          counts.rev <- counts
        }
      }
    }

    list(idxs.fwd=idxs.fwd, counts.fwd=counts.fwd,
         idxs.rev=idxs.rev, counts.rev=counts.rev)
  })

  for (cstats in cc) {
    if (length(cstats$idxs.fwd) > 0) {
      cs.counts[cstats$idxs.fwd] <- cstats$counts.fwd
    }
    if (length(cstats$idxs.rev) > 0) {
      cs.counts[cstats$idxs.rev] <- cstats$counts.rev
    }
  }

  cs.counts
}

flagEventsOnAnnotedGeneEnd <- function(events, annotated.genome,
                                       target.width=40) {
  anno <- as(annotated.genome, 'GRanges')
  anno <- anno[values(anno)$exon.anno == 'utr3']
  anno <- resize(anno, width=1, fix='end') + floor(target.width / 2)

  anno.dt <- as(anno, 'data.table')
  key(anno.dt) <- c('seqnames', 'entrez.id', 'strand', 'start')
  ends <- anno.dt[, {
    idx <- if (strand[1] == '+') which.max(end) else which.min(start)
    list(start=start[idx], end=end[idx], symbol=symbol[1])
  }, by=list(seqnames, entrez.id, strand)]

  countOverlaps(events, as(ends, "GRanges")) > 0
}
