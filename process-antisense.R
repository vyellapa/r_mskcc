## Look at events that span regions with reads from both strands
##   * What does the sense/antisense ratio look like?
##   * Can we automatically remove the "correct" artifact?
##
## chr1:563,654-572,031 is an amazing region:
##   what the hell is going on here?
##   looking at multimapped + 2mm reads, we have tags everywhere!
##
## This is an "easy" case of the "tricky" variety:
##   mcf10a : chr1:24,799,306-24,799,479
##
## A shifted antisense signal is caused by AAAA on the antisense strand?
##  chr1:32,830,599-32,831,321s
if (FALSE) {
  source('mayr/setup.R')
  expts <- loadMayr()
  expt <- expts[[1]]
  ip.dir <- 'mayr/diagnose/internal-priming'

  ## TODO: generate pictures of all "antisense areas" -- especially ones
  ##       that cover a wide range.
  ##       I think we're looking for a pattern that is a "wide event" where
  ##       the coverage on one strand is rather peaky/punctate, and the other
  ##       is wide/long range -- this is some signal of protocol "fail"
  anti.all <- senseAntisenseRatio(expt, unique.only=FALSE, max.mismatch=0L)
  anti.all <- transform(anti.all, width=end - start + 1L)
  anti.all <- transform(anti.all, log.ratio=log2(count.fwd) - log2(count.rev))

  anti.all.mid <- subset(anti.all, log.ratio > -1 & log.ratio < 1)
  anti.all.clear <- subset(anti.all, log.ratio < -1 | log.ratio > 1)
  high.all <- subset(anti.all, count.fwd + count.rev > 200)
  high.all.1 <- subset(high.all, ratio > 0.65 & ratio < 1.25)

  ##############################################################################
  ## Unique
  anti.u <- senseAntisenseRatio(expt, unique.only=TRUE, max.mismatch=0L)
  anti.u <- transform(anti.u, width=end - start + 1L)
  ## reads.all <- getReadsOnChromosome(expt, 'chr1', unique.only=FALSE, smooth.by=NULL)

  ## unique reads
  reads.u <- getReadsOnChromosome(expt, 'chr1', unique.only=TRUE, smooth.by=NULL)
}

##' @param expt The experiment to run the stats over
##' @param chrs The chromosomes to run this over (defaults to all in experiment)
##' @param dots send in values for unique.only and max.mismatch
senseAntisenseRatio <- function(expt, chrs=NULL, unique.only=FALSE, max.mismatch=0, ...) {
  if (is.null(chrs)) {
    chrs <- seqnames(expt)
  }
  ar <- foreach(chr=chrs, .packages="GenomicRanges") %dopar% {
    r <- getReadsOnChromosome(expt, chr, smooth.by=NULL, unique.only=unique.only,
                              max.mismatch=max.mismatch)
    covr <- coverage(r)[[1]]
    r.fwd <- r[strand(r) == '+']
    r.rev <- r[strand(r) == '-']

    event.bounds <- IntervalTree(slice(covr, 1, rangesOnly=TRUE))
    mm.fwd <- as.matrix(findOverlaps(ranges(r.fwd), event.bounds))
    mm.rev <- as.matrix(findOverlaps(ranges(r.rev), event.bounds))

    dt.fwd <- data.table(mm.fwd, key='subjectHits')
    dt.rev <- data.table(mm.rev, key='subjectHits')

    reads.per.fwd <- dt.fwd[, list(count=length(query)), by=key(dt.fwd)]
    reads.per.rev <- dt.rev[, list(count=length(query)), by=key(dt.rev)]

    antisense <- merge(reads.per.fwd, reads.per.rev, key='subject',
                       suffixes=c('.fwd', '.rev'))
    antisense$ratio <- with(antisense, count.fwd / count.rev)
    antisense <- transform(antisense, seqnames=chr,
                           start=start(event.bounds)[subject],
                           end=end(event.bounds)[subject])
    ## plot(density(log(antisense.regions$ratio)))
    antisense
  }
  do.call(rbind, ar)
}

##' @param x The object to get the reads from -- either a \code{GRanges} object
##' of reads already, or a SeqStore object to pull reads out of
##' @param chr Chromosome to reads off of (or a GRanges object that specifies
##' the bounds (chr,start,end)
##' @param start where to start to get reads off of on chr
##' @param end where to end getting reads off of chr
visualizeAntisense <- function(x, chr, start, end, bsg=NULL,
                               window.size=10, threshold=8,
                               kernel='normal', bandwidth=18, pad=50,
                               polya=kPolyA$dna.signal, unique.only=FALSE,
                               max.mismatch=NULL,
                               do.plot=c('smooth', 'peaks', 'prime', 'polya'),
                               ...) {
  in.bounds <- FALSE
  if (inherits(x, 'SeqStore')) {
    if (is.null(bsg)) {
      bsg <- getBsGenome(x)
    }
    x <- getReadsOnChromosome(x, chr, start=start, end=end, smooth.by=NULL,
                              unique.only=unique.only, max.mismatch=max.mismatch, ...)
    in.bounds <- TRUE
  }
  if (!inherits(x, 'GRanges')) {
    stop('Expecting GRanges object for `x`')
  }
  if (inherits(chr, 'GRanges')) {
    if (length(chr) > 1L) {
      stop("Bad GRanges object")
    }
    end <- end(chr)
    start <- start(chr)
    chr <- as.character(seqnames(chr))
  }

  if (!is.null(bsg)) {
    if (inherits(bsg, 'BSgenome')) {
      bsg <- unmasked(bsg[[chr]])
    }
    if (!inherits(bsg, 'DNAString')) {
      cat("Skipping nucelotide marking\n")
      bsg <- NULL
    }
  }

  bounds <- IRanges(start, end)
  if (!in.bounds) {
    x <- subsetByOverlaps(x, GRanges(chr[1], bounds, '*'))
  }
  if (length(x) == 0L) {
    cat("No reads in range\n")
    return(NULL)
  }

  fwd <- ranges(x[strand(x) == '+'])
  has.fwd <- length(fwd) > 0
  rev <- ranges(x[strand(x) == '-'])
  has.rev <- length(rev) > 0

  new.bounds <- IRanges(min(start(fwd), start(rev)) - pad,
                        max(end(fwd), end(rev)) + pad)
  shift.by <- start(new.bounds) - 1L

  if (shift.by > 0) {
    if (has.fwd) fwd <- shift(fwd, -shift.by)
    if (has.rev) rev <- shift(rev, -shift.by)
  }

  f.covr <- as.numeric(coverage(fwd))
  r.covr <- as.numeric(coverage(rev))
  covr <- coverage(shift(ranges(x), -shift.by))

  xlim <- c(start(new.bounds), end(new.bounds))
  ylim <- c(0, max(covr) * 1.1)

  .gen.x <- function(vals) {
    seq(xlim[1], length.out=length(vals))
  }

  plot(.gen.x(f.covr), f.covr, type='l', col='green', ylab="Coverage",
       xlim=xlim, ylim=ylim, main=paste("Coverage over", chr),
       xlab='Chromsome Coordinates')
  lines(.gen.x(r.covr), r.covr, col='red')
  lines(.gen.x(covr), covr, col='grey', lty='dashed')
  legend('topleft', legend=c('fwd', 'rev'), text.col=c('green', 'red'))
  if ('smooth' %in% do.plot && !is.null(kernel)) {
    if (has.fwd) {
      f.smooth <- kdeSmooth(f.covr, kernel=kernel, bandwidth=bandwidth, ...)
      lines(.gen.x(f.smooth), f.smooth, col='green', lty='dotted', lwd=2)

      if ('peaks' %in% do.plot) {
        ## turns
        peak.info <- findPeakBoundsInSignal(f.smooth)
        ## max
        idx <- values(peak.info)$peak.pos
        points(idx + xlim[1] - 1L, f.smooth[idx], col='green', pch=19)

        ## min
        idx <- c(start(peak.info)[1], end(peak.info))
        points(idx + xlim[1] - 1L, f.smooth[idx], col='green')
      }
    }

    if (has.rev) {
      r.smooth <- kdeSmooth(r.covr, kernel=kernel, bandwidth=bandwidth, ...)
      lines(.gen.x(r.smooth), r.smooth, col='red', lty='dotted', lwd=2)

      if ('peaks' %in% do.plot) {
        ## turns
        peak.info <- findPeakBoundsInSignal(r.smooth)
        ## max
        idx <- values(peak.info)$peak.pos
        points(idx + xlim[1] - 1L, r.smooth[idx], col='red', pch=19)

        ## min
        idx <- c(start(peak.info)[1], end(peak.info))
        points(idx + xlim[1] - 1L, r.smooth[idx], col='red')
      }
    }
  }

  if (!missing(bsg)) {
    vf <- Views(bsg, new.bounds)
    lf <- letterFrequencyInSlidingView(vf, window.size, 'A')[[1]]
    mark.f <- which(lf >= threshold)
    if ('prime' %in% do.plot && length(mark.f) > 0) {
      rug(mark.f + shift.by, col='green')
    }

    vr <- Views(bsg, shift(new.bounds, -window.size))
    lr <- letterFrequencyInSlidingView(vr, window.size, 'T')[[1]]
    mark.r <- which(lr >= threshold)
    if ('prime' %in% do.plot && length(mark.r) > 0) {
      rug(mark.r + shift.by, col='red', side=3)
    }

    if ('polya' %in% do.plot && !is.null(polya)) {
      hits.fwd <- markPolyA(vr, polya)
      hits.rev <- markPolyA(vr, polya, strand='-')
      pa.alpha <- seq(1, .3, length.out=length(polya))
      if (!is.null(hits.fwd)) {
        col.fwd <- rgb(1, 0, 0, alpha=pa.alpha[hits.fwd$rank])
        abline(v=hits.fwd$start, col=col.fwd)
      }
      if (!is.null(hits.rev)) {
        col.rev <- rgb(0, 1, 0, alpha=pa.alpha[hits.rev$rank])
        abline(v=hits.rev$start, col=col.rev)
      }
    }
  }

  ret <- list(fwd=f.covr, rev=r.covr, xlim=xlim)
  if ('smooth' %in% do.plot && !is.null(kernel)) {
    ret$f.smooth <- if (has.fwd) f.smooth else NULL
    ret$r.smooth <- if (has.rev) r.smooth else NULL
  }

  invisible(ret)
}

##' Returns a list of integer vectors indicating where the polya signal starts
##' for each signal in \code{polya}
##'
##' Assumes that sequence and polya are given in "forward" sequence space. The
##' matching is done on reverse complement of the signal if (i) do.rc is TRUE; or
##' (ii) do.rc is not set and strand == '-'
##'
##' Positions are always reported by their "left-most" position on the (+) strand.
##' The signals are assumed to be in "rank" order, a column of rank is returned
##' to indicate signal strength for potential downstream use.
##' The data.frame returned is also ordered by increasing hit position if
##' reorder==TRUE.
markPolyA <- function(sequence, polya=kPolyA$dna.signal, max.mismatch=0L,
                      strand='+', do.rc=NULL, reorder=TRUE) {
  if (is.null(polya) || length(polya) == 0) {
    return(NULL)
  }
  hits <- lapply(1:length(polya), function(idx) {
    pa <- as(polya[[idx]], 'DNAString')
    if (is.null(do.rc)) {
      if (strand == '-') {
        pa <- reverseComplement(pa)
      }
    } else {
      if (do.rc) {
        pa <- reverseComplement(pa)
      }
    }
    starts <- start(matchPattern(pa, sequence))
    if (length(starts) > 0L) {
      data.frame(start=starts, pattern=as.character(pa), rank=idx)
    } else {
      NULL
    }
  })

  hits <- do.call(rbind, hits)
  if (!is.null(hits) && reorder) {
    hits <- hits[order(hits$start),]
  }
  hits
}

if (FALSE) {
  ## chr1:182,856,526-182,856,897
  target <- IRanges(182856526, 182856897)
  fwd <- shift(subsetByOverlaps(ranges(r.fwd), target), -(start(target[1])-10))
  rev <- shift(subsetByOverlaps(ranges(r.rev), target), -(start(target[1])-10))
  fwd.covr <- coverage(fwd)
  rev.covr <- coverage(rev)
  plot(fwd.covr, type='l', col='green')
  lines(rev.covr, col='red')
}
