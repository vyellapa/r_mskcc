setGeneric("countTagsByRegion", signature=c("x"),
function(x, regions, ignore.strand=FALSE, assign.by=c('unique', 'all'),
         .parallel=TRUE, ...) {
  standardGeneric("countTagsByRegion")
})

setMethod("countTagsByRegion", c(x="GenomicRanges"),
function(x, regions, ignore.strand, assign.by, .parallel, ...) {
  assign.by <- match.arg(assign.by, c('unique', 'all'))
  if (assign.by == 'unique') {
    assignUniqueOverlaps(x, regions)
  } else {
    countOverlaps(regions, x)
  }
})

setMethod("countTagsByRegion", c(x="BamFile"),
function(x, regions, ignore.strand, assign.by, .parallel, ...) {
  assign.by <- match.arg(assign.by, c('unique', 'all'))
  .countTagsByRegion(x, regions, ignore.strand, assign.by, .parallel,
                     getReadsFromSequence, ...)
})

setMethod("countTagsByRegion", c(x="SeqStore"),
function(x, regions, ignore.strand, .parallel, ...) {
  assign.by <- match.arg(assign.by, c('unique', 'all'))
  .countTagsByRegion(x, regions, ignore.strand, assign.by, .parallel,
                     getReadsOnChromosome, ...)
})

## Inner engine for countTagsByRegion generic
.countTagsByRegion <- function(from, regions, ignore.strand,
                               assign.by, .parallel, fetchReads, ...) {
  stopifnot(inherits(regions, "GenomicRanges"))
  stopifnot(is.function(fetchReads))
  assign.by <- match.arg(assign.by, c('unique', 'all'))

  if (.parallel && length(chrs) > 1L) {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }

  chrs <- intersect(seqlevels(from), seqlevels(regions))
  ans <- integer(length(regions))

  counts <- foreach(chr=chrs, .packages=c("GenomicRanges"),
                    .options.multicore=list(preschedule=FALSE)) %loop% {
    reads <- fetchReads(from, chr, ...)
    if (ignore.strand) {
      strand(reads) <- '*'
    }
    o <- assignUniqueOverlaps(reads, regions, .subject.overlap.check=TRUE)
    o.counts <- tabulate(o)
    idx <- which(o.counts != 0)
    list(idx=idx, count=o.counts[idx])
  }

  for (count in counts) {
    if (length(count$idx) > 0L) {
      ans[count$idx] <- count$count
    }
  }

  ans
}

filter.unique.bwa <- function(x) {
  x[mcols(x)$X0 == 1]
}

filter.mapq <- function(min.score) {
  function(x) x[mcols(x)$mapq >= min.score]
}

unstrandedCountTagsByRegion <- function(from, regions, assign.by,
                                        .parallel=TRUE, tags='X0',
                                        filter.fn=NULL) {
  stop("Use GenomicCache::tabulateReads")
  if (is.character(from)) {
    from <- BamFile(from)
  }
  stopifnot(inherits(from, "BamFile"))
  stopifnot(file.exists(path(from)))
  stopifnot(file.exists(paste(path(from), 'bai', sep=".")))

  stopifnot(inherits(regions, "GenomicRanges"))

  if (!is.character(tags)) tags <- character()

  assign.by <- match.arg(assign.by, c('unique-quantify', 'unique-fix', 'all'))

  chrs <- intersect(seqlevels(from), seqlevels(regions))

  if (.parallel && length(chrs) > 1L) {
    '%loop%' <- getFunction('%dopar%')
  } else {
    '%loop%' <- getFunction('%do%')
  }


  values(regions)$.idx. <- 1:length(regions)
  fwd <- regions[strand(regions) == '+']
  rev <- regions[strand(regions) == '-']
  ans <- integer(length(regions))
  si <- seqinfo(from)

  str.regions <- list(fwd=fwd, rev=rev)

  counts <- foreach(chr=chrs, .packages=c("GenomicRanges", "Rsamtools"),
                    .options.multicore=list(preschedule=FALSE)) %loop% {
    param <- ScanBamParam(what='mapq',
                          which=GRanges(chr, IRanges(1, seqlengths(si)[chr])),
                          tag=tags, flag=scanBamFlag(isUnmappedQuery=FALSE))
    reads <- readGappedAlignments(path(from), param=param)
    if (is.function(filter.fn)) {
      reads <- filter.fn(reads)
    }
    strand(reads) <- '*'
    greads <- as(reads, "GRanges")

    cnts <- lapply(str.regions, function(these) {
      if (assign.by %in% c('unique.quantify', 'unique.fix')) {
        ass.by <- gsub('unique.', '', assign.by)
        suppressWarnings(o <- assignUniqueOverlaps(greads, these, ass.by, .subject.overlap.check=TRUE))
        tab <- tabulate(o)
        res <- integer(length(these))
        res[1L:length(tab)] <- tab
      } else {
        res <- countOverlaps(these, greads)
      }
      report <- res > 0
      list(idx=mcols(these)$.idx.[report], count=res[report])
    })

    res <- list(idx=unlist(sapply(cnts, '[[', 'idx')),
                count=unlist(sapply(cnts, '[[', 'count')))
  }

  for (count in counts) {
    if (length(count$idx) > 0L) {
      ans[count$idx] <- count$count
    }
  }

  ans
}

################################################################################
## Defunct

##' Counts all tags w/in the annotated tx bounds for a gene.
##' Also counts the number of disjoint read positions within the tx boundary.
##'
##' This is to summarize tag data at the "gene" level, however ill-defined a
##' "Gene" might be.
countTagsInGeneBounds <- function(x, gcache, which.chr=chromosomes(x),
                                  unique.only=FALSE, same.strand=FALSE,
                                  flank.up=0L, flank.down=flank.up,
                                  min.cover=5, order.by='npeaks') {
  .Defunct(msg='Use annotated reads for this')
  if (!is(x, 'RnaSeqStore')) {
    stop("Need an RnaSeqStore here.")
  }

  ## if (!is(gcache, 'GenomicCache')) {
  ##   stop("Need a GenomicCache to continue")
  ## }

  if (is(gcache, 'GenomicCache')) {
    gene.bounds <- transcriptsBy(gcache, 'gene')
  } else {
    gene.bounds <- gcache
  }

  if (!is(gene.bounds, "GRangesList")) {
    stop("Need a GRangesList to continue")
  }

  counts <- lapply(which.chr, function(chr) {
    cat(chr, "...\n")
    reads <- getReadsOnChromosome(x, chr, unique.only=unique.only,
                                  smooth.by=smooth.by, min.cover=min.cover)
    ## Extending flanks of genes in a GRangesList is really time consuming
    ## so manipulating reads
    ## DEBUG: Slicing with strand(reads) == '+' fails if the Rle has only 1 run?
    is.fwd <- as.vector(strand(reads) == '+')
    is.rev <- as.vector(!is.fwd)

    if (sum(is.fwd) > 0L) {
      if (flank.up > 0) {
        start(reads[is.fwd]) <- start(reads[is.fwd]) - flank.up + 1
      }
      if (flank.down > 0) {
        end(reads[is.fwd]) <- end(reads[is.fwd]) + flank.down - 1
      }
    }

    if (sum(is.rev) > 0L) {
      if (flank.up > 0) {
        end(reads[is.rev]) <- end(reads[is.rev]) + flank.up - 1
      }
      if (flank.down > 0) {
        start(reads[is.rev]) <- start(reads[is.rev]) - flank.down + 1
      }
    }

    if (!same.strand) {
      strand(reads) <- '*'
    }

    count <- t(countOverlaps(gene.bounds, reads))
    reads <- reduce(reads)
    npeaks <- t(countOverlaps(gene.bounds, reads))

    list(count=count, npeaks=npeaks)
  })

  npeaks <- do.call(rbind, lapply(counts, '[[', 'npeaks'))
  counts <- do.call(rbind, lapply(counts, '[[', 'count'))

  ## threeUTRsByTranscript don't ahve $symbol
  if ("symbol" %in% colnames(values(gene.bounds))) {
    symbols <- values(gene.bounds)$symbol
  } else {
    symbols <- names(gene.bounds)
  }

  df <- data.frame(symbol=symbols,
                   count=colSums(counts),
                   npeaks=colSums(npeaks),
                   stringsAsFactors=FALSE)

  if (!is.null(order.by) && order.by %in% colnames(df)) {
    df <- df[order(df[[order.by]], decreasing=TRUE),]
  }

  df
}
