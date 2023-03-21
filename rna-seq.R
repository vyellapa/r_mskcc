##' Tabules the reads over an AnnotatedGenome object.
##'
##' @param x AnnotatedGenome object
##' @param reads Likely a BamFile
##' @param ignore.strand Consider strand of reads significant. If FALSE,
##' the overlap of reads to regions is done one strand of a time -- this means
##' that one read can be counted for two features that are on different strands.
##' @param .parallel Parallelize this over chromosomes
##' @param scan.bam.what The "what" to include in the scaBam query. Defaults
##' to \code{mapq} (in addition to the other info returned in GappedAlignments)
##' @param scan.bam.tags Tags to include from the BAM file.
##' @param scan.bam.flag An integer flag (use \code{scanBamFlag}) to use
##' in the inner \code{scanBam} query.
##' @param filter.fn A function that accepts a GappedAlignments object
##' and returns the GappedAlignments to keep/tabulate. Typically you modify
##' what to include with \code{scan.bam.what} and \code{scan.bam.tags} and
##' use those to filter
##' @return An integer vector of counts for each region in \code{x}
tabulateReads <- function(x, from, assign.by='unique-quantify',
                          ignore.strand=FALSE,
                          scan.bam.param=bwaSamseParams(),
                          filter.fn=NULL, chrs=NULL,
                          min.cover=1L, ..., .parallel=TRUE) {
  if (is.character(from)) {
    from <- BamFile(from)
  }
  stopifnot(inherits(from, "BamFile"))
  stopifnot(file.exists(path(from)))
  stopifnot(file.exists(paste(path(from), 'bai', sep=".")))
  stopifnot(is(scan.bam.param, "ScanBamParam"))
  stopifnot(inherits(x, "GenomicRanges"))

  assign.by <- match.arg(assign.by, c('unique-quantify', 'unique-fix', 'all'))

  all.chrs <- intersect(seqlevels(from), seqlevels(x))
  if (is.null(chrs)) {
    chrs <- all.chrs
  }
  stopifnot(all(chrs %in% all.chrs))

  values(x)$.idx. <- 1:length(x)
  regions <- sapply(c("+", "-"), function(strnd) {
    x[strand(x) == strnd]
  }, simplify=FALSE)

  strands <- names(regions)
  si <- seqinfo(from)
  ans <- data.frame(count=integer(length(x)), percent=numeric(length(x)))

  if (any(mcols(x)$exon.anno %in% "intron") && assign.by != 'all') {
    warning("Reads that span introns do not work well with quantifyOverlaps",
            immediate.=TRUE)
  }

  fn <- if (.parallel && length(chrs) > 1L) "%dopar%" else "%do%"
  "%loop%" <- getFunction(fn)
  pkgs <- c("GenomicRanges", "Rsamtools")
  opts <- list(preschedule=FALSE)

  counts <- foreach(chr=chrs, .packages=pkgs, .options.multicore=opts) %loop% {
    param <- scan.bam.param
    bamWhich(param) <- GRanges(chr, IRanges(1, seqlengths(si)[chr]))
    reads <- readGAlignments(path(from), param=param)
    if (is.function(filter.fn)) {
      reads <- filter.fn(reads)
    }

    if (ignore.strand) {
      strand(reads) <- '*'
      cover <- coverage(reads)[[chr]]
      cover <- list("+"=cover, "-"=cover)
    } else {
      cover <- sapply(c("+", "-"), function(s) {
        coverage(reads[strand(reads) == s])[[chr]]
      }, simplify=FALSE)
    }

    greads <- as(reads, "GRanges")

    cnts <- sapply(strands, function(s) {
      these <- regions[[s]]
      if (assign.by %in% c('unique.quantify', 'unique.fix')) {
        ass.by <- gsub('unique.', '', assign.by)
        suppressWarnings({
          o <- assignUniqueOverlaps(greads, these, ass.by,
                                    .subject.overlap.check=TRUE)
        })
        tab <- tabulate(o)
        res <- integer(length(these))
        res[1L:length(tab)] <- tab
      } else {
        res <- countOverlaps(these, greads)
      }

      keep <- res > 0
      rc <- these[keep]
      cnt <- res[keep]

      ## Calculate the percent of each region that is covered by a read
      v <- Views(cover[[s]] >= min.cover, ranges(rc))
      list(idx=mcols(rc)$.idx., count=cnt, percent=viewSums(v) / width(v))
    }, simplify=FALSE)

    res <- list(idx=unlist(sapply(cnts, '[[', 'idx')),
                count=unlist(sapply(cnts, '[[', 'count')),
                percent=unlist(sapply(cnts, '[[', 'percent')))
  }

  for (count in counts) {
    if (length(count$idx)) {
      ans$count[count$idx] <- count$count
      ans$percent[count$idx] <- count$percent
    }
  }

  ans
}

## Parallelization is over each chromosome, so experiments are processed
## serially
tabulateIntoSummarizedExperiment <- function(x, input, colData=DataFrame(),
                                             exptData=SimpleList(),
                                             ..., .parallel=TRUE) {
  if (is.character(input) || is.character(unlist(input))) {
    input <- lapply(input, BamFile)
  }
  stopifnot(is.character(names(input)))
  stopifnot(all(sapply(input, is, "BamFile")))
  stopifnot(all(file.exists(sapply(input, path))))
  stopifnot(all(file.exists(sapply(input, function(f) paste(path(f), "bai", sep=".")))))

  info <- lapply(input, function(bf) tabulateReads(x, bf, ..., .parallel=.parallel))
  tabulated <- SimpleList(count=sapply(info, '[[', 'count'),
                          percent=sapply(info, '[[', 'percent'))
  SummarizedExperiment(tabulated, rowData=x, colData=colData, exptData=exptData)
}


##' Calculate RPKM from a SummarizedExperiment generated from
##' \code{tabulateIntoSummarizedExperiment}
##'
##' @param x A SummarizedExperiment generated from
##' \code{tabulateIntoSummarizedExperiment}
##' @param keys The \code{mcols} to use as keys -- indicate transcript units
##' @param min.count Minimum number of observed reads required for the exon
##' (across the atlas) to be considered for inclusion in the RPKM calc
##' @param with.percent Use \code{percent} column to calculate the "real"
##' length of the transcript (the K in RPKM).
##' (TODO: implement calcRPKM,with.percent=TRUE)
calcRPKM <- function(x, keys=c('seqnames', 'strand', 'entrez.id'),
                     min.count=1L, with.percent=FALSE, as.log=TRUE) {
  if (with.percent) {
    stop("with.percent not yet implemented")
  }
  stopifnot(is(x, "SummarizedExperiment"))
  stopifnot(all(setdiff(keys, c('seqnames', 'strand', 'start', 'end')) %in% names(mcols(x))))
  if (is.null(colnames(x))) {
    warning("No experiment names -- adding faux names", immediate.=TRUE)
    colnames(x) <- paste("expt", 1:ncol(x), sep=".")
  }

  if (is.numeric(min.count)) {
    x <- x[rowSums(assay(x)) >= min.count]
  }

  f <- cbind(as(rowData(x), 'data.table'), as.data.table(assay(x)))
  expt.cols <- colnames(x)
  setkeyv(f, keys)

  x.expr <- f[, {
    expr <- sapply(expt.cols, function(e) sum(.SD[[e]]), simplify=FALSE)
    c(list(symbol=symbol[1], len=sum(width)), expr)
  }, by=list(seqnames, strand, entrez.id)]

  for (e in expt.cols) {
    rpk <- paste(e, 'rpk', sep='.')
    rpkm <- paste(e, 'rpkm', sep='.')
    x.expr[[rpk]] <- 1e3 * (x.expr[[e]] / x.expr[['len']])
    x.expr[[rpkm]] <- 1e6 * (x.expr[[rpk]] / sum(x.expr[[e]]))
    if (as.log) {
      x.expr[[rpk]] <- log2(x.expr[[rpk]])
      x.expr[[rpkm]] <- log2(x.expr[[rpkm]])
    }
  }

  x.expr
}

