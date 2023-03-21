##' Creates a tag-count matrix from a list of experiments (or CompressedReads
##' generated from experiments) that is suitable to pass into things like edgeR
##' or DESeq
##'
##' @export
##' @seealso \code{\link{tagSummaryByGene}, \link{tagTable}}
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param expts list of \code{\linkS4class{TagSeq}} or \code{CompressedRead}
##' objects representing the data for each experiment.
##' @param by \code{"gene"} or \code{"tag"} indicating how summary counts
##' should be generated.
##' @param as.tpm A \code{logical} indicating if summary statistics should be
##' calculated as raw counts (\code{FALSE}) or tags per million \code{TRUE}.
##' @param lib.sizes A (named) integer vector as long as \code{expts} indicating
##' the library size of each experiment
##' @param ... Elements used to fetch annotated reads if \code{expts} is a list
##' of \code{TagSeq} objects. At a minimum, a \code{\linkS4class{GenomicCache}}
##' object is required to fetch the "default" annotated reads.
##'
##' @return A matrix of gene or tag counts. Each column of the matrix are the
##' counts for the corresponding experiment in \code{expts}. Each row of the
##' matrix is a gene/tag. The \code{rownames} of the matrix give the id of each
##' element.
createCountMatrix <- function(expts, by=c('gene', 'tag'), as.tpm=FALSE,
                              lib.size=NULL, ...) {
  ## Parameter beuracracy
  by <- match.arg(by)
  if (all(sapply(expts, inherits, 'TagSeq'))) {
    stopifnot(inherits(gcache, 'GenomicCache'))
    expts <- lapply(expts, getAnnotatedReads)
  }
  if (!all(sapply(expts, is, 'CompressedReads'))) {
    stop("Neead annotated/compressed reads")
  }
  if (is.null(lib.size)) {
    lib.size <- unlist(lapply(expts, function(x) sum(values(x)$count)))
  }

  ## Select appropriate summary statistic columns and function used to generate
  ## the tag table
  stat.f <- switch(by, gene=tagSummaryByGene, tag=tagTable)
  key.column <- switch(by, gene='entrez.id', tag='id')
  summary.column <- if (as.tpm) 'tpm' else 'count'

  ## Create the stats table over each experiment
  expt.stats <- lapply(1:length(expts), function(i) {
    e.stats <- stat.f(expts[[i]], lib.size=lib.size[[i]])
    e.stats <- e.stats[, c(key.column, summary.column), with=FALSE]
    setkeyv(e.stats, key.column)
    e.stats
  })
  names(expt.stats) <- names(expts)

  ## Specialized merge function to collapse experiment summaries using Reduce
  mergef <- function(x, y) {
    suffix1 <- ncol(x) - 1
    suffix2 <- suffix1 + 1
    suffixes <- as.character(c(suffix1, suffix2))
    suffixes <- paste(".", suffixes, sep="")
    merge(x, y, by=key.column, all=TRUE, suffixes=suffixes)
  }

  M <- as.data.frame(Reduce(mergef, expt.stats))
  row.names <- as.character(M[, 1])
  M <- as.matrix(M[, -1])
  colnames(M) <- names(expts)
  rownames(M) <- row.names
  M
}

##' Returns a data.frame with summary information for each gene.
##'
##' The number of tags over the gene + tpm (expression) is reported as well as
##' the distribution of tag "clusters" over different annotated regions of
##' the gene.
##'
##' @exportMethod tagSummaryByGene
##' @rdname tagSummaryByGene-methods
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param x A TagSeq experiment, or the full CompressedReads object like what
##' is returned from \code{\link{getAnnotatedReads}(x, ...)}
##' @param min.count The threshold to use to filter out reads
##'
##' @return A data.frame with entrez.id, symbol, gene region pile counts,
##' tag counts, and tpm columns for each gene (row)
setGeneric("tagSummaryByGene",
function(x, min.count=1L, lib.size=getReadCount(x), ...) {
  standardGeneric("tagSummaryByGene")
})

##' @nord
setMethod("tagSummaryByGene", c(x="TagSeq"),
function(x, min.count, lib.size, ...) {
  tagSummaryByGene(getAnnotatedReads(x, ...), min.count, lib.size, ...)
})

isAnnotatedGRanges <- function(x) {
  if (!inherits(x, 'GRanges')) return(FALSE)
  required <- c('count', 'entrez.id', 'symbol', 'exon.anno')
  if (!all(required %in% colnames(values(x)))) return(FALSE)
  TRUE
}

##' @nord
setMethod("tagSummaryByGene", c(x="GRanges"),
function(x, min.count, lib.size, ...) {
  stopifnot(isAnnotatedGRanges(x))
  if (missing(lib.size) || is.null(lib.size)) {
    lib.size <- sum(values(x)$count)
  }

  if (min.count > 1L) {
    x <- x[values(x)$count > min.count]
  }

  df <- subset(as.data.frame(x), !is.na(entrez.id))
  colnames(df) <- gsub('seqnames', 'chr', colnames(df))
  for (name in colnames(df)) {
    if (is.factor(df[[name]])) {
      df[[name]] <- as.character(df[[name]])
    }
  }
  mdf <- melt(df, id.vars=c("entrez.id", "symbol", "chr", "strand", "start",
                     "exon.anno"),
               measure.vars=c('count'))

  ## TODO: Why are these dcast's throwing warning messages for
  ##   In FUN(X[[2L]], ...) : NAs introduced by coercion
  ## Including + strand makes duplicate entrez.id combos for genes with
  ## antisense reads, since these have opposite strand
  genic.distro <- dcast(mdf, symbol + entrez.id + chr ~ exon.anno, length)
  genic.counts <- dcast(mdf, symbol + entrez.id + chr ~ exon.anno, sum)

  ## Remove antisense reads for the rest of the summary statistics
  mdf <- subset(mdf, exon.anno != 'antisense')

  ## Counts the number of reads over the gene (used for expression)
  gene.counts <- dcast(mdf, symbol + entrez.id  + chr ~ ., sum)
  colnames(gene.counts)[ncol(gene.counts)] <- 'total.tags'

  ## TODO: You should check that all(this[,4] == 1)
  ## strands <- dcast(subset(mdf, exon.anno != 'antisense'),
  ##                  symbol + entrez.id + chr ~ .,
  ##                  function(x) length(unique(x)),
  ##                  value_var='strand')

  ## Identify which strand the gene lies on, on a chromosome by chromosome
  ## basis
  ## strands <- dcast(subset(mdf, exon.anno != 'antisense'),
  strands <- dcast(mdf,
                   symbol + entrez.id + chr ~ ., function(x) x[1],
                   value_var='strand')
  colnames(strands)[ncol(strands)] <- 'strand'

  m <- merge(gene.counts, strands, by=c('symbol', 'entrez.id', 'chr'))
  m <- merge(m, genic.counts, by=c('symbol', 'entrez.id', 'chr'))
  m <- merge(m, genic.distro, by=c('symbol', 'entrez.id', 'chr'),
             suffixes=c("", ".peaks"))

  transform(m, tpm=TPM(total.tags, lib.size))
})

##' Creates a data.table with a unique identifier for a tag, and its counts.
##'
##' This prepares the data to be used for differential expression on a tag
##' by tag basis.
##'
##' @exportMethod tagTable
##' @rdname tagTable-methods
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param x A \code{\linkS4class{TagSeq}} or
##' \code{\linkS4class{CompressedReads}} object
##' @param ... parameters to pass down to \code{\link{getAnnotatedReads}} if
##' \code{x} is a \code{\linkS4class{TagSeq}} object
##'
##' @return A data.table with id and count columns
setGeneric("tagTable", function(x, ...) standardGeneric("tagTable"))
setMethod("tagTable", c(x="TagSeq"),
function(x, smooth.by='repackTags', ...) {
  tags <- lapply(chromosomes(x), function(chr) {
    as(getReadsOnChromosome(x, chr, smooth.by=smooth.by), 'CompressedReads')
  })
  tags <- do.call(c, tags)
  tagTalbe(tags, ...)
})

.resampleTags <- function(granges, proportion.resample, n.times=1L) {
  .compress <- FALSE
  if (is(granges, 'CompressedReads')) {
    granges <- rep(granges, values(granges)$count)
    .compress <- TRUE
  }
  if (proportion.resample > 1) {
    proportion.resample <- proportion.resample / 100
  }
  n <- floor(length(granges) * proportion.resample)
  re <- lapply(1:n.times, function(i) {
    gr <- granges[sample(1:length(granges), n, replace=TRUE)]
    if (.compress) {
      gr <- compressReads(gr)
    }
    gr
  })
  re
}

##' @nord
setMethod("tagTable", c(x="CompressedReads"),
function(x, lib.size=NULL, ...) {
  if (is.null(lib.size)) {
    lib.size <- sum(values(x)$count)
  }
  ids <- generateReadIdByPosition(x)
  tag.table <- data.table(id=ids, count=values(x)$count)
  tag.table$tpm <- TPM(tag.table$count, lib.size)
  tag.table
})

##' This does resampling -- you probably don't really want to use this
##' @nord
tagTable.resampling <- function(x, resample.times=0L, proportion.resample=.75,
                                expt.prefix='count', ...) {
  tag.tables <- list(x)
  if (resample.times > 0) {
    tt <- .resampleTags(x, proportion.resample, resample.times)
    tag.tables <- c(tag.tables, tt)
  }

  tag.tables <- lapply(tag.tables, function(tt) {
    ids <- generateReadIdByPosition(tt)
    tag.table <- data.table(id=ids, values(tt)$count)
    setkeyv(tag.table, 'id')
    tag.table
  })

  result <- tag.tables[[1]]
  if (length(tag.tables) > 1) {
    for (i in 2:length(tag.tables)) {
      result <- merge(result, tag.tables[[i]], by='id', all=TRUE)
    }
  }

  result <- as.data.frame(result)
  result[is.na(result)] <- 0L
  result <- data.table(result, key='id')
  colnames(result)[2] <- expt.prefix
  if (resample.times > 0L) {
    colnames(result)[(1:resample.times)+2] <- paste(expt.prefix, 'FR', 1:resample.times)
  }

  result
}

##' @nord
generateReadIdByPosition <- function(granges) {
  paste(seqnames(granges), start(granges), end(granges), strand(granges), sep='.')
}

##' @nord
tagIdToGRanges <- function(tag.id) {
  info <- strsplit(as.character(tag.id), '.', fixed=TRUE)
  seqs <- sapply(info, '[', 1)
  starts <- as.integer(sapply(info, '[', 2))
  ends <- as.integer(sapply(info, '[', 3))
  strands <- sapply(info, '[', 4)

  gr <- GRanges(seqs, IRanges(starts, ends), strands)
  gr
}
