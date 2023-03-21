.noReads <- function() {
  reads <- GRanges(seqnames=character(), ranges=IRanges(), strand=strand())
  values(reads) <- DataFrame(id=integer(), pair.id=integer())
  reads
}

## Assume `reads` is only from one chromosome
## reads that begin at same start position are assumed to be the same
splitMultimap <- function(seqstore, reads, tag='X0', max.split=10L, ...) {
  warning("Not splitting multimapped reads -- do the read/rescue")
  return(reads)

  stopifnot(inherits(reads, 'GRanges'))
  stopifnot(tag %in% colnames(values(reads)))
  stopifnot(length(unique(seqnames(reads))) == 1L)

  do.multimap <- function(dt) {
    if (nrow(dt) == 0L) {
      return(dt)
    }
    strands <- unique(dt$strand)
    stopifnot(length(strands) == 1L)
    col.names <- colnames(dt)

    ## Setting the key as either start,end is important for a SAGEseq
    ## experiment that didn't run through a repackTags smoothing. If it
    ## (or a TailSeq) experiment was smoothed thusly, then either start/end
    ## wouldn't make a difference
    pk <- if (strands == '-') 'end' else 'start'
    setkeyv(dt, c('seqnames', pk))
    bysub <- parse(text=paste("list(", pk, ")", sep=""))[[1]] ## UGH
    dt <- dt[, {
      n <- max(0L, floor(nrow(.SD) / max(.SD[[tag]])))
      if (n == 0L) {
        result <- .SD[1,]
        result$.keep <- FALSE
      } else {
        result <- .SD[sample(nrow(.SD), n), ]
        result$.keep <- rep(TRUE, n)
      }
      result
    }, by=eval(bysub)]
    dt <- dt[.keep]
    dt[, col.names, with=FALSE]
  }

  reads <- reads[values(reads)[[tag]] <= max.split]
  x <- as(reads, 'data.table')
  is.rev <- x$strand == '-'
  multimapped <- x$X0 > 1L

  x.fwd <- do.multimap(subset(x, multimapped & !is.rev))
  x.rev <- do.multimap(subset(x, multimapped & is.rev))

  result <- rbind(subset(x, !multimapped), x.fwd, x.rev)
  result <- as(result, 'GRanges')
  result[order(ranges(result))]
}

################################################################################
# These functions are no longer used
# Rely on SeqTools::getReadsFromSequence for this functionality. Even better
# is to use the "fetching reads" functionality from GenomicRanges and Rsamtools
# ------------------------------------------------------------------------------

# ##' Returns reads on this chromosome from a TailSeq experiment.
# ##'
# ##' This has to be treated slightly different than other tag/rna-seq experiments
# ##' since the protocol produces "flipped" reads. This method includes a
# ##' \code{reverse.strand} parameter to deal with this accordingly.
# ##'
# ##' The 3' TTTs of the transcript are on the 5' side of the reads as the come
# ##' off the sequencer.
# ##'
# ##' @param reverse.strand If \code{TRUE}, the strand of the reads are "flipped"
# ##' so that they correctly indicate which strand the original read came from
# ##' (and not the direction it was read off of the sequencer).
# setMethod("getReadsOnChromosome", c(x="TagSeq"),
# function(x, chromosome, start, end, strand, unique.only, max.mismatch,
#          split.multimap=FALSE, smooth.by=NULL, meta.what, opposite.strand,
#          tag=c('NM', 'MD'), ...) {
#   trace <- checkTrace(...)
#   if (trace) cat("getReadsOnChromosome::TagSeq\n")
#
#   meta <- metadata(x, ...)
#   if (!is.null(meta$post.processed) && meta$post.processed) {
#     if (missing(unique.only)) {
#       unique.only <- FALSE
#     }
#     if (missing(max.mismatch)) {
#       max.mismatch <- NULL
#     }
#     tag <- unique(c('Z0', 'ZZ', tag))
#   }
#
#   ## Gaurd against passing split.multimap=TRUE to upstream SeqStore class
#   ## since I don't know how to handle this similarly for the non-tagseq case.
#   ## split.mm <- split.multimap
#   ## split.multimap <- FALSE
#
#   reads <- callNextMethod(x, chromosome, start, end, strand, unique.only,
#                           max.mismatch, split.multimap, smooth.by, meta.what,
#                           opposite.strand=opposite.strand, tag=tag, ...)
#
#   ## if (!unique.only && split.mm) {
#   ##   reads <- splitMultimap(x, reads, tag='tag.X0', ...)
#   ## }
#
#   reads
# })
#
# setMethod("getReadsOnChromosome", c(x="TailSeq"),
# function(x, chromosome, start, end, strand, unique.only, max.mismatch,
#          split.multimap, smooth.by=NULL, meta.what, opposite.strand,
#          tag=c('NM', 'MD'), ...) {
#   trace <- checkTrace(...)
#   if (trace) cat("getReadsOnChromosome::TailSeq\n")
#
#   reads <- callNextMethod(x, chromosome, start, end, strand, unique.only,
#                           max.mismatch, split.multimap, smooth.by, meta.what,
#                           opposite.strand=opposite.strand, tag=tag, ...)
#
#   reads
# })
#
# setMethod("getReadsOnChromosome", c(x="SAGEseq"),
# function(x, chromosome, start, end, strand, unique.only, max.mismatch,
#          split.multimap,
#          ## smooth.by=c('anchorSite', 'repackTags'),
#          smooth.by=NULL, ## no mismatch in restriction specified during alignment
#          meta.what, opposite.strand=TRUE, ...) {
#   trace <- checkTrace(...)
#   if (trace) cat("getReadsOnChromosome::TailSeq\n")
#   chromosome <- as.character(chromosome)
#
#   if (any(c('anchorSite', 'smoother.anchorSite') %in% smooth.by)) {
#     bs.chr <- unmasked(getBsGenome(x)[[chromosome]])
#   } else {
#     bs.chr <- NULL
#   }
#
#   reads <- callNextMethod(x, chromosome, start, end, strand, unique.only,
#                           max.mismatch, split.multimap, smooth.by, meta.what,
#                           bs.chr=bs.chr, anchor.site=anchorSite(x),
#                           opposite.strand=opposite.strand, ...)
#   reads
# })
#
# setMethod("getReadsOnGene", c(x="TagSeq", gene="GFGene"),
# function(x, gene, flank.up, flank.down, same.strand, unique.only, max.mismatch,
#          split.multimap, smooth.by=c('repackTags'), meta.what, which.chr, ...) {
#   trace <- checkTrace(...)
#   if (trace) cat("getReadsOnGene::TagSeq\n")
#   callNextMethod(x, gene, flank.up, flank.down, same.strand, unique.only,
#                  max.mismatch, split.multimap, smooth.by, meta.what, which.chr,
#                  ...)
# })
#
# setMethod("getReadsOnGene", c(x="SAGEseq", gene="GFGene"),
# function(x, gene, flank.up, flank.down, same.strand, unique.only,
#          max.mismatch, split.multimap, smooth.by=c('anchorSite', 'repackTags'),
#          meta.what, which.chr, ...) {
#   trace <- checkTrace(...)
#   if (trace) cat("getReadsOnGene::SAGEseq\n")
#   callNextMethod(x, gene, flank.up, flank.down, same.strand, unique.only,
#                  max.mismatch, split.multimap, smooth.by, meta.what, which.chr,
#                  ...)
# })
#
# setMethod("getReadCount", c(x="CompressedReads"),
# function(x, strand, unique.only, split.multimap, smooth.by, ...) {
#   if (!is.null(strand)) {
#     x <- x[strand(x) == strand]
#   }
#   sum(values(x)$count)
# })


