## Centralized place to keep filename logic. When writing to or reading from
## a file on disk, call appropriate methods here instead of hard wiring the
## names of files in the code
##
## Functions should be vectorized to return results for several file queries
## at a time.

idealizedGeneFN <- function(annotation.source, gene.collapse='longest',
                            cds.cover='min', flank.up=1000L, flank.down=1000L,
                            ...) {
  message('cds.cover parameter not used yet')
  if (inherits(annotation.source, 'GenomicCache')) {
    annotation.source <- annotationSource(annotation.source)
  }
  stopifnot(is(annotation.source, 'character'))
  gene.collapse <- GenomicFeaturesX:::matchGFGeneCollapse(gene.collapse)

  fn <- paste(annotation.source, sprintf('collapse-%s', gene.collapse),
              sprintf('cdscover-%s', cds.cover),
              sprintf("up-%d.down-%d", flank.up, flank.down), sep=".")
  fn
}

## The way reads are returned determine most of the filenaming conventions
## of the files stored. The user can work on uniquely mapping reads, reads
## with some mismatches, etc.
readsCombinationFN <- function(unique.only=TRUE, split.multimap=FALSE,
                               max.mismatch=0L, max.split=10L, ...) {
  fn <- sprintf('reads.unique-%s', unique.only)
  if (!unique.only && split.multimap) {
    fn <- sprintf('%s.split-multimap-%d', fn, max.split)
  }
  if (is.null(max.mismatch)) {
    max.mismatch <- 'NULL'
  } else {
    max.mismatch <- as.character(max.mismatch)
  }
  fn <- sprintf('%s.maxmismatch-%s', fn, max.mismatch)
  fn
}

signalFN <- function(path, chr.name, kernel='normal', bandwidth=18,
                     unique.only=TRUE, max.mismatch=0L, ...) {
  fn.reads <- readsCombinationFN(unique.only, max.mismatch=max.mismatch, ...)
  if (inherits(path, 'SeqStore')) {
    path <- file.path(metaDirectory(path, ...), 'annotated.signal')
  }
  fn <- paste('annotated-signal', kernel, sprintf('bandwidth-%d', bandwidth),
              fn.reads, sep=".")
  if (!missing(chr.name) && !is.null(chr.name) &&
      all(sapply(chr.name, nchar) > 0)) {
    fn <- paste(fn, chr.name, sep=".")
  }
  fn <- paste(fn, 'rda', sep='.')
  file.path(path, fn)
}

annotatedReadsFN <- function(path, annotation.source, chr.name,
                             gene.collapse='longest', cds.cover='min',
                             flank.up=1000L, flank.down=flank.up,
                             unique.only=TRUE, split.multimap=FALSE,
                             max.mismatch=0L, stranded=TRUE,
                             nuke.antisense=TRUE, ...) {
  fn.gene <- idealizedGeneFN(annotation.source, gene.collapse, cds.cover,
                             flank.up, flank.down)
  fn.reads <- readsCombinationFN(unique.only, split.multimap, max.mismatch, ...)
  if (inherits(path, 'SeqStore')) {
    path <- file.path(metaDirectory(path), 'annotated.reads')
  }

  fn <- paste('annotated-reads', fn.gene, fn.reads, sep=".")
  if (!missing(chr.name) && !is.null(chr.name) &&
      all(sapply(chr.name, nchar) > 0)) {
    fn <- paste(fn, chr.name, sep=".")
  }
  fn <- paste(fn, 'rda', sep='.')
  file.path(path, fn)
}

geneSummaryTableFN <- function(path, annotation.source, gene.collapse='longest',
                               cds.cover='min', flank.up=1000L,
                               flank.down=flank.up, unique.only=TRUE,
                               split.multimap=FALSE, max.mismatch=0L, min.count=1,
                               dominant.killer=0.1, stranded=TRUE, ...) {
  message('cds.cover parameter not used yet')
  if (inherits(path, 'SeqStore')) {
    path <- file.path(metaDirectory(path))
  }
  fn.gene <- idealizedGeneFN(annotation.source, gene.collapse, cds.cover,
                             flank.up, flank.down)
  fn.reads <- readsCombinationFN(unique.only, split.multimap, max.mismatch, ...)

  fn <- paste('gene-summary-table', sprintf('mincount-%d', min.count),
              sprintf("dominant-%.2f", dominant.killer),
              fn.gene, fn.reads, "rda", sep=".")
  file.path(path, fn)
}

utr3StatsFN <- function(path, annotation.source, chr.name,
                        gene.collapse='longest', cds.cover='min',
                        flank.up=1000L, flank.down=1000L,
                        unique.only=TRUE, split.multimap=FALSE,
                        max.mismatch=0, ...) {
  fn.gene <- idealizedGeneFN(annotation.source, gene.collapse, cds.cover,
                             flank.up, flank.down)
  fn.reads <- readsCombinationFN(unique.only, split.multimap, max.mismatch, ...)
  fn <- paste('utr3stats', fn.gene, fn.reads, sep=".")
  if (!missing(chr.name) && !is.null(chr.name) &&
      all(sapply(chr.name, nchar) > 0)) {
    fn <- paste(fn, chr.name, sep=".")
  }
  fn <- paste(fn, 'rda', sep='.')
  file.path(path, fn)
}

apaSummaryTableFN <- function(path, annotation.source, gene.collapse='longest',
                              cds.cover='min', flank.up=1000L,
                              flank.down=flank.up, unique.only=TRUE,
                              split.multimap=FALSE, max.mismatch=0L, min.count=1L,
                              ...) {
  if (inherits(path, 'SeqStore')) {
    path <- file.path(metaDirectory(path))
  }
  fn.gene <- idealizedGeneFN(annotation.source, gene.collapse, cds.cover,
                             flank.up, flank.down)
  fn.reads <- readsCombinationFN(unique.only, split.multimap, max.mismatch, ...)
  min.count <- as.integer(min.count)
  fn <- paste('apa-summary-table',  sprintf('mincount-%d', min.count),
              fn.gene, fn.reads, 'rda', sep='.')
  file.path(path, fn)
}
