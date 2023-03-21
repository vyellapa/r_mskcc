## New post processing
##
## == Removing Antisense Reads ==
## 1. Look at sense/antisense ratio in row data over non overlapping gene
##    regions.
##    a. fencepost regions into ones w/ *any* reads -- no threshold
##    b. find ratio of sense/antisense in these regions.
## 2. Remove any events that have overlapping sense/antisense reads
##    (not just ones in annotated regions).
##    This can lose info over regions that have genes on opposite strands
##    Is there some "shift" pattern that we see in artifact vs. normal?
##
## == Remove Internal Priming ==
## ...

generateUniqueLawn <- function(x, opposite.strand=FALSE) {
  if (inherits(x, 'SeqStore')) {
    x <- BamFile(bamFile(x))
  }
  if (is.character(x)) {
    x <- BamFile(x)
  }
  stopifnot(inherits(x, 'BamFile'))

  chrs <- as.character(seqnames(x))
  cat("=== Generating Unique Lawn ===\n")
  unique.lawn <- foreach(chr=chrs, .packages="SeqTools") %dopar% {
    r <- getReadsFromSequence(x, chr, unique.only=TRUE, max.mismatch=0L,
                              smooth.by=NULL, opposite.strand=opposite.strand)
    bounds <- lapply(c("+", "-"), function(.strand) {
      .ranges <- ranges(r[strand(r) == .strand])
      covr <- coverage(.ranges)
      sbounds <- slice(covr, lower=1, rangesOnly=TRUE)
      if (length(sbounds) == 0) {
        return(GRanges())
      }
      GRanges(chr, sbounds, .strand, count=countOverlaps(sbounds, .ranges))
    })
    do.call(c, bounds)
  }

  suppressWarnings(do.call(c, unname(unique.lawn)))
}

##' Assign multimapped reads to their "most probable" unique location
##' if one exists.
##'
##' We establish a "unique lawn" of reads where we try to assign multimapped
##' reads back to.
##'
##' Requires that the experiment has been annotated using unique.only=TRUE and
##' max.mismatch=0
##'
##' NOTE: tag.XA is also NA when the number of "best" matches is less than
##' \code{max.mm.consider}, but the number of "non optimal" matches pushes
##' it over the limit.
##'
##' For multimapped reads that are "remapped" back to the unique lawn, I
##' reuse the bwa XT tag and set it to L (for lawn(?))
postProcessAlignments <- function(expt, annotated.genome, max.mismatch=2L,
                                  max.mm.consider=5L, do.annotate=TRUE,
                                  remove.internal.priming=TRUE,
                                  chrs=NULL, nuke.antisense=FALSE,
                                  min.lawn.coverage=iTPM(1, expt),
                                  unique.lawn=NULL, unique.transcribed=TRUE,
                                  kernel='normal', bandwidth=18,
                                  bam.file=NULL, polya.rescue=kPolyA$dna.signal[1:5],
                                  path=directory(expt),
                                  outfile='realign',
                                  opposite.strand=metadata(expt)$opposite.strand,
                                  ...) {
  if (is.null(chrs)) {
    chrs <- as.character(seqnames(expt))
  }
  if (is.null(bam.file)) {
    bam.file <- bamFile(expt)
  }
  if (!inherits(annotated.genome, 'AnnotatedChromosome')) {
    stop("need an annotated.genome")
  }

  header <- scanBamHeader(bam.file)[[1]]

  if (nuke.antisense) {
    do.annotate <- TRUE
  }

  ## Write header
  final.fn <- file.path(path, paste(outfile, 'sam', sep="."))
  df <- data.frame(names(header$text), sapply(header$text, '[[', 1),
                   sapply(header$text, '[[', 2))
  write.table(df, final.fn, append=FALSE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

  ##############################################################################
  ## Generate the unique lawn
  if (is.null(unique.lawn)) {
    unique.lawn <- generateUniqueLawn(bam.file, opposite.strand=opposite.strand)
  }

  if (min.lawn.coverage > 1) {
    unique.lawn <- unique.lawn[values(unique.lawn)$count >= min.lawn.coverage]
  }

  ##############################################################################
  ## Remap the reads back to the lawn
  if (is.null(annotated.genome)) {
    annotated.genome <- getAnnotatedGenome(gcache, ...)
  }

  if (unique.transcribed) {
    is.transcribed <- annotated.genome[values(annotated.genome)$exon.anno != 'intergenic']
    unique.lawn <- subsetByOverlaps(unique.lawn, is.transcribed)
  }

  cat("=== Remapping Reads ===\n")
  remap.chrs <- mclapply(chrs, mc.preschedule=FALSE, function(chr) {
    cat("------", chr, "------\n")
    all.bam.what <- SeqTools:::all.bam.what
    reads <- getReadsFromSequence(BamFile(bam.file), chr, unique.only=FALSE,
                                  max.mismatch=max.mismatch, smooth.by=NULL,
                                  meta.what=all.bam.what, with.sequence=TRUE,
                                  opposite.strand=opposite.strand)
    if (length(reads) == 0L) {
      return(NULL)
    }

    cat('... realigning.\n')

    ## realignReads can return reads that land on a different chromosome
    reads <- realignReads(reads, unique.lawn, max.mm.consider, max.mismatch)
    ## cat("(", chr, " [realigned]): will save these tmp files:",
        ## paste(as.character(unique(seqnames(reads))), collapse=','), "\n")
    ## if (!is.null(ignore)) {
    ##   reads <- reads[!seqnames(reads) %in% ignore]
    ## }

    if (do.annotate) {
      cat('... annotating.\n')
      ## anno.reads <- tryCatch({
      ##   annotateReads(reads, annotated.genome, nuke.antisense=nuke.antisense)
      ## }, error=function(e) NULL)
      reads <- annotateReads(reads, annotated.genome,
                             nuke.antisense=nuke.antisense)
      anno.column <- which(colnames(values(reads)) == 'exon.anno')
      colnames(values(reads))[anno.column] <- 'tag.ZZ'
      ## cat("(", chr, "[annotated]): will save these tmp files:",
          ## paste(as.character(unique(seqnames(reads))), collapse=','), "\n")
    }

    ## Save result temporarily until all reads are remmapped. Each group of reads
    ## will be saved to separate chromosome files.
    ## After processing, all of the reads will be processed for each chromosome
    ## at a time to do "hump finding," then internal-priming axe, etc.
    to.chrs <- as.character(unique(seqnames(reads)))
    cat("(", chr, "): will save these tmp files:", paste(to.chrs, collapse=','), "\n")

    for (schr in to.chrs) {
      fn <- file.path(path, sprintf('%s-from-%s-to-%s.rda', outfile, chr, schr))
      cat("... (", chr, ") saving: ", fn, "\n", sep="")
      sreads <- reads[seqnames(reads) == schr]
      save(sreads, file=fn)
    }

    to.chrs
  })
  remap.chrs <- unique(unlist(remap.chrs))

  cat("=== Reloading Reads for Hump Finding ===\n")
  ## Load realigned reads by chromosome
  ## Do the "hump-finding" thing
  ## (optionally) Axe internal events.
  ## Serialize a new SAM file

  bsg <- getBsGenome(expt)
  sam.files <- mclapply(remap.chrs, mc.preschedule=FALSE, function(chr) {
    cat("...", chr, "...\n")
    rda.files <- list.files(path, sprintf('%s.*to-%s.rda', outfile, chr),
                            full.names=TRUE)
    parts <- lapply(rda.files, function(fn) {
      creads <- load.it(fn)
      meta <- values(creads)
      values(creads) <- NULL
      list(granges=creads, meta=meta)
    })
    reads <- do.call(c, lapply(parts, '[[', 'granges'))
    values(reads) <- combineDataFrames(lapply(parts, '[[', 'meta'))

    signal <- summarizeSignal(expt, kernel, bandwidth, chrs=chr,
                              cached.reads=reads, with.readcount=FALSE,
                              do.save=FALSE, .parallel=FALSE)[[1L]]
    sig.file <- file.path(path, sprintf('smooth-signal.%s.rda', chr))
    save(signal, file=sig.file)

    signal <- as(signal, 'GRanges') ## splits events into "sub-peaks"

    ## from process-internal-priming.R
    values(signal)$is.primed <-
      flagInternallyPrimedPeakAnnotations(signal, bsg, a.count=8L, ip.window=10L,
                                          ip.distal.downstream=0L,
                                          ip.internal.downstream=5L,
                                          flag.by='position')

    ## Attempt to rescue events that look like they are internally primed
    ## by looking for poly(A) signals in them?
    if (!is.null(polya.rescue) && length(polya.rescue) > 0L) {
      ## From process-polya.R
      pa.stats <- collectPolyAPositionStatistics(signal, bsg, polya.rescue,
                                                 .parallel=FALSE)
      pa.good <- isRealPolyA(pa.stats, signal)
      pa.conf <- pa.stats[pa.good,]
      values(signal)$is.real <- flagRealEventsByPolyA(signal, pa.conf)
    } else {
      values(signal)$is.real <- !values(signal)$is.primed
    }

    fn <- file.path(path, sprintf('%s-signal-%s.rda', outfile, chr))
    cat("... saving", fn, "\n")
    save(signal, file=fn)

    ## map2signal <- assignUniqueOverlaps(reads, signal)
    ## na.map <- is.na(map2signal)
    ## if (any(na.map)) {
    ##   message(sum(na.map), "NA's in map from read to signal/event -- why?\n")
    ## }

    ############################################################################
    ## Dropping reads that aren't "real"
    ## TODO: Do we want to flag reads as belonging to internally primed events
    ##       but still serialize them? We can remove them later using an
    ##       appropriate filter rule with bamtools?
    ## real.reads <- reads[!na.map & values(signal)$is.real[map2signal]]
    ## cat("(", chr, "): dropping", length(reads) - length(real.reads), "reads\n")

    sam <- toSamTable(reads)
    fn <- file.path(path, paste(chr, outfile, 'sam', sep="."))
    cat("... saving", nrow(sam), "rows to", fn, "\n")
    write.table(sam, fn, sep="\t", col.names=FALSE, row.names=FALSE,
                quote=FALSE, append=FALSE)
    fn
  })


  ## cat("=== combine files manually and generate BAM file manually ===\n")
  ## cat("for f in chr*sam; do echo $f; cat $f >> realign.sam; done\n")

  real.final.fn <- gsub('\\.sam$', '-final.sam', final.fn)
  cat("\n\n=== combining sam files to", real.final.fn, "===\n")

  system(paste("cp", final.fn, real.final.fn))
  for (sf in sam.files) {
    cmd <- paste("cat", sf, ">>", real.final.fn)
    cat("...", cmd, "\n")
    system(cmd)
  }

  cat("Creating BAM file\n")
  asBam(real.final.fn, gsub("\\.sam$", "", real.final.fn))

  invisible(NULL)
}

realignReads <- function(x, destination, annotated.genome, max.mismatch=2L,
                         max.mm.consider=5L, do.annotate=TRUE, chrs=NULL,
                         nuke.antisense=FALSE, min.lawn.coverage=1,
                         unique.lawn=NULL, unique.transcribed=TRUE,
                         bam.file=NULL, path=directory(expt),
                         outfile='realign-only', opposite.strand, ...) {
  verbose <- checkVerbose(...)
  if (inherits(x, 'TagSeq')) {
    bam.file <- bamFile(x)
  }
  if (is.character(x)) {
    bam.file <- bamFile(x)
  }
  if (file.exists(destination)) {
    stop(destination, " already exists.")
  }
  if (missing(opposite.strand)) {
    opposite.strand <- metadata(expt)$opposite.strand
  }

  out.path <- dirname(destination)
  out.name <- basename(destination)
  out.name <- gsub("\\.bam", "", out.name)
  out.name <- gsub("\\.sam", "", out.name)
  out.file <- file.path(out.path, paste(out.name, 'sam', sep="."))

  ## write.header
  header <- scanBamHeader(x)[[1]]
  df <- data.frame(names(header$text), sapply(header$text, '[[', 1),
                   sapply(header$text, '[[', 2))
  write.table(df, out.file, append=FALSE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

  if (is.null(chrs)) {
    chrs <- as.character(df[[2]])
    chrs <- gsub("SN:", "", chrs[grep("SN:", chrs)])
  }
  if (!inherits(annotated.genome, 'AnnotatedChromosome')) {
    stop("need an annotated.genome")
  }

  #############################################################################
  ## TODO: CONTINUE FROM HERE
  stop("This needs refactoring")
  if (nuke.antisense) {
    do.annotate <- TRUE
  }


  header <- scanBamHeader(bam.file)[[1]]
  df <- data.frame(names(header$text), sapply(header$text, '[[', 1),
                   sapply(header$text, '[[', 2))


  ## Write header
  final.fn <- file.path(path, paste(outfile, 'sam', sep="."))
  write.table(df, final.fn, append=FALSE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

  ##############################################################################
  ## Generate the unique lawn
  if (is.null(unique.lawn)) {
    unique.lawn <- generateUniqueLawn(bam.file, opposite.strand=opposite.strand)
  }

  if (min.lawn.coverage > 1) {
    unique.lawn <- unique.lawn[values(unique.lawn)$count >= min.lawn.coverage]
  }

  ##############################################################################
  ## Remap the reads back to the lawn
  if (is.null(annotated.genome)) {
    annotated.genome <- getAnnotatedGenome(gcache, ...)
  }

  if (unique.transcribed) {
    ignore.lawn <- c('intergenic', 'utr3*', 'utr5*')
    keep.lawn <- !(values(annotated.genome)$exon.anno %in% ignore.lawn)
    is.transcribed <- annotated.genome[keep.lawn]
    unique.lawn <- subsetByOverlaps(unique.lawn, is.transcribed)
  }

  cat("=== Remapping Reads ===\n")
  remap.chrs <- mclapply(chrs, mc.preschedule=FALSE, function(chr) {
    cat("------", chr, "------\n")
    all.bam.what <- SeqTools:::all.bam.what
    reads <- getReadsFromSequence(BamFile(bam.file), chr, unique.only=FALSE,
                                  max.mismatch=max.mismatch, smooth.by=NULL,
                                  meta.what=all.bam.what, with.sequence=TRUE,
                                  opposite.strand=opposite.strand)
    if (length(reads) == 0L) {
      return(NULL)
    }

    cat('... realigning.\n')

    ## realignReads can return reads that land on a different chromosome
    reads <- realignReads(reads, unique.lawn, max.mm.consider, max.mismatch)

    if (do.annotate) {
      cat('... annotating.\n')
      reads <- annotateReads(reads, annotated.genome,
                             nuke.antisense=nuke.antisense)
      anno.column <- which(colnames(values(reads)) == 'exon.anno')
      colnames(values(reads))[anno.column] <- 'tag.ZZ'
    }

    ## Save result temporarily until all reads are remmapped. Each group of reads
    ## will be saved to separate chromosome files.
    to.chrs <- as.character(unique(seqnames(reads)))

    if (verbose) {
      cat("(", chr, "): will save these tmp files:",
          paste(to.chrs, collapse=','), "\n")
    }

    for (schr in to.chrs) {
      fn <- file.path(path, sprintf('%s-from-%s-to-%s.rda', outfile, chr, schr))
      cat("... (", chr, ") saving: ", fn, "\n", sep="")
      sreads <- reads[seqnames(reads) == schr]
      save(sreads, file=fn)
    }

    to.chrs
  })
  remap.chrs <- unique(unlist(remap.chrs))

  cat("=== Serializing tmp files into final *.sam")
  sam.files <- mclapply(remap.chrs, mc.preschedule=FALSE, function(chr) {
    cat("...", chr, "...\n")
    rda.files <- list.files(path, sprintf('%s.*to-%s.rda', outfile, chr),
                            full.names=TRUE)
    parts <- lapply(rda.files, function(fn) {
      creads <- load.it(fn)
      meta <- values(creads)
      values(creads) <- NULL
      list(granges=creads, meta=meta)
    })
    reads <- do.call(c, lapply(parts, '[[', 'granges'))
    values(reads) <- combineDataFrames(lapply(parts, '[[', 'meta'))

    sam <- toSamTable(reads)
    fn <- file.path(path, paste(chr, outfile, 'sam', sep="."))
    cat("... saving", nrow(sam), "rows to", fn, "\n")
    write.table(sam, fn, sep="\t", col.names=FALSE, row.names=FALSE,
                quote=FALSE, append=FALSE)
    fn
  })

  real.final.fn <- gsub('\\.sam$', '-final.sam', final.fn)
  cat("\n\n=== combining sam files to", real.final.fn, "===\n")

  system(paste("cp", final.fn, real.final.fn))
  for (sf in sam.files) {
    cmd <- paste("cat", sf, ">>", real.final.fn)
    cat("...", cmd, "\n")
    system(cmd)
  }

  cat("Creating BAM file\n")
  asBam(real.final.fn, gsub("\\.sam$", "", real.final.fn))
}

##' Flags events as real (TRUE) or not.
##'
##' A real event is either not internally-primed or it is primed,
##' but has a "valid" pa.signal which is taken from \code{pa.good}
##'
##' @param events A \code{GRanges} object with a \code{values()$is.primed}
##' column.
##' @param pa.good A poly(a) stats data.frame like you would get from
##' isRealPolyA(collectPolyAPositionStatistics(events, ...), events)
flagRealEventsByPolyA <- function(events, pa.good) {
  is.real <- logical(length(events))
  is.primed <- as.logical(values(events)$is.primed)

  is.real[!is.primed] <- TRUE
  colnames(pa.good) <- gsub('.event', '', colnames(pa.good))
  pa.good <- as(pa.good, 'GRanges')
  pa.count <- countOverlaps(events, pa.good)
  is.real[pa.count > 0L] <- TRUE
  is.real
}

## Move multimap reads back to the unique lawn if it makes sense to do so.
realignReads <- function(reads, unique.lawn, max.mm.consider=5L,
                         max.mismatch=2L, best.match.tag='X0',
                         alt.alignment.tag='XA', ...) {
  utag <- paste('tag', best.match.tag, sep=".")
  alt.align.tag <- paste('tag', alt.alignment.tag, sep=".")
  ## automatically set uniquely mapping reads as valid
  unique.reads <- reads[values(reads)[[utag]] == 1]
  values(unique.reads)$tag.Z0 <- 1L

  ## tease out the reads we're interested in.
  reads <- reads[values(reads)[[utag]] > 1]

  too.many <- values(reads)[[utag]] > max.mm.consider
  cat("Dumping", sum(too.many), "reads -- map to too many places.\n")
  reads <- reads[!too.many]

  ## Reads that map to more than one place, but map to more than
  ## `bwa samse -n` number of places have an XA field that is blank.
  ## NOTE: Unfortunately, bwa adds the X0 + X1 alignments, and if they
  ## they are more than the `-n` parameter, the XA tag is left blank!
  unknown.multimap <- is.na(values(reads)[[alt.align.tag]])
  if (any(unknown.multimap)) {
    cat("Removing", sum(unknown.multimap),
        "multimapped reads that have no remapping info\n")
    reads <- reads[!unknown.multimap]
  }

  expanded <- convertBwaMultimap(reads, keep='best',
                                 multimap.column=alt.align.tag)
  if (length(expanded) == 0L) {
    return(unique.reads)
  }

  expanded <- expanded[values(expanded)$edit.distance <= max.mismatch]

  ## Seqinfo's get screwed somehow
  ## fix them for all objects
  all.seqinfo <- merge(seqinfo(unique.lawn), seqinfo(expanded))
  seqinfo(unique.lawn, match(seqlevels(all.seqinfo), seqlevels(unique.lawn))) <- all.seqinfo
  seqinfo(expanded, match(seqlevels(all.seqinfo), seqlevels(expanded))) <- all.seqinfo
  seqinfo(unique.reads, match(seqlevels(all.seqinfo), seqlevels(unique.reads))) <- all.seqinfo

  ## I make a new lawn here in order to "rescue" the case where a remapping
  ## *technically* aligns to more than one place, but really it's the same
  ## place. Imagine R is the remapping positions, and L is the lawn::
  ##
  ##                RRRRRRRRRRRRRRRRRRRRRRRRRR
  ##        LLLLLLLLLLLLLLLLLLLL      LLLLLLLLLLLLLLLLLLLLLL
  ##
  ## The new lawn therefore looks like::
  ##
  ##        LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  new.lawn <- unique.lawn
  xp <- expanded
  values(new.lawn) <- NULL
  values(xp) <- NULL

  new.lawn <- reduce(c(new.lawn, xp))
  new.lawn <- subsetByOverlaps(new.lawn, unique.lawn)

  o <- findOverlaps(expanded, new.lawn)

  if (length(o) == 0L) {
    ## No remapping to unique regions -- let's bail!
    return(unique.reads)
  }

  ## Only keep 1 remapped positions per read, and only if that position
  ## overlaps with some part of our "lawn"
  mm <- as.matrix(o)
  dt <- data.table(query=mm[, 1], subject=mm[, 2],
                   reference=values(expanded)$reference[mm[,1]],
                   key='reference')
  xref <- dt[, list(count=length(query), remap.idx=query[1]), by='reference']
  xref.keep <- xref[count == 1]

  remapped <- expanded[xref.keep$remap.idx]
  r.meta <- values(remapped)
  original <- reads[r.meta$reference]

  ## Deal with new mappings to opposite strand
  flip <- as.logical(strand(remapped) != strand(original))
  r.meta$seq[flip] <- reverseComplement(r.meta$seq[flip])

  ## Clean meta information for the remapped reads
  r.meta[[utag]] <- values(original)[[utag]]
  r.meta$tag.XT <- 'L'
  r.meta$tag.NM <- r.meta$edit.distance
  colnames(r.meta)[which(colnames(r.meta) == 'edit.distance')] <- 'tag.XM'
  r.meta$tag.Z0 <- 1L
  ## Store where this read was remapped FROM
  r.meta$tag.ZA <- toBwaMultimapString(original)

  r.meta <- r.meta[, !colnames(r.meta) %in% c('reference')]

  ## Construct final GRanges object
  ## TODO: Make combineDataFrames faster
  new.meta <- combineDataFrames(values(unique.reads), r.meta)
  values(unique.reads) <- NULL
  values(remapped) <- NULL
  new.reads <- c(unique.reads, remapped)
  values(new.reads) <- new.meta
  ## cat("[combined]", paste(as.character(unique(seqnames(new.reads))), collapse=','), "\n")
  new.reads
}

## Run this function on the the "processed signal" files generated during
## the postProcessAlignments step to generate one "signal file" for the
## entire experiment that goese into datasets/SOMETHING/meta/annotated.signal
##
## After reads have been postprocessed, they'll never be "flipped" -- even if
## the original experiment sequenced 3' -> 5'
reprocessProcessedSignals <- function(expt, annotated.genome,
                                      dataset.name='hump-rescue-1',
                                      seqname.pattern='chr') {
  if (missing(dataset.name)) {
    stop("dataset.name required")
  }
  dataset(expt) <- dataset.name

  signal.path <- file.path(metaDirectory(expt), 'annotated.signal')
  signal.pattern <- sprintf(".*%s.+\\.rda$", seqname.pattern)
  signal.files <- list.files(signal.path, signal.pattern, full.names=TRUE)

  signals <- lapply(signal.files, load.it)
  signals <- suppressWarnings(do.call(c, signals))
  signals <- annotateReads(signals, annotated.genome, nuke.antisense=FALSE)

  ## Quantify reads over events
  chrs <- as.character(unique(seqnames(signals)))
  processed <- mclapply(chrs, mc.preschedule=FALSE, function(chr) {
    reads <- getReadsOnChromosome(expt, chr, unique.only=FALSE, max.mismatch=NULL,
                                  smooth.by=NULL, sequence.flipped=FALSE)
    take <- as.logical(seqnames(signals) == chr)
    if (sum(take) == 0L) {
      return(NULL)
    }
    regions <- signals[take]
    count <- countUniqueOverlaps(regions, reads)
    list(take=take, count=count)
  })

  processed <- processed[sapply(processed, is.list)]

  count <- integer(length(signals))
  for (wut in processed) {
    count[as.logical(wut$take)] <- wut$count
  }
  values(signals)$count <- count
  values(signals)$tpm <- TPM(count, sum(count))
  fn <- file.path(signal.path, 'all.signals.rda')
  cat(".. saving", fn, "\n")
  save(signals, file=fn)
}

## Pass covr.breaks, start.quantile, end.quantile in through ...
## if you want to change the defaults when flagging internal priming
## events
reprocessSmoothedSignals <- function(expt, annotated.genome,
                                     kernel='normal', bandwidth=18L,
                                     chrs=NULL, dataset.name='realign-only',
                                     seqname.pattern='chr',
                                     unique.only=FALSE, max.mismatch=NULL,
                                     sequence.flipped=FALSE,
                                     check.xover=FALSE,
                                     bsg=getBsGenome(expt), ...) {
  verbose <- checkVerbose(...)
  if (missing(dataset.name)) {
    stop("dataset.name required")
  }
  stopifnot(isValidSeqstoreDataset(expt, dataset.name))
  in.path <- file.path(metaDirectory(expt, .dataset=dataset.name),
                       'smooth.signal')
  out.path <- file.path(metaDirectory(expt, .dataset=dataset.name),
                        'annotated.signal')
  if (!checkOrCreateDirectory(out.path, create=TRUE)) {
    stop("Problems creating output directory: ", out.path)
  }
  if (is.null(chrs)) {
    chrs <- extractSeqnames(list.files(in.path))
  }

  expt.read.count <- getReadCount(expt, .dataset=dataset.name)
  if (!is.numeric(expt.read.count)) {
    stop("Experiment readcount not found -- is your metadata.yaml file OK?\n")
  }
  bsg <- getBsGenome(expt)

  processed <- mclapply(chrs, mc.preschedule=FALSE, function(chr) {
    cat(chr, "\n")
    reads <- getReadsOnChromosome(expt, chr, unique.only=unique.only,
                                  max.mismatch=max.mismatch, smooth.by=NULL,
                                  .dataset=dataset.name,
                                  sequence.flipped=sequence.flipped)
    signal <- tryCatch({
      getSmoothSignal(expt, chr, kernel=kernel, bandwidth=bandwidth,
                      .dataset=dataset.name)
    }, error=function(e) NULL)

    if (is.null(signal)) {
      cat("Couldn't retrieve signal for", chr, "\n")
      return(NULL)
    }

    dna <- unmasked(bsg[[chr]])

    expanded <- annotateReads(as(signal, 'GRanges'), annotated.genome,
                              nuke.antisense=FALSE)

    counts <- countUniqueOverlaps(expanded, reads, assign.by='quantify')

    if (check.xover) {
      ## TODO: Handle xover like spp handles pcr?
      covr <- values(expanded)$covr
      covr.count.ratio <- (counts - covr) / covr
      ## 0.5 was just taken empirically
      likely.xover <- abs(covr.count.ratio) > 0.5 & !values(expanded)$is.distal

      if (any(likely.xover)) {
        tx.bounds <- IRanges(values(expanded)$event.start,
                             values(expanded)$event.end)
        likely.xover <- which(likely.xover)
        if (likely.xover[1L] == 1L) {
          consider.prev <- likely.xover[-1L]
        } else {
          consider.prev <- likely.xover
        }
        if (likely.xover[length(likely.xover)] == length(expanded)) {
          consider.next <- likely.xover[-length(likely.xover)]
        } else {
          consider.next <- likely.xover
        }

        cleaned <- expanded
        peak.points <- values(cleaned)$peak.pos

        ## expand endpoint of previous peak if it lies within the same event.bound
        expand.previous <- tx.bounds[consider.prev] == tx.bounds[consider.prev - 1L]
        if (any(expand.previous)) {
          reset <- consider.prev[expand.previous]
          end(cleaned[reset - 1L]) <- peak.points[reset] - 1L
        }

        ## expand start of next peak backwards if its within same event.bound
        expand.next <- tx.bounds[consider.next] == tx.bounds[consider.next + 1L]
        if (any(expand.next)) {
          reset <- consider.next[expand.next]
          start(cleaned[reset + 1L]) <- peak.points[reset] + 1L
        }

        cleaned <- cleaned[-likely.xover]
        counts <- countUniqueOverlaps(cleaned, reads, assign.by='quantify')
        expanded <- cleaned
      }
    }

    values(expanded)$tpm <- TPM(counts, expt.read.count)
    values(expanded)$count <- counts

    ## Pass in
    if (verbose) cat(chr, "calculating quantile bounds\n")
    qcounts <- pmax(values(expanded)$count, values(expanded)$covr)
    bbq <- boundEventsWithSignalByQuantile(expanded, covr=qcounts,
                                           signal.fwd=signal@signal.fwd,
                                           signal.rev=signal@signal.rev, ...)
    values(expanded)$quantile.start <- start(bbq)
    values(expanded)$quantile.end <- end(bbq)

    if (verbose) cat(chr, "calling internal priming\n")
    priming <- flagInternallyPrimedEventsByBoundedQuantiles(expanded, dna, ...)
    values(expanded) <- cbind(values(expanded), as(priming, 'DataFrame'))

    ## save the reannotated object (with signal) into annotated.signals
    ## directory
    reanno <- expanded
    class(reanno) <- 'AnnotatedSignal'
    slot(reanno, 'signal.fwd') <- signal@signal.fwd
    slot(reanno, 'signal.rev') <- signal@signal.rev

    fn <- signalFN(out.path, chr, kernel, bandwidth, unique.only=FALSE,
                   max.mismatch=NULL, ...)
    cat("Saving to", fn, "\n")
    save(reanno, file=fn)
    expanded
  })

  ## Compile the granges inton an uber one and stash into
  ## annotated.signals/all.signals.rda
  processed <- processed[sapply(processed, function(p) !is.null(p))]
  processed <- unlist(GRangesList(processed))

  fn <- file.path(out.path, 'all.signals.rda')
  cat("Saving combined thing to", fn, "\n")
  save(processed, file=fn)
}

################################################################################
## Somehow general utility functions
## -----------------------------------------------------------------------------
combineDataFrames <- function(...) {
  dfs <- unlist(list(...), recursive=FALSE)

  col.names <- unlist(lapply(dfs, colnames))
  col.names <- col.names[!duplicated(col.names)]

  .data.type <- function(colname) {
    for (idx in 1:length(dfs)) {
      if (colname %in% colnames(dfs[[idx]])) {
        column <- dfs[[idx]][[colname]]
        if (is.integer(column)) {
          return(NA_integer_)
        }
        if (is.numeric(column)) {
          return(NA_real_)
        }
        if (is(column)[1] == 'PhredQuality') {
          return(PhredQuality('x'))
        }
        return(NA_character_)
      }
    }
  }

  dfs2 <- lapply(dfs, function(df) {
    columns <- lapply(col.names, function(col) {
      if (col %in% colnames(df)) {
        return(df[[col]])
      } else {
        return(rep(.data.type(col), nrow(df)))
      }
    })
    names(columns) <- col.names
    ## df <- do.call(DataFrame, columns)
    ## df <- as(columns, 'DataFrame')
    ## colnames(df) <- col.names
    ## df
    as(columns, 'DataFrame')
  })
  do.call(rbind, dfs2)
}
