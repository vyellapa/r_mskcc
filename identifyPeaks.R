#' Identifies peaks and quantifies them from a list of bam files
#'
#' To make this run in parallel, make sure to setup and register the
#' foreach backend you are using, eg. `registerDoParallel(...)`, shutting
#' down the parallel functionality (eg. `stopCluster(...)`) is the caller's
#' responsibility
#'
#' @param bamfiles A character vector listing the absolute paths to the BAM
#' files to use for the analysis, or a BamFlieList that serves the same purpose.
#' The *.bai files should be "sister files" to each BAM file, respectively.
#' @param outdir The directory to dump the analysis into. Some sub-dirs will
#' be created as well as temp files.
#' @param chrs If you want to restrict the analysis that is run to a subset
#' of the chromosomes found in the `seqinfo` of `bamfiles`, list the names
#' of the chromosomes here.
#' @param stranded logical indicating whether the peaks should be identified
#' separately per strand.
#' @param bam.params.peak.locs ScanBamParam object that has the appropriate
#' filters set to select for the reads used during peak identification.
#' @param bam.params.peak.quantify ScanBamParam object that has the appropriate
#' filters set to select for the reads used during peak quantification.
summarizeTagClustersFromExperiments <-
  function(bamfiles, outdir='.', chrs=NULL, stranded=TRUE,
           bam.param.peaks=bwaSamseParams(),
           bam.flag.peaks=scanBamFlag(isUnmappedQuery=FALSE),
           bam.filter.peaks=bwaSamseFilter(unique.only=FALSE),
           bam.param.quantify=bwaSamseParams(),
           bam.flag.quantify=scanBamFlag(isUnmappedQuery=FALSE),
           bam.filter.quantify=bwaSamseFilter(unique.only=TRUE),
           assign.by.quantify='unique-quantify',
           bandwidth=34L, resize.reads=FALSE, resize.fix='end', resize.is='max',
           min.count.event.boundary=min(10, length(bamfiles)),
           trim.pos=c(.05, .98), trim.neg=c(.02, .95), min.width=21L,
           col.data=NULL, expt.data=NULL, annotation=NULL,
           .options.multicore=list(preschedule=FALSE), ...) {
  bamfiles <- checkBamFileList(bamfiles)
  si <- unifySeqinfo(bamfiles)

  all.chrs <- si$seqnames

  # Prepare subdirectories for analysis
  subdirs <- c(peaks='peaks')
  all.paths <- sapply(subdirs, function(x) {
    x <- file.path(outdir, x)
    checkOrCreateDirectory(x, create=TRUE)
  }, simplify=FALSE)
  stopifnot(all(is.character(unlist(all.paths))))

  pdir <- all.paths['peaks']
  message("Identifying peaks. Saving temp results in: ", pdir)

  pkgs <- c('biosignals', 'GenomicRanges', 'Rsamtools')
  peaks <- foreach(chr=chrs, .packages=pkgs,
                   .options.multicore=.options.multicore) %dopar% {
    message("===", chr, "===")
    these <-
      identifyPeaksOnChromosome(bamfiles, chr, stranded=stranded,
                                bandwidth=bandwidth, resize.reads=resize.reads,
                                resize.fix=resize.fix, resize.is=resize.is,
                                min.count.event.boundary=min.count.event.boundary,
                                trim.pos=trim.pos, trim.neg=trim.neg,
                                bam.param=bam.param.peaks,
                                bam.flag=bam.flag.peaks,
                                bam.filter=bam.filter.peaks)
    saveRDS(these, file.path(pdir, paste(chr, 'peaks.rds', sep='.')))
    these
  }

  errs <- sapply(peaks, is, 'simpleError')
  if (any(errs)) {
    for (x in peaks[errs]) {
      message("Error during peak calling: ", as.character(x))
    }
    stop("... no point in continuing, bye bye")
  }

  all.peaks <- suppressWarnings(do.call(c, unname(peaks[!sapply(peaks, is.null)])))
  all.peaks <- all.peaks[order(all.peaks)]
  saveRDS(all.peaks, file.path(outdir, "peaks.rds"))

  ## Did we get an annotated genome?
  if (!is.null(annotation)) {
    if (isValidAnnotatedGenome(annotation)) {
      all.peaks <- annotateReads(all.peaks, annotation)
    } else {
      warning("The annotate genome that was passed in looked fishy, ",
              "skipping annotation step", immediate.=TRUE)
    }
  }

  # ----------------------------------------------------------------------------
  # Calculate expression per peak across bamfiles
  message("Calculating expression of peaks per bamfile")
  peak.counts <-
    quantifyPeakExpression(all.peaks, bamfiles, chrs=chrs, stranded=stranded,
                           bam.param=bam.param.quantify,
                           bam.filter=bam.filter.quantify,
                           assign.by=assign.by.quantify)

  col.data <- DataFrame(row.names=colnames(peak.counts[[1]]),
                        library.size=colSums(peak.counts[['count']]))
  se <- SummarizedExperiment(peak.counts, rowData=all.peaks, colData=col.data)
  se
}

#' Note that this will probably do very weird things to spliced reads
identifyPeaksOnChromosome <-
  function(bamfiles, chr, stranded=TRUE, bandwidth=35L,
           resize.reads=NA, resize.fix=c('end', 'start', 'center', 'none'),
           resize.is=c('min', 'max', 'absolute'),
           min.count.event.boundary=min(10, length(bamfiles)),
           trim.pos=c(.05, .98), trim.neg=c(.02, .95),
           bam.param=bwaSamseParams(),
           bam.flag=scanBamFlag(isUnmappedQuery=FALSE),
           bam.which=NULL,
           bam.filter=bwaSamseFilter(unique.only=FALSE)) {
  if (is.numeric(resize.reads) && !is.na(resize.reads)) {
    resize.reads <- as.integer(resize.reads)
    stopifnot(resize.reads > 0)
    if (resize.reads < 20) {
      warning("The width you are resizing reads to is rather small",
              immediate.=TRUE)
    }
    resize.fix <- match.arg(resize.fix)
    resize.is <- match.arg(resize.is)
  } else {
    resize.reads <- NA
  }
  if (resize.fix == 'none') {
    resize.reads <- NA
  }
  if (is.null(bam.param)) {
    bam.param <- bwaSamseParams()
  }
  stopifnot(is(bam.param, 'ScanBamParam'))
  if (is.null(bam.flag)) {
    bam.flag <- scanBamFlag(isUnmappedQuery=FALSE)
  }
  bamfiles <- checkBamFileList(bamfiles, clean.names=FALSE)
  si <- unifySeqinfo(bamfiles)
  chr.info <- subset(si, seqnames == chr)
  if (nrow(chr.info) != 1) {
    stop("Unknown chromosome: ", chr)
  }

  if (is.null(bam.filter)) {
    bam.filter <- function(x) x
  } else if (!is.function(bam.filter)) {
    stop("bam.filter needs to be a function (or NULL) that works on ",
         "gapped alignments, see SeqTools::alignFilter file")
  }

  # Get the reads over the entire chromosome specified
  if (is.null(bam.which)) {
    bam.which <- GRanges(chr, IRanges(1, chr.info$seqlengths))
  }
  stopifnot(is(bam.which, "GenomicRanges"))
  stopifnot(all(1 <= start(bam.which)))
  stopifnot(all(end(bam.which) <= chr.info$seqlengths))

  # Setup flags for each "stranded" pass over the data. Doing it this way
  # maintains the other parameters that were set in the bam.flag passed
  # into this function
  bfm <- bamFlagAsBitMatrix(bam.flag)
  if (stranded) {
    fwd.flags <- bfm
    fwd.flags[, 'isMinusStrand'] <- c(1L, 0L)
    rev.flags <- bfm
    rev.flags[, 'isMinusStrand'] <- c(0L, 1L)
    bam.flags <- list('+'=fwd.flags, '-'=rev.flags)
  } else {
    all.flags <- bfm
    all.flags <- all.flabs[, 'isMinusStrand'] <- c(1L, 1L)
    bam.flags <- list('*'=all.flags)
  }
  bam.flags <- lapply(bam.flags, function(bf) {
    flags <- S4Vectors:::implodeIntBits(bf)
    names(flags) <- c('keep0', 'keep1')
    flags
  })

  # Let it rip
  peaks <- lapply(names(bam.flags), function(.strand) {
    flag <- bam.flags[[.strand]]
    # Fetch reads from all input files and manipulate them to the correct length
    xranges <- lapply(bamfiles, function(bf) {
      # Fetch and manipulate reads to appropriate size
      bamWhich(bam.param) <- bam.which
      bamFlag(bam.param) <- flag
      r <- readGAlignments(bf, param=bam.param)
      if (length(r) == 0) {
        return(NULL)
      }
      r <- bam.filter(r)
      if (length(r) == 0) {
        return(NULL)
      }
      r <- unname(ranges(r))
      if (!is.na(resize.reads)) {
        widths <- switch(resize.is,
                         max=pmin(width(r), resize.reads),
                         min=pmax(width(r), resize.reads),
                         rep(resize.reads, length(r)))
        fix <- resize.fix
        if (.strand == '-' && fix != 'center') {
          fix <- c(start='end', end='start')[fix]
        }
        r <- resize(r, width=widths, fix=fix)
      }
      r
    })

    xranges <- do.call(c, unname(xranges[!sapply(xranges, is.null)]))
    xranges <- xranges[order(xranges)]

    # Do peak calling over a coverage vector
    signal <- coverage(xranges)
    trim <- if (.strand == '-') trim.neg else trim.pos
    info <- detectPeaksByEdges(signal, bandwidth=bandwidth,
                               min.height=min.count.event.boundary, trim=trim)
    meta <- mcols(info)
    mcols(info) <- NULL
    GRanges(chr, info, .strand, meta)
  })
  peaks <- peaks[!sapply(peaks, is.null)]

  if (length(peaks)) {
    peaks <- do.call(c, unname(peaks))
  } else {
    peaks <- GRanges()
  }

  peaks
}

#' Note: This functions is parallelized over chromosome, so each BAM file is
#' processed serially.
quantifyPeakExpression <-
  function(peaks, bamfiles, chrs=NULL, stranded=TRUE,
           assign.by='unique-quantify', bam.param=bwaSamseParams(),
           bam.filter=bwaSamseFilter(unique.only=TRUE)) {
  stopifnot(inherits(peaks, "GenomicRanges"))
  bamfiles <- checkBamFileList(bamfiles, clean.names=TRUE)

  si <- unifySeqinfo(bamfiles)
  ## TODO: check that peaks are within the bounds of the chromosome as
  ##       specified in the seqinfo() of the BAM files

  tabulated <- lapply(bamfiles, function(bam) {
    tabulateReads(peaks, bam, assign.by=assign.by, ignore.strand=!stranded,
                  scan.bam.param=bam.param, filter.fn=bam.filter, chrs=chrs,
                  .parallel=TRUE)
  })

  SimpleList(count=sapply(tabulated, '[[', 'count'),
             percent=sapply(tabulated, '[[', 'percent'))
}

# ------------------------------------------------------------------------------
# Helper functions for parameter checking and whatnot

#' Handles parameter checking beuracracy for a list of "things" that should
#' be a BamFileList
checkBamFileList <- function(bamfiles, clean.names=TRUE) {
  if (is.character(bamfiles)) {
    bamfiles <- BamFileList(bamfiles)
  }
  if (!is(bamfiles, 'BamFileList')) {
    stop("BamFileList required")
  }
  if (!all(file.exists(path(bamfiles)))) {
    lost <- !file.exist(path(bamfiles))
    stop(paste(path(bamfiles[lost]), collapse=","), " BAMs not found")
  }

  if (clean.names) {
    nms <- names(bamfiles)
    if (is.null(nms)) {
      nms <- path(bamfiles)
    }
    nms <- gsub("\\.bam$", "", basename(nms))
    names(bamfiles) <- make.unique(nms)
  }

  bamfiles
}

unifySeqinfo <- function(bamfiles) {
  stopifnot(is(bamfiles, "BamFileList"))
  # Ensure that all the seqlengths from the bam headers agree with each other
  all.si <- rbindlist(lapply(names(bamfiles), function(x) {
    si <- as.data.frame(seqinfo(bamfiles[[x]]))
    si <- transform(si, seqnames=rownames(si), source=x)
    as.data.table(si)
  }))
  setkeyv(all.si, 'seqnames')
  si <- all.si[, {
    lens <- unique(seqlengths)
    n.obs <- length(lens)
    circular <- (isCircular[!is.na(isCircular)])[1L] # will be NA if all
    genome <- (genome[!is.na(genome)])[1L]           # elements are NA
    list(seqlengths=lens[1L], isCircular=circular, genome=genome, n.obs=n.obs)
  }, by='seqnames']
  fail <- si$n.obs > 1
  if (any(fail)) {
    stop("mismatched seqlengths for ", paste(si$seqnames[fail], collapse=","))
  }
  si
}

if (FALSE) {
  ## Dicer 3'utr
  ## chr14:95,552,442-95,556,977
  bampath <- "/Users/stavros/cBio/papers/apa-1/data/bams/multimap/clean"
  bamfiles <- file.path(bampath, c('bcells-1.bam', 'brain.bam'))
  bamfiles <- BamFileList(bamfiles)
}
