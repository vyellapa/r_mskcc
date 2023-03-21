## TODO: Building flank annotation is still screwed up! Look at RefSeq EIF4A1.
## there are many INTERNAL utr3* locations!
## setAs("GRanges", "AnnotatedChromosome", function(from) {
##   class(from) <- "AnnotatedChromosome"
##   from
## })

#' Checks a GenomicRanges-like object if it looks like an annotated genome
isValidAnnotatedGenome <- function(x, check.extensions=FALSE,
                                   key.by=c('seqnames','strand','entrez.id')) {
  if (!inherits(x, "GenomicRanges")) {
    return("This does not inherit from GenomicRanges")
  }
  errs <- character()
  ##############################################################################
  ## Ensure object is densely annotated and has no overlaps
  x.gaps <- gaps(x)
  x.gaps <- x.gaps[strand(x.gaps) != '*']
  has.gaps <- length(x.gaps) > 0
  has.overlaps <- length(reduce(x, min.gap = 0L)) != length(x)
  if (has.gaps) {
    errs <- "Object has gaps in the annotation"
  }
  if (has.overlaps) {
    errs <- c(errs, "Object has overlapping annotations")
  }

  ###############################################################################
  ## Check that there is a max of 1 utr3* and utr5* for entrez.id
  if (check.extensions) {
    ag.dt <- as(anno, 'data.table')
    setkeyv(ag.dt, key.by)
    n.ext <- ag.dt[, {
      list(utr5e=sum(exon.anno == 'utr5*'), utr3e=sum(exon.anno == 'utr3*'))
    }, by=key.by]
    if (any(n.ext$utr5e > 1)) {
      errs <- c(errs, "utr5* annotations are not consistent")
    }
    if (any(n.ext$utr3e > 1)) {
      errs <- c(errs, "utr3* annotations are not consistent")
    }
  }

  if (length(errs) > 0) {
    return(errs)
  } else {
    return(TRUE)
  }
}


##' Utility method to check if the GRanges-like object is 'annotated'
isAnnotated <- function(x) {
  if (inherits(x, 'GenomicRanges') || inherits(x, 'Ranges')) {
    col.names <- colnames(values(x))
  } else if (inherits(x, 'DataFrame') || inherits(x, 'data.frame')) {
    col.names <- colnames(x)
  } else {
    stop("Don't know how to deal with `x`")
  }

  required <- c('exon.anno', 'entrez.id', 'utr3.index')
  all(required %in% col.names)
}


##' Adds an exon.anno column to a \code{\linkS4class{CompressedReads}} object
##' indicating which part of the genome it lands in.
##'
##' This function assumes that no reads in \code{x} overlaps with eachother,
##' such as reads returned from \code{\link{smoother.repackTags}}. Furthermore,
##' all reads are resized to be of width 1 so that no reads span annotation
##' boundaries (set the \code{fix} parameter accordingly).
##'
##' @exportMethod
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param x A \code{\linkS4class{GRanges}} object
##' @param annotation An \code{\linkS4class{AnnotatedChromosome}} object used as
##' the reference for the annotations.
##' @param fix \code{'start'}, \code{'end'} indicating where to anchor the read
##' when shrinking its width to 1.
##' @param nuke.antisense If \code{TRUE}, removes all potential antisense reads.
##'
##' @return An augmented \code{\linkS4class{CompressedReads}} object with an
##' \code{exon.anno} column in its \code{elementMetadata} indicating the read's
##' annotation.
setGeneric("annotateReads",
function(x, annotation, fix=c('none', 'center', 'end', 'start'), fix.width=1L,
         ignore.strand=FALSE, assign.by=c('quantify', 'fix'), ...) {
  standardGeneric("annotateReads")
})

setMethod("annotateReads", c(x="GRanges"),
function(x, annotation, fix, fix.width, ignore.strand, assign.by, ...) {
  stopifnot(isAnnotated(annotation))
  assign.by <- match.arg(assign.by)
  if (ignore.strand) {
    stop("ignore.strand not implemnted yet")
  }

  if (length(x) == 0) {
    ## meta <- DataFrame(exon.anno=character(), entrez.id=)
    return(x)
  }

  fix <- match.arg(fix)
  if (fix != 'none') {
    assign.by <- 'fix'
    orig.start <- start(x)
    orig.end <- end(x)
    if ('peak.pos' %in% colnames(values(x))) {
      cat("Annotating at peak positions\n")
      start(x) <- values(x)$peak.pos
      end(x) <- values(x)$peak.pos
    } else {
      x <- resize(x, width=fix.width, fix=fix)
    }
  }

  annos <- matchToAnnotation(x, annotation, assign.by)

  if (!ignore.strand) {
    ##############################################################################
    ## Check antisense
    anti.annos <- c('intergenic', 'utr5*', 'utr3*')
    maybe.anti <- annos$exon.anno %in% anti.annos
    base.ga <- annotation[!values(annotation)$exon.anno %in% anti.annos]
    suppressWarnings({
      anti.xref <- matchToAnnotation(x, swapStrand(base.ga), assign.by)
    })
    is.anti <- !anti.xref$exon.anno %in% c('unknown', 'intergenic')
    is.anti <- is.anti & maybe.anti

    anti.xref$exon.anno[is.anti] <- 'antisense'
    annos[is.anti,] <- anti.xref[is.anti,]
  } else {
    ## TODO: annotateReads,ignore.strand=TRUE
    op.anno <- matchToAnnotation(x, swapStrand(annotation), assign.by)
    ## ifelse(anno$exon.anno == 'intergenic', op.anno, anno)
  }

  if (is.null(values(x)) || ncol(values(x)) == 0) {
    values(x) <- annos
  } else {
    for (anno in colnames(annos)) {
      values(x)[[anno]] <- annos[[anno]]
    }
  }

  ## Resize the reads back to the original size they were sent in as.
  if (fix != 'none') {
    start(x) <- orig.start
    end(x) <- orig.end
  }

  x
})

##' Returns the DataFrame of annotation info for the ranges in `ranges`.
##' Helper function to annotateReads
##' @nord
matchToAnnotation <- function(ranges, annotation,
                              assign.by=c('quantify', 'fix')) {
  stopifnot(is(ranges, 'GRanges'))
  stopifnot(isAnnotated(annotation))
  assign.by <- match.arg(assign.by)

  anno.df <- values(annotation)
  if (any(width(ranges) > 1L)) {
    uo <- assignUniqueOverlaps(ranges, annotation, assign.by=assign.by)
  } else {
    o <- findOverlaps(ranges, annotation)
    uo <- rep(NA_integer_, length(ranges))
    uo[queryHits(o)] <- subjectHits(o)
  }

  no.anno <- is.na(uo)
  n.no <- sum(no.anno)
  if (n.no > 0L) {
    ## This can happen if reads land in regions of the genome that we don't
    ## have annotation information for. chrM has no annotation information
    ## using RefSeq annotations, for instance.
    warning(n.no, " reads have no matching annotation!\n",
            "exon.anno is set to 'unknown' -- are they from chrM?")
    dummy.anno <- tryCatch({
      lapply(IRanges::as.list(values(annotation)[1,]), function(anno.col) {
        ## How else can we make an "empty copy" of the DataFrame with
        ## the appropriate columns? (anno.df[0,] would wwork, but then
        ## we have 0 rows)
        as(NA, class(anno.col))
      })
    }, error=function(e) NULL)

    if (is.null(dummy.anno)) {
      stop("Do you have a factor column in your values(annotation)?")
    }
    dummy.anno$exon.anno <- 'unknown'
    dummy.anno <- as(dummy.anno, 'DataFrame')
    anno.df <- rbind(anno.df, dummy.anno)
    uo[is.na(uo)] <- nrow(anno.df)
  }

  anno.df[uo,]
}

if (FALSE) {
##' Annotate the reads in the experiment against an annotated genome/chromosome.
##'
##' @exportMethod
##' @rdname annotateExperiment-methods
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param x The TagSeq experiment
##' @param gcache \code{\link{GenomicCache} object used to annotate the tags.
##' @param gene.collapse The policy used to collapse genes with multiple
##' isoforms into one "universal" gene, through \code{\link{idealized}}.
##' @param flank.up The 5' extension to UTRs for annotation purposes
##' @param flank.down The 3' extension to UTRs.
##' @param smooth.by Strategy used to smooth the reads
##' @param unique.only Annotated uniquely mapping reads, or not.
##' @param chrs Optional list of chromosomes used to subset the tags annotated
##' @param stranded Whether or not to annotated the tags against a stranded, or
##' unstranded reference.
##' @param fix.read Tags are shrunk to be 1 bp. This decides which direction to
##' do the shrinking.
##' @param annotate.ds.rsite Annotate the downstream restriction site?
##' @param do.save A logical indicating whether the annotations should be saved.
##' @param save.path The directory the reads should be saved into. Defaults
##' to \code{x}'s \code{meta/annotated.reads} directory.
##' @param create.dir Create the \code{save.path} if it doesn't already exist?
##'
##' @return Invisibly returns the annotated reads.
setMethod("annotateExperiment", c(x="TagSeq"),
function(x, gcache, gene.collapse='longest', flank.up=1000L, flank.down=1000L,
         smooth.by='repackTags', unique.only=TRUE, split.multimap=FALSE,
         max.split=10L, max.mismatch=0L, chrs=NULL, stranded=TRUE,
         fix.read=c('center', 'end', 'start', 'none'),
         nuke.antisense=FALSE, annotate.ds.rsite=FALSE, do.save=TRUE,
         save.path=NULL, create.dir=TRUE, as.peaks=FALSE,
         annotation.source=NULL, ...) {
  if (missing(gcache)) {
    stop("GenomicCache object needed")
  }
  if (!(inherits(gcache, 'GenomicCache') || inherits(gcache, 'GenomicRanges'))) {
    stop("gcache needs to be a genomic cache or annotatedchromosome object")
  }
  if (is.null(save.path)) {
    save.path <- file.path(metaDirectory(x, ...), 'annotated.reads')
  }
  gene.collapse <- GenomicFeaturesX:::matchGFGeneCollapse(gene.collapse)
  verbose <- checkVerbose(...)
  args <- list(...)
  fix.read <- match.arg(fix.read)

  if (is.null(max.mismatch)) {
    max.mismatch <- metadata(x)$max.mismatch
    if (is.null(max.mismatch)) {
      stop("max.mismatch not found in metadata")
    }
  }

  if (do.save) {
    checkOrCreateDirectory(save.path, create.dir)
  }

  ## Tweak the smoothing?
  if (is(x, 'SAGEseq')) {
    if (missing(smooth.by)) {
      ## Testing for missing(smooth.by) because we may already have a BAM
      ## file that has no anchorSite mismatches by ensuring there are 0
      ## mismatches in the seed during alignment. In this case, the caller
      ## can explicitly set the smooth.by to whatever they like, such as
      ## simply 'repackTags' and skip the anchorSite filter
      smooth.by <- unique(c('anchorSite', smooth.by))
    }
  } else {
    if (as.peaks) {
      smooth.by <- NULL
    }
  }

  if (is.null(chrs)) {
    chrs <- chromosomes(x, unique.only=unique.only)
  }
  chrs <- as.character(chrs)
  bad.chr <- !(chrs %in% as.character(chromosomes(x, unique.only=unique.only)))

  if (any(bad.chr)) {
    stop("Unknown chromosomes: ", chrs[bad.chr])
  }
  bsg <- getBsGenome(x)
  if (is.null(annotation.source)) {
    annotation.source <- annotationSource(gcache)
  }
  annotated <- foreach(chr=chrs, .packages=c("GenomicFeaturesX"),
                       .verbose=verbose) %dopar% {
    if (verbose) {
      cat("====", chr, "====\n")
    }
    if (inherits(annotated.genome, 'GenomicRanges')) {
      chr.anno <- gcache[seqnames(gcache) == chr]
    } else {
      .gc <- duplicate(gcache)
      on.exit(dispose(.gc))
      chr.anno <- tryCatch({
        getAnnotatedChromosome(.gc, chr, gene.collapse, flank.up, flank.down,
                               stranded)
      }, error=function(e) NULL)

      if (is.null(chr.anno)) {
        cat("Annotation not found for", chr, "skipping ...\n")
        return(NULL)
      }
    }

    if (verbose) cat("Getting reads ...\n")

    reads <- getReadsOnChromosome(x, chr, unique.only=unique.only,
                                  split.multimap=split.multimap,
                                  smooth.by=smooth.by, max.split=max.split,
                                  max.mismatch=max.mismatch, verbose=verbose,
                                  ...)
    if (length(reads) == 0) {
      cat("... no reads found for", chr, "\n")
      return(NULL)
    }

    ## values(reads) <- values(reads)[,c('tag.X0', 'tag.X1')]
    values(reads) <- NULL

    if (as.peaks) {
      if (verbose) cat("Peak-ifying...\n")
    } else {
      if (verbose) cat("Compressing reads ...\n")
      cr <- compressReads(reads, with.values=FALSE)
    }

    if (verbose) cat("Annotating reads ...\n")
    .annotated <- annotateReads(cr, chr.anno, fix=fix.read,
                                nuke.antisense=nuke.antisense)

    if (is(x, 'SAGEseq') && annotate.ds.rsite) {
      ## Find the downstream restriction site and annotate that.
      ds.sites <- locateDownstreamRestrictionSite(.annotated, chr.anno, bsg,
                                                  restrictionSite(x))
      ## Reannotate the downstream site
      redo <- cr
      ranges(redo) <- IRanges(ds.sites$ds.site, width=1)
      redo <- annotateReads(redo, chr.anno, fix=fix.read)
      ds.meta <- values(redo)[, -1]
      colnames(ds.meta) <- paste('ds', colnames(ds.meta), sep=".")

      new.anno <- cbind(values(.annotated), ds.sites, ds.meta)

      values(.annotated) <- new.anno
    }

    if (do.save) {
      fn <- annotatedReadsFN(save.path, annotation.source, chr,
                             gene.collapse=gene.collapse, flank.up=flank.up,
                             flank.down=flank.down, unique.only=unique.only,
                             split.multimap=split.multimap, max.mismatch=max.mismatch,
                             stranded=stranded)
      cat("Saving to:", fn, "\n")
      save(.annotated, file=fn)
    }
    .annotated
  }

  suppressWarnings({
    ## Don't warn me about combining GRanges object that have different levels
    ## in the seqlengths
    annotated <- do.call(c, annotated[!sapply(annotated, is.null)])
  })

  invisible(annotated)
})


##' Loads cached annotated reads objects
getAnnotatedReads <- function(x, anno.source, gene.collapse="longest",
                              cds.cover='min', flank.up=1000L, flank.down=1000L,
                              unique.only=TRUE, split.multimap=FALSE,
                              max.mismatch=0L, chrs=NULL,
                              meta.cols=c('count', 'exon.anno', 'symbol',
                                'entrez.id', 'utr3.index'), ...) {
  gene.collapse <- GenomicFeaturesX:::matchGFGeneCollapse(gene.collapse)
  if (is(anno.source, "GenomicCache")) {
    anno.source <- annotationSource(anno.source)
  }
  stopifnot(is(anno.source, 'character'))

  if (is.null(max.mismatch)) {
    max.mismatch <- 0L
  }

  adir <- file.path(metaDirectory(x), 'annotated.reads')
  all.chrs <- is.null(chrs)
  if (all.chrs) {
    chrs <- 'chr\\w+'
  }
  files <- annotatedReadsFN(adir, anno.source, chrs, gene.collapse=gene.collapse,
                            cds.cover=cds.cover, flank.up=flank.up,
                            flank.down=flank.down, unique.only=unique.only,
                            split.multimap=split.multimap,
                            max.mismatch=max.mismatch, stranded=stranded)

  if (all.chrs) {
    files <- list.files(dirname(files), basename(files), full.names=TRUE)
  }
  is.file <- sapply(files, file.exists)

  if (!any(is.file)) {
    stop("No annottated reads found.")
  }

  if (!all.chrs && sum(is.file) != length(chrs)) {
    warning("Chromosomes found are not same as asked for")
  }

  anno <- lapply(files, function(file) {
    a <- load.it(file)
    ## Earlier versions of compressReads let the values() attached
    ## to the internal IRanges object slip through with a DataFrame
    ## of its counts. An update to GenomicRanges made this an illegal
    ## GenomicRanges object. So lets nullify this thing.
    elementMetadata(ranges(a)) <- NULL
    if (!is.null(meta.cols)) {
      values(a) <- values(a)[, colnames(values(a)) %in% meta.cols]
    }
    a
  })

  suppressWarnings({
    ## Don't warn me about combining GRanges object that have different levels
    ## in the seqlengths
    anno <- do.call(c, anno)
  })

  values(anno)$tpm <- TPM(values(anno)$count, x)
  anno
}

annotatedTagDistribution <- function(x, expression.weighted=TRUE, do.plot=TRUE,
                                     title=NULL, expt.name=NULL, min.count=1L,
                                     ...) {
  if (inherits(x, 'SeqStore')) {
    expt.name <- experimentName(x, ...)
    x <- getAnnotatedReads(x)
  }
  if (inherits(x, 'GRanges')) {
    x <- as(x, 'data.table')
    setkeyv(x, c('exon.anno'))
  }
  if (is(x, 'data.frame')) {
    x <- as.data.table(x)
  }
  if (!is.data.table(x)) {
    stop("Need a data.table by now")
  }
  if (!'count' %in% colnames(x)) {
    x$count <- 1L
  }
  if (min.count > 1L) {
    x <- subset(x, count >= min.count)
  }

  if (expression.weighted) {
    counts <- x[, list(count=sum(count)), by='exon.anno']
  } else {
    counts <- x[, list(count=length(count)), by='exon.anno']
  }

  if (do.plot) {
    if (is.null(title)) {
      if (expression.weighted) {
        title <- paste("Expression Weighted Annotated Tag Distribution\n",
                       formatC(sum(counts$count), big.mark=","), " tags",
                       sep="")
      } else {
        title <- paste("Unique Genomic Positions with Aligned Reads\n",
                       formatC(sum(counts$count), big.mark=","), " positions",
                       sep="")
      }
      if (!is.null(expt.name)) {
        title <- paste(title, sprintf("[%s]", expt.name))
      }
    }
    g <- ggplot(as.data.frame(counts), aes(exon.anno, count)) + theme_bw() +
      geom_bar(aes(fill=exon.anno), stat='identity') +
        ylab("Count") + xlab("Exon Annotation") +
          opts(axis.text.x=theme_text(angle=-45, hjust=0, vjust=1),
               title=title)


    print(g)
  }

  invisible(counts)
}
}
