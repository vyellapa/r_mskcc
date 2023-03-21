##' Takes a SummarizedExperiment and annotates it up for ApaDiff usage
##'
##' TODO: To support intronic ApA analysis, add exon.number to cleavage site as
##'       well as total number of exons for the locus.
##'
##' @export
##' @param SummarizedExperiment The experimental data to annotate
##' @param by.key Columns required to isolate all peaks for 1 transcript
##' @param order.key Order the peaks this way, this is different from by.key
##' @param filter.count Minimum number of reads per expt required to pass a filter
##' @param filter.fraction Minimum number of expts require to pass `filter.count`
##' to pass through
##' @param rm.single Remove genes w/ single peaks
##' @param assay.idx The assay index in the SummarizedExperiment that peak counts
##' are held in
##'
##' @return A SummarizedExperiment with annotated rowRanges() required to convert
##' it to an ApaCountSet
prepSE4ApaTesting <- function(x, by.key=c('seqnames', 'strand', 'entrez.id'),
                              order.key=c(by.key, 'start'),
                              ## filter.count=10, filter.fraction=0.05,
                              ## filter.count=NULL, filter.fraction=NULL,
                              rm.single=FALSE, assay.idx=1L, ...) {
  stopifnot(inherits(x, 'SummarizedExperiment'))

  dt <- as(rowRanges(x), 'data.table')
  dt <- dt[,row.idx:=1:nrow(x)]
  setkeyv(dt, order.key)

  ## Annotate peaks with some ordering information before we filter them
  add.pa.order <- is.null(dt$pa.order)
  if (add.pa.order) {
    dt[, pa.order := calculateUtr3CleavageOrder(dt, by.key=by.key, order.key=order.key)]
  }

  add.peak.id <- is.null(dt$peak.id)
  if (add.peak.id) {
    #dt[, peak.id := ApADiff:::calculatePeakId(dt, by.key, order.key)]
    peak.id <- calculatePeakIdLocal(dt, by.key, order.key)
    dt[, peak.id := peak.id]
  }

  add.exon.id <- is.null(dt$exonID)
  if (add.exon.id) {
    dt[, exonID := peak.id]
  }

  add.group.id <- is.null(dt$group.id)
  if (add.group.id) {
    dt$group <- 'genic'
    dt$group <- ifelse(dt$exon.anno %in% c("utr3", "utr3*"), "utr3", dt$group)
    dt$group <- ifelse(dt$exon.anno == 'intron.utr3', 'intron.utr3', dt$group)
  }

  add.geneID <- is.null(dt$geneID)
  if (add.geneID) {
    dt$geneID <- paste(dt$entrez.id, dt$group, sep=".")
  }

  if (rm.single) {
    dt <- dt[, if (.N == 1) .SD[0] else .SD, by=c(by.key, 'group')]
  }

  ## ---------------------------------------------------------------------------
  ## Filter out peaks in different ways
  cnts <- assay(x, assay.idx)[dt$row.idx,]

  ## Always filter out 0 count events
  keep <- rowSums(cnts) > 0
  dt <- dt[keep,]
  cnts <- cnts[keep,]

  ## Perhaps the user wanted somethig more elaborate?
  # if (is.numeric(filter.count) && is.numeric(filter.fraction)) {
  #   if (filter.fraction < 1) {
  #     ns2pass <- max(2, floor(filter.fraction * ncol(x)))
  #   } else {
  #     ns2pass <- filter.fraction
  #   }
  #   ns.passed <- rowSums(cnts > filter.count) >= ns2pass
  #   if (!all(ns.passed)) {
  #     message("Removing ", sum(!ns.passed), " cleavage events")
  #     dt <- dt[ns.passed,]
  #   }
  # }

  dt <- dt[order(dt$row.idx),]

  ans <- x[dt$row.idx,]
  if (add.pa.order) {
    mcols(rowRanges(ans))$pa.order <- dt$pa.order
  }
  if (add.peak.id) {
    mcols(rowRanges(ans))$peak.id <- dt$peak.id
  }
  if (add.exon.id) {
    mcols(rowRanges(ans))$exonID <- dt$exonID
  }
  if (add.group.id) {
    mcols(rowRanges(ans))$group <- dt$group
  }
  if (add.geneID) {
    mcols(rowRanges(ans))$geneID <- dt$geneID
  }

  ans
}



###################################################################
#This function takes the summarized experiment object of the peaks 
#with the atlas and removes the blacklisted regions and particular
#exon.anno(eg: antisense) or vector of exon.anno(from count and granges)
#It adjusts the library size accordingly
## se.peaks <- summarized experiment object
## blackisted <- logical(indicates whether it should remove blacklister 
##  peaks or not)
## remove.exon.anno <- vector, tells which exon.anno has to be removed
###################################################################
removePeaks <- function(se.peaks, blacklisted = TRUE, remove.exon.anno = NULL, 
  ip.axe = TRUE, change.library.size = TRUE, genome = NULL){

  if(genome %in% c("hg19", "hg38", "mm9", "mm10")){
    blacklist.file <- system.file("extdata", sprintf("%s.blacklist.bed", genome),
                            package="TagSeq")
    black <- import(con = blacklist.file, format = "bed")
    gr.black <- as(black, "RangedData")    
  } else {
    stop("No blacklist bed file present")
  }

  gr.peaks <- rowRanges(se.peaks)
  vec.enames <- rownames(colData(se.peaks))
  vec.remove.all <- c()
  ##' Finding the overlapping peaks in the blacklisted regions
  if(blacklisted){
    vec.remove <- subjectHits(findOverlaps(query = gr.black, subject = gr.peaks))
    vec.remove.all <- c(vec.remove.all, vec.remove)
  }

  if(length(remove.exon.anno) > 0){
    vec.remove <- which(gr.peaks$exon.anno %in% remove.exon.anno)
    vec.remove.all <- c(vec.remove.all, vec.remove)
  }

  if(ip.axe){
    vec.remove <- which(gr.peaks$ip.axe & !gr.peaks$annotated.end)
    vec.remove.all <- c(vec.remove.all, vec.remove)
  }

  vec.remove.all <- unique(vec.remove.all)
  
  ## new summarized experiment object
  se.new.peaks <- se.peaks[-vec.remove.all]

  ##' Changing the library size
  if(change.library.size){
    m.count <- assays(se.peaks)$count
    vec.blacklisted.lib <- colSums(m.count[vec.remove.all, vec.enames])
    vec.original.lib <- colData(se.peaks)[vec.enames, "library.size"]
    names(vec.original.lib) <- vec.enames
    vec.new.lib <- vec.original.lib - vec.blacklisted.lib
  } else {
    vec.new.lib <- colData(se.peaks)[vec.enames, "library.size"]
  }

  ## Changing the library size
  colData(se.new.peaks)[vec.enames, "library.size"] <- vec.new.lib

  se.new.peaks

}

##' Cluster peaks that are close together, optinaly add reads together.
##'
##' @param x The SummarizedExperiment to do this to
##' @param max.distance The max distance to consider combining peaks together
##' @param peak.pos Where to measure the distance between peaks from? Currently
##' only support end,start,center. edges is to use 'end' to measure dist
##' downstream, and start to measure dist upstream
##' @param combine.counts add the read count from each merged peak togther?
clusterPeaks <- function(x, max.distance=50L,
                         peak.pos=c('end', 'start', 'center', 'edges'),
                         requantify=TRUE) {
  peak.pos <- match.arg(peak.pos)
  if (peak.pos == 'edges') {
    stop("`edges` not supported yet")
  }
  stopifnot(inherits(x, 'SummarizedExperiment'))
  stopifnot(is.numeric(max.distance) && max.distance > 0L)

  xevts <- resize(rowRanges(x), width=max.distance, fix=peak.pos)
  mcols(xevts) <- NULL
  mask <- reduce(xevts[order(xevts)])

  mm <- as.matrix(findOverlaps(xevts, mask))
  if (any(duplicated(mm[,1]))) {
    stop("There shouldn't be any duplicate queryHits")
  }

  ## Reasign fenceposts from original peak definitions
  xx <- rowRanges(x)
  mcols(xx) <- NULL
  xx$n.signal <- rowRanges(x)$n.signal
  xx <- as(xx, 'data.table')

  mask.idx <- integer(nrow(mm))
  mask.idx[mm[,1]] <- mm[,2]
  xx[, mask.idx := mask.idx]
  setkeyv(xx, 'mask.idx')

  bounds <- xx[, {
    list(seqnames=seqnames[1L], start=min(start), end=max(end),
         strand=strand[1L], npeaks=.N, n.signal=sum(n.signal))
  }, by="mask.idx"]

  out <- as(bounds[, mask.idx := NULL], 'GRanges')
  if (requantify) {
    out <- newSEfromRegions(out, x)
  }
  out
}

## expression
getTPM <- function(se){
  enames <- colnames(se)
  m.count <- assays(se)$count
  vec.lib <- colData(se)[enames, "library.size"]
  names(vec.lib) <- enames
  m.tpm <- sapply(enames, function(name) m.count[ ,name]*10^6/vec.lib[name])
  return(m.tpm)
}

##' Calculate the usage of a given site in the tx unit
calculateUTR3AndIntronUsage <- function(x, fraction.of=c('total', 'max'),
                           with.events=FALSE, defactor=TRUE,
                           group.by=c('seqnames', 'strand', 'entrez.id'),
                           assay.idx=1L, usage.on = "group") {
  stopifnot(is(x, "SummarizedExperiment"))

  fraction.of <- match.arg(fraction.of)
  agg <- switch(fraction.of, total=sum, max=max)

  dt <- cbind(as(rowRanges(x), 'data.table'), data.table(assay(x, assay.idx)))
  dt[, .idx. := seq(nrow(x))]
  setkeyv(dt, group.by)

  if (!usage.on %in% colnames(dt)) {
    warning("All events in the tx-unit are considered together -- call `addProcessingGroup`")
    dt[, usage.on := 'unk']
  }

  keep <- complete.columns(dt, group.by)
  if (!all(keep)) {
    nuke <- !keep
    warning("There are", sum(nuke), "NA keyed events", immediate.=TRUE)
  }

  expt.cols <- colnames(x)

  ans <- matrix(NA_real_, nrow=nrow(x), ncol=ncol(x))
  colnames(ans) <- colnames(x)

  ##This step is crucial because it will make sure that the intronic usage
  ##is calculated before the utr3 usage
  groups <- sort(unique(dt[[usage.on]]))

  for (pgroup in groups) {
    if(grepl("intronic", pgroup)){
      dsub <- subset(dt, dt[[usage.on]] %in% c(pgroup, "utr3"))
    } else {
      dsub <- subset(dt, dt[[usage.on]] == pgroup)
    }
    u <- dsub[, {
      out <- sapply(expt.cols, function(e) .SD[[e]] / agg(.SD[[e]]), simplify=FALSE)
      out$.idx. <- .idx.
      out
    }, by=group.by]
    ans[u$.idx., expt.cols] <- as.matrix(u[, expt.cols, with=FALSE])
  }
  

  ans[is.nan(ans)] <- NA

  assay(x, 'usage') <- ans
  x
}

#' Get assay data and stick it to the right event
#'
#' @param x A data.frame like object w/ seqnames,start,end,strand info to
#' identify the correct event in `atlas`
#' @param atlas The atlas to get the data from
#' @param assay.idx The assay to fetch from the atlas
#' @param x.melted Indicates if each event has a unique row, and the data for
#' each sample must be appended as columns (FALSE), or if TRUE each event is
#' represented multiple times and each row represents the value for the assay
#' in a particular sample indicated by the `melted.column`.
#' @param melted.column The name of the column in `x` that indicates the
#' experiment/condition/whatever that each stat applies to.
#'
#' @return An augmented version of `x` with the assay data attached to each
#' event.
alignAssayData <- function(x, atlas, assay.idx, x.melted=TRUE,
                           melted.column="sample",
                           combine.by=c('mean', 'median', 'sum', 'none'),
                           key.by=c('seqnames', 'start', 'end', 'strand')) {
  combine.by <- match.arg(combine.by)
  atlas.conds <- unique(as.character(colData(atlas)$condition))
  if (!missing(x.melted) && isTRUE(x.melted)) {
    stopifnot(melted.column %in% names(x))
  }
  x <- as.data.table(x)
  changed.seqnames.col <- !"seqnames" %in% names(x) && "chr" %in% names(x)
  if (changed.seqnames.col) {
    setnames(x, "chr", "seqnames")
  }
  setkeyv(x, key.by)

  looks.melted <- local({
    tmp <- head(x, 100)
    dup.cnt <- tmp[, list(N=.N), by=key.by]
    any(dup.cnt$N > 1)
  })

  if (!x.melted == looks.melted) {
    stop("'melted' status of `x` doesn't jibe w/ the value the x.melted arg")
  }

  if (combine.by == "none") {
    stop("You need to combine replicates")
  }

  ad <- as.data.table(getAssay(atlas, assay.idx, combine.by=combine.by))
  assay.name <- if (is.numeric(assay.idx)) {
    names(assays(atlas))[assay.idx]
  } else {
    assay.idx
  }

  if (x.melted) {
    ## Need to explicitly call into melt.data.frame for the setting of
    ## variable.name to work correctly
    mad <- reshape2:::melt.data.frame(as.data.frame(ad), key.by,
                                      names(ad)[names(ad) %in% atlas.conds],
                                      variable.name=melted.column)
    mad <- as.data.table(mad)
    ###The unique line added
    mad <- unique(mad)
    setkeyv(mad, key.by)
    mad[[melted.column]] <- as.character(mad[[melted.column]])
    setnames(mad, 'value', assay.name)
    out <- merge(x, mad, by=c(key.by, melted.column), all.x=TRUE)
  } else {
    ad <- ad[, names(ad) %in% c(key.by, atlas.conds), with=FALSE]
    setnames(ad, atlas.conds, paste(atlas.conds, assay.name, sep='.'))
    out <- merge(x, ad, by=key.by, all.x=TRUE)
  }

  if (nrow(out) != nrow(x)) {
    warning("The number of rows coming out here is wrong", immediate.=TRUE)
  }

  if (changed.seqnames.col) {
    setnames(out, 'seqnames', 'chr')
  }

  out
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
annotateReads <- function(x, annotation, fix=c('none', 'center', 'end', 'start'), fix.width=1L,
         ignore.strand=FALSE, assign.by=c('quantify', 'fix'), ...) {
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
    if ('peak.pos' %in% colnames(rowRanges(x))) {
      cat("Annotating at peak positions\n")
      start(x) <- rowRanges(x)$peak.pos
      end(x) <- rowRanges(x)$peak.pos
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
    base.ga <- annotation[!mcols(annotation)$exon.anno %in% anti.annos]
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

  if (is.null(mcols(x)) || ncol(mcols(x)) == 0) {
    mcols(x) <- annos
  } else {
    for (anno in colnames(annos)) {
      mcols(x)[[anno]] <- annos[[anno]]
    }
  }

  ## Resize the reads back to the original size they were sent in as.
  if (fix != 'none') {
    start(x) <- orig.start
    end(x) <- orig.end
  }

  x
}

##' Returns the DataFrame of annotation info for the ranges in `ranges`.
##' Helper function to annotateReads
##' @nord
matchToAnnotation <- function(ranges, annotation,
                              assign.by=c('quantify', 'fix')) {
  stopifnot(is(ranges, 'GRanges'))
  stopifnot(isAnnotated(annotation))
  assign.by <- match.arg(assign.by)

  anno.df <- mcols(annotation)
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
      lapply(IRanges::as.list(mcols(annotation)[1,]), function(anno.col) {
        ## How else can we make an "empty copy" of the DataFrame with
        ## the appropriate columns? (anno.df[0,] would wwork, but then
        ## we have 0 rows)
        as(NA, class(anno.col))
      })
    }, error=function(e) NULL)

    if (is.null(dummy.anno)) {
      stop("Do you have a factor column in your rowRanges(annotation)?")
    }
    dummy.anno$exon.anno <- 'unknown'
    dummy.anno <- as(dummy.anno, 'DataFrame')
    anno.df <- rbind(anno.df, dummy.anno)
    uo[is.na(uo)] <- nrow(anno.df)
  }

  anno.df[uo,]
}

##' Cluster peaks that are close together, optinaly add reads together.
##'
##' @param x The SummarizedExperiment to do this to
##' @param max.distance The max distance to consider combining peaks together
##' @param peak.pos Where to measure the distance between peaks from? Currently
##' only support end,start,center. edges is to use 'end' to measure dist
##' downstream, and start to measure dist upstream
##' @param combine.counts add the read count from each merged peak togther?
clusterPeaks <- function(x, max.distance=50L,
                         peak.pos=c('end', 'start', 'center', 'edges'),
                         requantify=TRUE) {
  peak.pos <- match.arg(peak.pos)
  if (peak.pos == 'edges') {
    stop("`edges` not supported yet")
  }
  stopifnot(inherits(x, 'SummarizedExperiment'))
  stopifnot(is.numeric(max.distance) && max.distance > 0L)

  xevts <- resize(rowRanges(x), width=max.distance, fix=peak.pos)
  mcols(xevts) <- NULL
  mask <- reduce(xevts[order(xevts)])

  mm <- as.matrix(findOverlaps(xevts, mask))
  if (any(duplicated(mm[,1]))) {
    stop("There shouldn't be any duplicate queryHits")
  }

  ## Reasign fenceposts from original peak definitions
  xx <- rowRanges(x)
  mcols(xx) <- NULL
  xx$n.signal <- rowRanges(x)$n.signal
  xx <- as(xx, 'data.table')

  mask.idx <- integer(nrow(mm))
  mask.idx[mm[,1]] <- mm[,2]
  xx[, mask.idx := mask.idx]
  setkeyv(xx, 'mask.idx')

  bounds <- xx[, {
    list(seqnames=seqnames[1L], start=min(start), end=max(end),
         strand=strand[1L], npeaks=.N, n.signal=sum(n.signal))
  }, by="mask.idx"]

  out <- as(bounds[, mask.idx := NULL], 'GRanges')
  if (requantify) {
    out <- newSEfromRegions(out, x)
  }
  out
}


newSEfromRegions <- function(regions, original, change = c("count", "count.nomismatch", "tpm")) {
  stopifnot(inherits(regions, 'GenomicRanges'))
  stopifnot(inherits(original, 'SummarizedExperiment'))

  mm <- as.matrix(findOverlaps(regions, rowRanges(original)))

  sl.new.assays <- SimpleList()
  for(i in names(assays(original))){
    m.assay <- assays(original)[[i]]
    dt <- data.table(region.id=mm[,1], as.matrix(m.assay)[mm[,2], ], key='region.id')
    if( i %in% change){
      dt.new <- dt[, lapply(.SD, sum), by='region.id']
    } else {
      dt.new <- dt[, .SD[1L, ], by='region.id']
    }
    new.regions <- regions[dt.new[["region.id"]]]
    m.new.assay <- as.matrix(dt.new[ , region.id := NULL])
    sl.new.assays[[i]] <- m.new.assay
  }

  se <- SummarizedExperiment(sl.new.assays, rowRanges=new.regions,
                             colData=colData(original))
  se

}

##' Utility method to check if the GRanges-like object is 'annotated'
isAnnotated <- function(x) {
  if (inherits(x, 'GenomicRanges') || inherits(x, 'Ranges')) {
    col.names <- colnames(mcols(x))
  } else if (inhierts(x, 'DataFrame') || inherits(x, 'data.frame')) {
    col.names <- colnames(x)
  } else {
    stop("Don't know how to deal with `x`")
  }

  required <- c('exon.anno', 'entrez.id', 'utr3.index')
  all(required %in% col.names)
}

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
    ## Earlier versions of compressReads let the rowRanges() attached
    ## to the internal IRanges object slip through with a DataFrame
    ## of its counts. An update to GenomicRanges made this an illegal
    ## GenomicRanges object. So lets nullify this thing.
    elementMetadata(ranges(a)) <- NULL
    if (!is.null(meta.cols)) {
      rowranges(a) <- rowRanges(a)[, colnames(rowRanges(a)) %in% meta.cols]
    }
    a
  })

  suppressWarnings({
    ## Don't warn me about combining GRanges object that have different levels
    ## in the seqlengths
    anno <- do.call(c, anno)
  })

  rowranges(anno)$tpm <- TPM(rowRanges(anno)$count, x)
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

##' Assigns each range in query to a unique range in subject.
##'
##' The ranges in \code{subject} cannot be overlapping!
##'
##' @param query *Ranges object
##' @param subject *Ranges object
##' @param assign.by how do you want the unique assignment to be caluclated?
##' @param fix The param to pass to resize if \code{assign.by == 'fix'}
##'
##' @return An integer vector giving the index into \code{subject} that
##' each element in \code{query} is assigned to. ranges that cannot be
##' assigned are \code{NA}.
assignUniqueOverlaps <- function(query, subject, assign.by=c('quantify', 'fix'),
                                 maxgap=0L, minoverlap=1L,
                                 fix=c('start', 'end', 'center'),
                                 .subject.overlap.checked=TRUE, ...) {
  if (!.subject.overlap.checked) {
    o <- findOverlaps(subject, subject, type='any', ignoreSelf=TRUE,
                      ignoreRedundant=TRUE)
    if (length(o) > 0L) {
      stop("subject ranges cannot overlap.")
    }
  }
  fix <- match.arg(fix)
  assign.by <- match.arg(assign.by)
  ans <- rep(NA_integer_, length(query))

  if (assign.by == 'fix') {
    query <- resize(query, width=1L, fix=fix)
    o <- findOverlaps(query, subject, maxgap, minoverlap)
    ans[queryHits(o)] <- subjectHits(o)
  } else {
    qo <- quantifyOverlaps(query, subject, maxgap, minoverlap)
    qo <- data.table(qo, key="queryHits") ## key by query

    qo.unique <- qo[, list(subjectHits=subjectHits[which.max(p.overlap)]),
                    by=key(qo)]

    ans[qo.unique$queryHits] <- as.integer(qo.unique$subjectHits)
  }

  ans
}

quantifyOverlaps <- function(query, subject, maxgap, minoverlap, ...) {

  ##############################################################################
  ## TODO: Fix quantifyOverlaps to work on circular chromosomes
  ##       Had a hell of a time chasing a bug down from a read that went
  ##       "circular" around chrM ... fix ends of query to be no longer
  ##       than the end of the chromosome
  ##
  ##       I think this error only ever happened because the reads object
  ##       didn't have any proper rowRanges for seqlengths, so the "overhang"
  ##       wasn't noticed.
  seqlvls <- intersect(seqlevels(query), seqlevels(subject))
  is.circular <- union(names(isCircular(subject))[isCircular(subject)],
                       names(isCircular(query))[isCircular(query)])
  is.circular <- is.circular[is.circular %in% seqlvls]
  if (length(is.circular) > 0L) {
    seqlens <- pmin(seqlengths(query)[is.circular],
                    seqlengths(subject)[is.circular], na.rm=TRUE)
    for (sname in names(seqlens)) {
      len <- as.integer(seqlens[sname])
      end(query[seqnames(query) == sname]) <-
        pmin.int(end(query[seqnames(query) == sname]), len)
      end(subject[seqnames(subject) == sname]) <-
        pmin.int(end(subject[seqnames(subject) == sname]), len)
    }
  }
  ## ---------------------------------------------------------------------------

  o <- findOverlaps(query, subject, maxgap, minoverlap)
  overlapped <- ranges(o, ranges(query), ranges(subject))
  p.overlap <- width(overlapped) / width(query)[queryHits(o)]
  p.overlap <- ifelse(p.overlap > 1, 1, p.overlap)
  cbind(as.matrix(o), p.overlap)
}

#########
annotateSummarizedExperiment <- function(x, annotations, nuclear=TRUE) {
  stopifnot(inherits(x, "SummarizedExperiment"))
  stopifnot(inherits(annotations, "GenomicRanges"))

  regions <- rowRanges(x)
  xmeta <- mcols(regions)
  mcols(regions) <- NULL

  anno.cols <- names(mcols(annotations))

  if (nuclear) {
    for (wut in anno.cols) {
      if (wut %in% names(xmeta)) {
        xmeta[[wut]] <- NULL
      }
    }
  }

  mcols(regions) <- xmeta
  suppressWarnings(regions <- annotateReads(regions, annotations))
  rowRanges(x) <- regions
  x
}

###Making  more groupls than that defined by steve
calculateUtr3CleavageOrder <- function(x, by.key=c('seqnames', 'strand', 'entrez.id'),
           order.key=c(by.key, 'start'), groups=NULL) {
  if (is.null(groups)) {
    groups <- list(utr3=c('utr3', 'utr3*'), intron.utr3=c('intron.utr3'))
  }

  class <- class(x)
  if (inherits(x, 'TagAtlas')) {
    xx <- as(x@events, 'data.table')
  } else if (inherits(x, 'SummarizedExperiment')) {
    xx <- as(rowRanges(x), 'data.table')
  } else if (inherits(x, "GenomicRanges")) {
    xx <- as(x, 'data.table')
  } else if (is(x, 'data.frame')) {
    xx <- as.data.table(x)
  }
  stopifnot(is(xx, "data.table"))

  xx$.order. <- 1:nrow(x)
  setkeyv(xx, order.key)

  pa.order <- rep("unk", nrow(x))

  for (group in groups) {
    #x.sub <- subset(xx, exon.anno %in% group)
    pos <- which(xx$exon.anno %in% group)
    x.sub <- xx[pos, ]
    if (nrow(x.sub) > 0) {
      x.pa <- x.sub[, {
        if (length(start) == 1L) {
          pao <- 'single'
        } else {
          pao <- rep('internal', .N)
          if (strand[1] == '+') {
            pao[1L] <- 'proximal'
            pao[.N] <- 'distal'
          } else {
            pao[.N] <- 'proximal'
            pao[1L] <- 'distal'
          }
        }
        list(pa.order=pao, .order.=.order.)
      }, by=by.key]
      pa.order[x.pa$.order.] <- x.pa$pa.order
    }
  }

  pa.order
}


calculatePeakIdLocal <- function(x, by.key=c('seqnames', 'strand', 'entrez.id'),
                            order.key=c(by.key, 'start'), peak.prefix='pas') {
  if (inherits(x, 'SummarizedExperiment')) {
    x <- as(rowRanges(x), 'data.frame')
  }
  if (inherits(x, "GenomicRanges")) {
    x <- as.data.frame(x)
  }
  if (!inherits(x, 'data.frame')) {
    stop("data.frame required by now")
  }

  new.x <- as.data.table(x)
  x <- new.x[ ,.order. := 1:nrow(x)]

 
  setkeyv(x, order.key)
  ans <- x[, {
    list(peak.id=if (strand == '+') 1:.N else .N:1, .order.=.order.)
  }, by=by.key]

  paste(peak.prefix, ans[order(ans$.order.), ]$peak.id, sep='.')
}

###########################################
callDifferentialUsedIntronicEvents <- function(se.atlas, enames, df.design = NULL, vs = NULL,
  comparison.name, utr3.exprs = TRUE,
  peak.tpm.cutoff = 5L, utr3.tpm.cutoff = 5.5, percent.exprs = 0.75, sf.for.all = T){

    
  ##Subset the experiments to look at
  se.prep <- se.atlas[ ,enames]
  rowRanges(se.prep)$idx <- 1:length(rowRanges(se.prep))
  
  se.intron <- se.prep[rowRanges(se.prep)$cleavage.unit == "intronic"]
    
  ##The intronic peak should be expressed in atleast one condition with the  given cutoff
  all.cond <- as.character(unique(colData(se.intron)$condition))

  m.exprs.logical <- sapply(all.cond, function(x) {
    pos <- which(colData(se.intron)$condition == x)
    m.tpm <- as.matrix(assays(se.intron)$tpm[ ,pos])
    l <- apply(m.tpm, 1, function(idx) {
      length(which(idx > peak.tpm.cutoff)) >= (percent.exprs*length(idx))
      })
    return(l)
    })
  exprs.intron <- apply(m.exprs.logical, 1, function(x) any(x))
  se.exprs.intron <- se.intron[exprs.intron]

  if(utr3.exprs){
    ######
    ##We can test only the intronic peaks that have atleast intron or utr3 expression 
    ##in all conditions
    m.exprs.all.logical <- sapply(all.cond, function(x) {
      pos <- which(colData(se.exprs.intron)$condition == x)
      m.intron.tpm <- as.matrix(assays(se.exprs.intron)$tpm[ ,pos] > peak.tpm.cutoff)
      m.utr3.tpm <- as.matrix(assays(se.exprs.intron)$utr3.tpm[ ,pos] > utr3.tpm.cutoff)
      exp <- apply(m.intron.tpm | m.utr3.tpm, 1, function(idx){
        length(which(idx)) >= (percent.exprs*length(idx))
        })
      })
    exprs.intron.all <- apply(m.exprs.all.logical, 1, function(x) all(x))
    se.intron.to.test <- se.exprs.intron[exprs.intron.all]
      
  } else {

    se.intron.to.test <- se.exprs.intron

  }

  ###Number of introns to be tested
  print(length(rowRanges(se.intron.to.test)))
  
  ##Now select the utr3 and intron eventss for the genes that have intronic peaks
  se.exprs.utr3.intron <- se.prep[(rowRanges(se.prep)$cleavage.unit == "utr3" &
    rowRanges(se.prep)$entrez.id %in% unique(rowRanges(se.intron.to.test)$entrez.id)) |
    (rowRanges(se.prep)$idx %in% rowRanges(se.intron.to.test)$idx)]

  
  ##Since evry intronic peak is tested against its full length so we
  ##need to change the gene id for the full length
  hg19.annotation.file <- "/Users/singhi1/ifs/annotation_data/hg19/hg19_annotation_pfam.rds"
  hg19.annotation <- readRDS(hg19.annotation.file)
  se.exprs.utr3.intron <- preparingSEForIntronicUsageTesting(se.exprs.utr3.intron, hg19.annotation)

  ##Get the relevant data and perform the analysis
  se.exprs.utr3.intron <- se.exprs.utr3.intron[which(rowRanges(se.exprs.utr3.intron)$cleavage.unit %in% c("utr3", "intronic"))]

  colData(se.exprs.utr3.intron)$condition <- df.design[colnames(se.exprs.utr3.intron), ]
  colData(se.exprs.utr3.intron)$prep.by <- 'peggy'
  acs.utr3.intron <- se2ACS(se.exprs.utr3.intron)

  all.expts <- unique(as.character(design(acs.utr3.intron)$condition))
  all.pairs <- combn(all.expts, 2)

  if(is.null(vs)){
    reqd.pairs <- all.pairs
    all.denominator <- NULL
    all.ref <- all.expts
    all.vs <- all.expts
  } else if(length(vs) == 1){
    these.pairs <- apply(all.pairs, 2, function(x) any(x %in% vs))
    reqd.pairs <- as.matrix(all.pairs[ ,these.pairs])
    all.denominator <- vs
    all.ref <- setdiff(all.expts, vs)
    all.vs <- vs
  }

  chisq.df <- 2
  fit.qtile <- 0.55
  pwise.keep <- PairwiseTestForDEU(ecs = acs.utr3.intron, nCores=1, min.tx.read.count=10, 
    all.pairs = reqd.pairs, all.denominator = all.denominator,
    chisq.df = chisq.df, fit.qtile = fit.qtile, adjust.p = TRUE,
    sf.for.all = sf.for.all)
  
  step1.result.file <- paste(comparison.name, "-diff-intronic-usage-all-results.rds", sep = "")
  saveRDS(pwise.keep, file = step1.result.file)
  
  if("chr" %in% colnames(pwise.keep)){
    setnames(pwise.keep, old = "chr", new = "seqnames")
  }
  
  pwise <- alignAssayData(pwise.keep, se.exprs.utr3.intron, 'usage')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'LI')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'DUI')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'count')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'tpm')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'utr3.tpm')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'transcription.tpm')
  dt.gene.info <- addGeneInfo(se.atlas, annotation.file.location = NULL)
  dt.gene.info$geneID <- NULL
  
  at <- as(rowRanges(se.atlas), "data.frame")
  at$geneID <- NULL
  at$exonID <- NULL
  xx <- merge(at, dt.gene.info, by = c("seqnames", "strand", "entrez.id"), 
    all.x = T, suffixes = c("", ""))
  xx <- xx[ , unique(colnames(xx))]
  yy <- merge(pwise, xx, by = c("seqnames", "strand", "start", "end"), all.x = T,
    suffixes = c("", ""))
  yy <- yy[ , unique(colnames(yy)), with = FALSE]
  
  pwise <- yy
  setnames(pwise, old = c("tpm", "utr3.tpm"), 
    new = c('events.tpm.qq', 'utr3.tpm.qq'))
  
  ##Since the focus is on intronic events so select those
  pwise <- subset(pwise, exon.anno %in% c("intron", "intron.utr3"))
  
  all.pwise <- mclapply(all.ref, function(x) {
      results <- lapply(all.vs, function(y) {
        print(x)
        print(y)
        if(x %in% pwise$sample && (y %in% pwise$sample1 || y %in% pwise$sample2) && x!=y){
          e <- extractPairwiseApA(pwise, reference = x, vs = y, readjust.p=FALSE)
        }
      })
      rbindlist(results)
  }, mc.cores=24L)
  all.pwise <- rbindlist(all.pwise)
  all.pwise <- as.data.frame(all.pwise)
  all.pwise$diff.usage <- abs(all.pwise$ref.usage - all.pwise$vs.usage)
  
  all.pwise$ref.is <- ifelse(all.pwise$ref.B.raw < 0, 
    "longer", "shorter")

  step2.result.file <- paste(comparison.name, "-diff-intronic-usage-final-results.rds", sep = "") 
  saveRDS(all.pwise, file = step2.result.file)
}


###########################################
callDifferentialUsedUTR3Events <- function(se.atlas, enames, df.design = NULL, vs = NULL,
  comparison.name,
  peak.tpm.cutoff = 3L, utr3.tpm.cutoff = 5.5, percent.exprs = 0.75, sf.for.all = T,
  nCores = 1){

  ##Subset the experiments to look at
  se.prep <- se.atlas[ ,enames]

  ##The gene should be expressed in all samples
  se.prep <- se.prep[rowRanges(se.prep)$exon.anno %in% c("utr3", "utr3*")]
  se.prep <- se.prep[apply(assays(se.prep)$utr3.tpm > utr3.tpm.cutoff, 1, all)]
  
  
  ##The peak should be expressed in atleast in 75% samples of any condition
  all.cond <- as.character(unique(colData(se.prep)$condition))
  m.exprs.logical <- sapply(all.cond, function(x) {
    pos <- which(colData(se.prep)$condition == x)
    m.tpm <- as.matrix(assays(se.prep)$tpm[ ,pos])
    l <- apply(m.tpm, 1, function(idx) {
      length(which(idx > peak.tpm.cutoff)) >= (percent.exprs*length(idx))
      })
    return(l)
    })
  exprs.utr3 <- apply(m.exprs.logical, 1, function(x) any(x))
  se.exprs.utr3 <- se.prep[exprs.utr3]

  ##Get the multi-usage ids and flag the two peaks that should be testeed
  filter.for.these.ids <- names(which(table(rowRanges(se.exprs.utr3)$entrez.id) > 1))
  se.multi.utr3 <- se.exprs.utr3[rowRanges(se.exprs.utr3)$entrez.id %in% filter.for.these.ids]
  se.multi.utr3 <- addTestOfUsage(se.multi.utr3)
  ##genes that go into analysis
  length(unique(filter.for.these.ids))

  ##Prepare for differential usage testing
  colData(se.multi.utr3)$condition <- df.design[colnames(se.multi.utr3), ]
  colData(se.multi.utr3)$prep.by <- 'peggy'
  
  acs.utr3 <- se2ACS(se.multi.utr3, 
    additional.cols = "diff.usage.test.candidate")

  all.expts <- unique(as.character(design(acs.utr3)$condition))
  all.pairs <- combn(all.expts, 2)

  if(is.null(vs)){
    reqd.pairs <- all.pairs
    all.denominator <- NULL
    all.ref <- all.expts
    all.vs <- all.expts
  } else if(length(vs) == 1){
    these.pairs <- apply(all.pairs, 2, function(x) any(x %in% vs))
    reqd.pairs <- as.matrix(all.pairs[ ,these.pairs])
    all.denominator <- vs
    all.ref <- setdiff(all.expts, vs)
    all.vs <- vs
  }
  
  chisq.df <- 2
  fit.qtile <- 0.55
  pwise.keep <- PairwiseTestForDEU(acs.utr3, nCores=nCores, min.tx.read.count=10, 
    all.pairs = reqd.pairs, all.denominator = all.denominator,
    chisq.df = chisq.df, fit.qtile = fit.qtile, adjust.p = FALSE, 
    test.peaks = c("proximal.test.peak", "distal.test.peak"), 
    test.peak.col = "diff.usage.test.candidate", sf.for.all = sf.for.all)

  if("chr" %in% colnames(pwise.keep)){
    setnames(pwise.keep, old = "chr", new = "seqnames")
  }

  
  #############
  pwise <- alignAssayData(pwise.keep, se.multi.utr3, 'usage')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'LI')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'DUI')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'count')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'tpm')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'utr3.tpm')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'transcription.tpm')
  dt.gene.info <- addGeneInfo(se.atlas, annotation.file.location = NULL)
  pwise <- annotateCleavageSites(pwise, se.atlas, dt.gene.info)

  setnames(pwise, old = c("tpm", "utr3.tpm"), 
    new = c('events.tpm.qq', 'utr3.tpm.qq'))

  step1.result.file <- paste(comparison.name, "diff-utr3-usage-all-results.rds", sep = "-")
  saveRDS(pwise, file = step1.result.file)
  
  ###############
  all.pwise <- mclapply(all.ref, function(x) {
      results <- lapply(all.vs, function(y) {
        print(x)
        print(y)
        if(x %in% pwise$sample && (y %in% pwise$sample1 || y %in% pwise$sample2) && x!=y){
          e <- extractPairwiseApA(pwise, reference = x, vs = y, readjust.p=FALSE)
        }
      })
      rbindlist(results)
  }, mc.cores=24L)
  all.pwise <- rbindlist(all.pwise)
  all.pwise <- as.data.frame(all.pwise)

  step2.result.file <- paste(comparison.name, "diff-utr3-usage-all-pwise-results.rds", sep = "-")
  saveRDS(all.pwise, file = step2.result.file)
  ######################################
  all.genewise <- genewiseApAForTwoPeaks(all.pwise)

  step3.result.file <- paste(comparison.name, "diff-utr3-usage-all-genewise-results.rds", sep = "-")
  saveRDS(all.genewise, file = step3.result.file)
}

preparingSEForIntronicUsageTesting <- function(se, hg19.annotation){

  se.intron <- se[which(rowRanges(se)$cleavage.unit == "intronic")]
  se.utr3 <- se[which(rowRanges(se)$cleavage.unit == "utr3")]

  se.clustered.utr3 <- clusterPeaksByAnnotation(se.utr3, ann.col="cleavage.unit", 
    cluster.unit = "utr3", requantify=TRUE, 
    requantify.these = c("count", "tpm", "usage"))
  se.clustered.utr3 <- annotateSummarizedExperimentClusterUTR3(se.clustered.utr3, 
    hg19.annotation)
  rowRanges(se.clustered.utr3)$exon.anno <- "utr3"
  
  common.cols <- intersect(colnames(mcols(se.clustered.utr3)), 
    colnames(mcols(se.intron)))

  rowRanges(se.intron) <- rowRanges(se.intron)[ ,common.cols]
  rowRanges(se.clustered.utr3) <- rowRanges(se.clustered.utr3)[ ,common.cols]

  se <- rbind(se.intron, se.clustered.utr3)
  se.prep <- prepSE4ApaTesting(se)

  rowRanges(se.prep)$new.group <-  "genic.cleavage"
  rowRanges(se.prep)$new.group <- ifelse(rowRanges(se.prep)$exon.anno %in% c("intron", "intron.utr3"),
    paste("intron", rowRanges(se.prep)$unit.position, sep = ""), rowRanges(se.prep)$new.group)
  
  rowRanges(se.prep)$new.group <-  ifelse((rowRanges(se.prep)$exon.anno %in% 
    c("utr3", "utr3*")) | (rowRanges(se.prep)$exon.anno == "intron" & rowRanges(se.prep)$flanking.up.unit == "utr3"),
    "utr3.cleavage", rowRanges(se.prep)$new.group)
  
  rowRanges(se.prep)$cleavage.unit <- "genic"
  rowRanges(se.prep)$cleavage.unit <- ifelse(rowRanges(se.prep)$new.group == "utr3.cleavage", "utr3", 
    rowRanges(se.prep)$cleavage.unit)
  rowRanges(se.prep)$cleavage.unit <- ifelse(grepl("intron", rowRanges(se.prep)$new.group), "intronic", 
    rowRanges(se.prep)$cleavage.unit)

  ##As we want to test every intronic peak vs full length utr3
  ##so I am repeating the 3'utr with different geneIDs, one for every
  ##intronic peak using the exon id of each intrnic peak
  dt <- as(rowRanges(se.prep), "data.table")
  dt$row.idx <- 1:nrow(dt)
  yy <- dt[ ,{
      entrez.id <- .BY[1]
      n.intron <- sum(cleavage.unit == "intronic")
      idx.utr3 <- which(cleavage.unit == "utr3")
      idx.intronic <- which(cleavage.unit == "intronic")
      intronic.exon.id <- exonID[idx.intronic]
      gene.id <- paste(entrez.id, intronic.exon.id, sep =".")
      x <- {}
      for(i in 1:n.intron){
        x <- rbind(x, .SD[idx.utr3, ])
      }
      x$geneID <- gene.id
      y <- .SD[idx.intronic, ]
      y$geneID <- gene.id
      rbind(x, y)
    }, by = "entrez.id"]

  ##using the row.ids to recreate the summarized experiment object
  new.se <- se.prep[yy$row.idx, ]
  rowRanges(new.se) <- as(yy, "GRanges")
  rowRanges(new.se)$row.idx <- NULL

  new.se
  
}

#####################
#This clusters peaks by the annotation unit
#####################
clusterPeaksByAnnotation <- function(x, ann.col="cleavage.unit", 
  cluster.unit = "utr3", requantify=TRUE, 
  requantify.these = c("count", "tpm", "usage")) {

  stopifnot(inherits(x, 'SummarizedExperiment'))
  stopifnot(ann.col %in% colnames(rowRanges(x)))

  xx <- rowRanges(x)
  mcols(xx) <- NULL

  y <- as(rowRanges(x), "data.table")
  mask <- y[ ,{
    if(.BY[[2]] == cluster.unit){
      y.start <- min(start)
      y.end <- max(end)
      y.strand <- strand[1L]
      y.seq <- seqnames[1L]
      yy <- list(seqnames = y.seq, start = y.start, end = y.end, strand = y.strand)
    } else {
      yy <- .SD[ ,list(seqnames, start, end, strand)]
    }
     yy
    }, by = c("entrez.id", "cleavage.unit")]
  mask <- as(mask, "GRanges")

  mm <- as.matrix(findOverlaps(xx, mask))
  if (any(duplicated(mm[,1]))) {
    stop("There shouldn't be any duplicate queryHits")
  }

  ## Reasign fenceposts from original peak definitions
  xx$n.signal <- rowRanges(x)$n.signal
  xx <- as(xx, 'data.table')

  mask.idx <- integer(nrow(mm))
  mask.idx[mm[,1]] <- mm[,2]
  xx[, mask.idx := mask.idx]
  setkeyv(xx, 'mask.idx')

  bounds <- xx[, {
    list(seqnames=seqnames[1L], start=min(start), end=max(end),
         strand=strand[1L], npeaks=.N, n.signal=sum(n.signal))
  }, by="mask.idx"]

  out <- as(bounds[, mask.idx := NULL], 'GRanges')
  if (requantify) {
    out <- newSEfromRegions(out, x, change = requantify.these)
  }
  out
}

##############
#This annotation uses assign = "fix" to make sure that the clustered utr3 peaks
#get the right annotation
##############
annotateSummarizedExperimentClusterUTR3 <- function(x, annotations, nuclear=TRUE, assign = "fix") {
  stopifnot(inherits(x, "SummarizedExperiment"))
  stopifnot(inherits(annotations, "GenomicRanges"))

  regions <- rowRanges(x)
  xmeta <- mcols(regions)
  mcols(regions) <- NULL

  anno.cols <- names(mcols(annotations))

  if (nuclear) {
    for (wut in anno.cols) {
      if (wut %in% names(xmeta)) {
        xmeta[[wut]] <- NULL
      }
    }
  }

  rowRanges(regions) <- xmeta
  suppressWarnings(regions <- annotateReads(regions, annotations, assign = assign))
  rowRanges(x) <- regions
  x
}



###########################################
callDifferentialUsedIntronicEvents <- function(se.atlas, enames, df.design = NULL, vs = NULL,
  comparison.name, utr3.exprs = TRUE,
  peak.tpm.cutoff = 5L, utr3.tpm.cutoff = 5.5, percent.exprs = 0.75, sf.for.all = T){

    
  ##Subset the experiments to look at
  se.prep <- se.atlas[ ,enames]
  rowRanges(se.prep)$idx <- 1:length(rowRanges(se.prep))
  
  se.intron <- se.prep[rowRanges(se.prep)$cleavage.unit == "intronic"]
    
  ##The intronic peak should be expressed in atleast one condition with the  given cutoff
  all.cond <- as.character(unique(colData(se.intron)$condition))

  m.exprs.logical <- sapply(all.cond, function(x) {
    pos <- which(colData(se.intron)$condition == x)
    m.tpm <- as.matrix(assays(se.intron)$tpm[ ,pos])
    l <- apply(m.tpm, 1, function(idx) {
      length(which(idx > peak.tpm.cutoff)) >= (percent.exprs*length(idx))
      })
    return(l)
    })
  exprs.intron <- apply(m.exprs.logical, 1, function(x) any(x))
  se.exprs.intron <- se.intron[exprs.intron]

  if(utr3.exprs){
    ######
    ##We can test only the intronic peaks that have atleast intron or utr3 expression 
    ##in all conditions
    m.exprs.all.logical <- sapply(all.cond, function(x) {
      pos <- which(colData(se.exprs.intron)$condition == x)
      m.intron.tpm <- as.matrix(assays(se.exprs.intron)$tpm[ ,pos] > peak.tpm.cutoff)
      m.utr3.tpm <- as.matrix(assays(se.exprs.intron)$utr3.tpm[ ,pos] > utr3.tpm.cutoff)
      exp <- apply(m.intron.tpm | m.utr3.tpm, 1, function(idx){
        length(which(idx)) >= (percent.exprs*length(idx))
        })
      })
    exprs.intron.all <- apply(m.exprs.all.logical, 1, function(x) all(x))
    se.intron.to.test <- se.exprs.intron[exprs.intron.all]
      
  } else {

    se.intron.to.test <- se.exprs.intron

  }

  ###Number of introns to be tested
  print(length(rowRanges(se.intron.to.test)))
  
  ##Now select the utr3 and intron eventss for the genes that have intronic peaks
  se.exprs.utr3.intron <- se.prep[(rowRanges(se.prep)$cleavage.unit == "utr3" &
    rowRanges(se.prep)$entrez.id %in% unique(rowRanges(se.intron.to.test)$entrez.id)) |
    (rowRanges(se.prep)$idx %in% rowRanges(se.intron.to.test)$idx)]

  
  ##Since evry intronic peak is tested against its full length so we
  ##need to change the gene id for the full length
  hg19.annotation.file <- "/Users/singhi1/ifs/annotation_data/hg19/hg19_annotation_pfam.rds"
  hg19.annotation <- readRDS(hg19.annotation.file)
  se.exprs.utr3.intron <- preparingSEForIntronicUsageTesting(se.exprs.utr3.intron, hg19.annotation)

  ##Get the relevant data and perform the analysis
  se.exprs.utr3.intron <- se.exprs.utr3.intron[which(rowRanges(se.exprs.utr3.intron)$cleavage.unit %in% c("utr3", "intronic"))]

  colData(se.exprs.utr3.intron)$condition <- df.design[colnames(se.exprs.utr3.intron), ]
  colData(se.exprs.utr3.intron)$prep.by <- 'peggy'
  acs.utr3.intron <- se2ACS(se.exprs.utr3.intron)

  all.expts <- unique(as.character(design(acs.utr3.intron)$condition))
  all.pairs <- combn(all.expts, 2)

  if(is.null(vs)){
    reqd.pairs <- all.pairs
    all.denominator <- NULL
    all.ref <- all.expts
    all.vs <- all.expts
  } else if(length(vs) == 1){
    these.pairs <- apply(all.pairs, 2, function(x) any(x %in% vs))
    reqd.pairs <- as.matrix(all.pairs[ ,these.pairs])
    all.denominator <- vs
    all.ref <- setdiff(all.expts, vs)
    all.vs <- vs
  }

  chisq.df <- 2
  fit.qtile <- 0.55
  pwise.keep <- PairwiseTestForDEU(ecs = acs.utr3.intron, nCores=1, min.tx.read.count=10, 
    all.pairs = reqd.pairs, all.denominator = all.denominator,
    chisq.df = chisq.df, fit.qtile = fit.qtile, adjust.p = TRUE,
    sf.for.all = sf.for.all)
  
  step1.result.file <- paste(comparison.name, "-diff-intronic-usage-all-results.rds", sep = "")
  saveRDS(pwise.keep, file = step1.result.file)
  
  if("chr" %in% colnames(pwise.keep)){
    setnames(pwise.keep, old = "chr", new = "seqnames")
  }
  
  pwise <- alignAssayData(pwise.keep, se.exprs.utr3.intron, 'usage')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'LI')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'DUI')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'count')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'tpm')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'utr3.tpm')
  pwise <- alignAssayData(pwise, se.exprs.utr3.intron, 'transcription.tpm')
  dt.gene.info <- addGeneInfo(se.atlas, annotation.file.location = NULL)
  dt.gene.info$geneID <- NULL
  
  at <- as(rowRanges(se.atlas), "data.frame")
  at$geneID <- NULL
  at$exonID <- NULL
  xx <- merge(at, dt.gene.info, by = c("seqnames", "strand", "entrez.id"), 
    all.x = T, suffixes = c("", ""))
  xx <- xx[ , unique(colnames(xx))]
  yy <- merge(pwise, xx, by = c("seqnames", "strand", "start", "end"), all.x = T,
    suffixes = c("", ""))
  yy <- yy[ , unique(colnames(yy)), with = FALSE]
  
  pwise <- yy
  setnames(pwise, old = c("tpm", "utr3.tpm"), 
    new = c('events.tpm.qq', 'utr3.tpm.qq'))
  
  ##Since the focus is on intronic events so select those
  pwise <- subset(pwise, exon.anno %in% c("intron", "intron.utr3"))
  
  all.pwise <- mclapply(all.ref, function(x) {
      results <- lapply(all.vs, function(y) {
        print(x)
        print(y)
        if(x %in% pwise$sample && (y %in% pwise$sample1 || y %in% pwise$sample2) && x!=y){
          e <- extractPairwiseApA(pwise, reference = x, vs = y, readjust.p=FALSE)
        }
      })
      rbindlist(results)
  }, mc.cores=24L)
  all.pwise <- rbindlist(all.pwise)
  all.pwise <- as.data.frame(all.pwise)
  all.pwise$diff.usage <- abs(all.pwise$ref.usage - all.pwise$vs.usage)
  
  all.pwise$ref.is <- ifelse(all.pwise$ref.B.raw < 0, 
    "longer", "shorter")

  step2.result.file <- paste(comparison.name, "-diff-intronic-usage-final-results.rds", sep = "") 
  saveRDS(all.pwise, file = step2.result.file)
}


###########################################
callDifferentialUsedUTR3Events <- function(se.atlas, enames, df.design = NULL, vs = NULL,
  comparison.name,
  peak.tpm.cutoff = 3L, utr3.tpm.cutoff = 5.5, percent.exprs = 0.75, sf.for.all = T,
  nCores = 1){

  ##Subset the experiments to look at
  se.prep <- se.atlas[ ,enames]

  ##The gene should be expressed in all samples
  se.prep <- se.prep[rowRanges(se.prep)$exon.anno %in% c("utr3", "utr3*")]
  se.prep <- se.prep[apply(assays(se.prep)$utr3.tpm > utr3.tpm.cutoff, 1, all)]
  
  
  ##The peak should be expressed in atleast in 75% samples of any condition
  all.cond <- as.character(unique(colData(se.prep)$condition))
  m.exprs.logical <- sapply(all.cond, function(x) {
    pos <- which(colData(se.prep)$condition == x)
    m.tpm <- as.matrix(assays(se.prep)$tpm[ ,pos])
    l <- apply(m.tpm, 1, function(idx) {
      length(which(idx > peak.tpm.cutoff)) >= (percent.exprs*length(idx))
      })
    return(l)
    })
  exprs.utr3 <- apply(m.exprs.logical, 1, function(x) any(x))
  se.exprs.utr3 <- se.prep[exprs.utr3]

  ##Get the multi-usage ids and flag the two peaks that should be testeed
  filter.for.these.ids <- names(which(table(rowRanges(se.exprs.utr3)$entrez.id) > 1))
  se.multi.utr3 <- se.exprs.utr3[rowRanges(se.exprs.utr3)$entrez.id %in% filter.for.these.ids]
  se.multi.utr3 <- addTestOfUsage(se.multi.utr3)
  ##genes that go into analysis
  length(unique(filter.for.these.ids))

  ##Prepare for differential usage testing
  colData(se.multi.utr3)$condition <- df.design[colnames(se.multi.utr3), ]
  colData(se.multi.utr3)$prep.by <- 'peggy'
  
  acs.utr3 <- se2ACS(se.multi.utr3, 
    additional.cols = "diff.usage.test.candidate")

  all.expts <- unique(as.character(design(acs.utr3)$condition))
  all.pairs <- combn(all.expts, 2)

  if(is.null(vs)){
    reqd.pairs <- all.pairs
    all.denominator <- NULL
    all.ref <- all.expts
    all.vs <- all.expts
  } else if(length(vs) == 1){
    these.pairs <- apply(all.pairs, 2, function(x) any(x %in% vs))
    reqd.pairs <- as.matrix(all.pairs[ ,these.pairs])
    all.denominator <- vs
    all.ref <- setdiff(all.expts, vs)
    all.vs <- vs
  }
  
  chisq.df <- 2
  fit.qtile <- 0.55
  pwise.keep <- PairwiseTestForDEU(acs.utr3, nCores=nCores, min.tx.read.count=10, 
    all.pairs = reqd.pairs, all.denominator = all.denominator,
    chisq.df = chisq.df, fit.qtile = fit.qtile, adjust.p = FALSE, 
    test.peaks = c("proximal.test.peak", "distal.test.peak"), 
    test.peak.col = "diff.usage.test.candidate", sf.for.all = sf.for.all)

  if("chr" %in% colnames(pwise.keep)){
    setnames(pwise.keep, old = "chr", new = "seqnames")
  }

  
  #############
  pwise <- alignAssayData(pwise.keep, se.multi.utr3, 'usage')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'LI')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'DUI')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'count')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'tpm')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'utr3.tpm')
  pwise <- alignAssayData(pwise, se.multi.utr3, 'transcription.tpm')
  dt.gene.info <- addGeneInfo(se.atlas, annotation.file.location = NULL)
  pwise <- annotateCleavageSites(pwise, se.atlas, dt.gene.info)

  setnames(pwise, old = c("tpm", "utr3.tpm"), 
    new = c('events.tpm.qq', 'utr3.tpm.qq'))

  step1.result.file <- paste(comparison.name, "diff-utr3-usage-all-results.rds", sep = "-")
  saveRDS(pwise, file = step1.result.file)
  
  ###############
  all.pwise <- mclapply(all.ref, function(x) {
      results <- lapply(all.vs, function(y) {
        print(x)
        print(y)
        if(x %in% pwise$sample && (y %in% pwise$sample1 || y %in% pwise$sample2) && x!=y){
          e <- extractPairwiseApA(pwise, reference = x, vs = y, readjust.p=FALSE)
        }
      })
      rbindlist(results)
  }, mc.cores=24L)
  all.pwise <- rbindlist(all.pwise)
  all.pwise <- as.data.frame(all.pwise)

  step2.result.file <- paste(comparison.name, "diff-utr3-usage-all-pwise-results.rds", sep = "-")
  saveRDS(all.pwise, file = step2.result.file)
  ######################################
  all.genewise <- genewiseApAForTwoPeaks(all.pwise)

  step3.result.file <- paste(comparison.name, "diff-utr3-usage-all-genewise-results.rds", sep = "-")
  saveRDS(all.genewise, file = step3.result.file)
}

assignCleavageUnit <- function(se){
  ## Add another column which clubs the utr3, utr3*
  ## as one cleavge event
  rowRanges(se)$new.group <-  "genic.cleavage"
  rowRanges(se)$new.group <- ifelse(rowRanges(se)$exon.anno %in% c("intron", "intron.utr3"),
    paste("intron", rowRanges(se)$unit.position, sep = ""), rowRanges(se)$new.group)
  
  rowRanges(se)$new.group <-  ifelse((rowRanges(se)$exon.anno %in% 
    c("utr3", "utr3*")) | (rowRanges(se)$exon.anno == "intron" & rowRanges(se)$flanking.up.unit == "utr3"),
    "utr3.cleavage", rowRanges(se)$new.group)
  
  rowRanges(se)$cleavage.unit <- "genic"
  rowRanges(se)$cleavage.unit <- ifelse(rowRanges(se)$new.group == "utr3.cleavage", "utr3", 
    rowRanges(se)$cleavage.unit)
  rowRanges(se)$cleavage.unit <- ifelse(grepl("intron", rowRanges(se)$new.group), "intronic", 
    rowRanges(se)$cleavage.unit)

  rowRanges(se)$new.group <- NULL

  se
}


## for Rle or vector
##For swapping the strand
swapStrand <- function(x, ...) {
  strands <- as.character(strand(x))
  consider <- strands != '*'
  to.neg <- strands == '+' & consider
  to.pos <- strands == '-' & consider
  if (sum(to.neg) > 0) {
    strands[to.neg] <- '-'
  }
  if (sum(to.pos) > 0) {
    strands[to.pos] <- '+'
  }
  strand(x) <- strands
  x
}

calculateCI <- function(x, at='proximal', fraction.of='total',
                        group.by=c('seqnames', 'strand', 'entrez.id'),
                        groups=list(
                          utr3=c('utr3', 'utr3*'),
                          intron.utr3='intron.utr3',
                          genic=c('utr5', 'cds')),
                        groups.col='exon.anno') {
  stopifnot(is(x, "SummarizedExperiment"))
  fraction.of <- match.arg(fraction.of)
  if (is.character(at)) {
    at <- match.arg(at, c('proximal', 'distal'))
  } else {
    stop("Numeric `at` not supported")
  }

  if (is.null(groups) || length(groups) == 0) {
    groups <- list(all=unique(mcols(rowRanges(x))[[groups.col]]))
  } else {
    stopifnot(all(unlist(groups) %in% mcols(rowRanges(x))[[groups.col]]))
  }

  things <- assays(x)
  if (!"usage" %in% names(things)) {
    x <- calculateUsage(x)
  }

  u <- cbind(as(rowRanges(x), 'data.table'), data.table(assay(x, "usage")))
  u[, .idx. := seq(nrow(u))]
  setkeyv(u, c(group.by, 'start'))

  if (!'group' %in% colnames(u)) {
    warning("All events in the tx-unit are considered together -- call `addProcessingGroup`")
    u[, group := 'unk']
  }

  expt.names <- colnames(x)

  ans <- matrix(NA_real_, nrow=nrow(x), ncol=ncol(x))
  colnames(ans) <- colnames(x)

  for (pgroup in unique(u$group)) {
    usub <- subset(u, group == pgroup)
    if (at == 'proximal' || at == 1) {
      ci <- usub[, {
        idx <-  if (strand == '+') 1L else .N
        cbind(.SD[rep(idx, .N), expt.names, with=FALSE], .SD[, list(.idx.)])
      }, by=group.by]
    } else if (at == "distal") {
      ci <- usub[, {
        idx <-  if (strand == '+') .N else 1L
        cbind(.SD[rep(idx, .N), expt.names, with=FALSE], .SD[, list(.idx.)])
      }, by=group.by]
    }
    ans[ci$.idx., expt.names] <- as.matrix(ci[, expt.names, with=FALSE])
  }

  wut <- if (at == 'proximal') "SUI" else "DUI"
  assay(x, wut) <- ans
  x
}


##' Returns TRUE for all rows in `x` that have non-NA columns
complete.columns <- function(x, columns, complete.test=Negate(is.na)) {
  if (is.character(columns)) {
    stopifnot(all(columns %in% colnames(x)))
  }
  if (is.numeric(columns)) {
    stopifnot(all(columns >= 1) & all(columns <= ncol(x)))
  }
  is.complete <- rep(TRUE, nrow(x))
  for (i in columns) {
    is.complete <- is.complete & complete.test(x[, i, with = F])
  }
  is.complete
}

setAs("GRanges", "data.table", function(from) {
  if (length(from) == 0L) {
    return(data.table())
  }
  as.data.table(as(from, 'data.frame'))
})

setAs("GRanges", "data.frame", function(from) {
  x <- lapply(as.data.frame(from), function(xx) {
    if (is.factor(xx)) {
      xx <- as.character(xx)
    }
    xx
  })
  as.data.frame(x, stringsAsFactors=FALSE)
})

setAs("IRanges", "data.table", function(from) {
  if (length(from) == 0L) {
    return(data.table())
  }
  as.data.table(cbind(as.data.frame(from), as.data.frame(values(from))))
})

setAs("data.table", "GRanges", function(from) {
  as(as.data.frame(from), "GRanges")
})

setAs("data.frame", "GRanges", function(from) {
  if (nrow(from) == 0L || all(is.na(from))) {
    return(GRanges())
  }
  if (!'seqnames' %in% colnames(from)) {
    stop("seqnames required")
  }
  gr.meta.take <- colnames(from) %in% c('seqnames', 'strand')
  gr.meta <- from[, gr.meta.take, drop=FALSE]
  from <- from[, !gr.meta.take, drop=FALSE]

  .ranges <- as(from, 'IRanges')
  DF <- elementMetadata(.ranges)
  elementMetadata(.ranges) <- NULL

  .strand <- if ('strand' %in% colnames(gr.meta)) gr.meta$strand else '*'
  gr <- GRanges(seqnames=gr.meta$seqnames, ranges=.ranges, strand=.strand)
  values(gr) <- DF
  gr
})

setAs("data.table", "IRanges", function(from) {
  as(as.data.frame(from), "IRanges")
})

setAs("data.frame", "IRanges", function(from) {
  if (nrow(from) == 0L || all(is.na(from))) {
    return(IRanges())
  }
  need <- c('start', 'end', 'width')
  have <- colnames(from)[colnames(from) %in% need]
  if (length(have) < 2) {
    stop("Need any two of 'start', 'end', or 'width'")
  }
  ## prefer start/end
  if (all(c('start', 'end') %in% have)) {
    iranges <- IRanges(from$start, from$end)
  } else {
    iranges <- do.call(IRanges, as.list(from[, have]))
  }

  meta.cols <- setdiff(colnames(from), have)

  if (length(meta.cols) > 0) {
    DF <- DataFrame(from[, meta.cols, drop=FALSE])
    colnames(DF) <- meta.cols
    values(iranges) <- DF
  }

  iranges
})

GroupFilterOnUsageAndTPM <- function(se, enames = NULL,
  vec.group = c("utr3", "utr3*"), unit = NULL,
  tpm.cutoff = 3, ui.cutoff = 0.1, iqr.filter = FALSE, iqr.cutoff = 3){

  if(is.null(enames)){
      enames <- colnames(se)
    }


  dt.tpm <- getAssay(se[ ,enames], "tpm")
  setnames(dt.tpm, old = enames, new = paste(enames, "tpm", sep = "."))
  dt.tpm <- dt.tpm[ ,paste(enames, "tpm", sep = "."), with = FALSE]

  dt.ui <- getAssay(se[ ,enames], "usage")
  setnames(dt.ui, old = enames, new = paste(enames, "ui", sep = "."))
  dt.ui <- dt.ui[ ,paste(enames, "ui", sep = "."), with = FALSE]

  m.ui.logical <- apply(dt.ui, 2, function(x) x >= ui.cutoff)
  m.tpm.logical <- apply(dt.tpm, 2, function(x) x >= tpm.cutoff)

  if(iqr.filter){
    dt.iqr <- getAssay(se[ ,enames], "start.iqr")
    setnames(dt.iqr, old = enames, new = paste(enames, "iqr", sep = "."))
    dt.iqr <- dt.iqr[ ,paste(enames, "iqr", sep = "."), with = FALSE]
    m.iqr.logical <- apply(dt.iqr, 2, function(x) x >= iqr.cutoff)

    vec.cutoff.filter <- apply(m.ui.logical & m.tpm.logical & m.iqr.logical, 1, any)
  } else {
    vec.cutoff.filter <- apply(m.ui.logical & m.tpm.logical, 1, any)
  }

  ## The usage and tpm filter is only applicable to group
  ## So select rows that qualify all the 3 filter

  if(is.null(unit)){
    vec.is.group <- values(se)$exon.anno %in% vec.group
    } else {
      vec.is.group <- values(se)$cleavage.unit %in% unit
    }
  
  vec.group.filter <- apply(cbind(vec.cutoff.filter, vec.is.group), 1, all)
  vec.is.not.group <- !vec.is.group
  vec.logical <- apply(cbind(vec.group.filter, vec.is.not.group), 1, any )
  se[which(vec.logical)]
}

## Gets the data crammed in the assays() of the atlas SE

#' Gets a names experiment statistic associated to each event
#'
#' @param atlas the atlas
#' @param assay.idx The name (or integer idx) of the assay you want to get
#' @param group only used as a convenience to subset out only values for
#' events that are assigned to the particular `group`. We are using `utr3` for
#' this analysis.
getAssay <- function(x, assay.idx, group=NULL, meta.cols=NULL,
                     combine.by=c('none', 'mean', 'median', 'sum'), ...) {
  stopifnot(is(x, "SummarizedExperiment"))
  if (is.character(combine.by)) {
    combine.by <- match.arg(combine.by)
    combine.by <- switch(combine.by, mean=rowMeans, median=rowMedians,
                         sum=rowSums, none=NA)
  } else {
    stopifnot(is.function(combine.by))
  }

  if (is.numeric(assay.idx)) {
    stopifnot(as.integer(assay.idx) == assay.idx)
    stopifnot(assay.idx >= 1 && assay.idx <= length(assays(x)))
  }
  if (is.character(assay.idx)) {
    stopifnot(assay.idx %in% names(assays(x)))
  }
  if (!is.null(group)) {
    stopifnot(group %in% mcols(rowData(x))[['group']])
  }

  meta.cols <- c('seqnames', 'strand', 'start', 'end', colnames(values(x)))

  dat <- assay(x, assay.idx)
  if (is.function(combine.by)) {
    adesign <- defactor(as.data.frame(colData(x)))
    xref <- split(rownames(adesign), adesign$condition)
    dat <- sapply(xref, function(col.names) {
      combine.by(dat[, col.names, drop=FALSE], na.rm=TRUE)
    })
  }
  dat[!is.finite(dat)] <- NA_real_

  ans <- cbind(as.data.frame(rowRanges(x))[, meta.cols], dat)
  ans <- as.data.table(ans)
  #commenting this in local
  #setnames(ans, 'peak.id', 'exonID')

  if (!is.null(group)) {
    ans <- ans[ans$group == group,]
  }

  ans[, `:=`(seqnames=as.character(seqnames), strand=as.character(strand))]
  ans
}
