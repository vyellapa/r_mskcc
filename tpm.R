##' Calculates the TPM for the number of reads \code{x} in a sample with a total
##' number of reads \code{y}
setGeneric("TPM", signature=c('x', 'normer'),
function(x, normer, ...) {
  standardGeneric("TPM")
})

setMethod("TPM", c(x="numeric", normer="numeric"),
function(x, normer, ...) {
  (1e6 / normer) * x
})

setMethod("TPM", c(x="numeric", normer="TagSeq"),
function(x, normer, ...) {
    TPM(x, getReadCount(normer))
})

setMethod("TPM", c("RangedSummarizedExperiment", "numeric"),
function(x, normer, do.group=!missing(group.cols),
         group.cols=c('seqnames', 'strand', 'entrez.id'),
         meta.cols=c(group.cols, "symbol"),
         assay.idx=1, na.rm.keys=TRUE, ...) {
  expt.names <- colnames(x)
  if (is.null(expt.names) || length(expt.names) != ncol(x)) {
    stop("Experiment names are required for the SummarizedExperiment")
  }
  if (!all(expt.names %in% names(normer))) {
    stop("Library sizes are not present for all expt.names")
  }

  meta <- data.table( as(rowRanges(x), 'data.frame') )
  if (!all(meta.cols %in% names(meta))) {
    stop("Requested meta.cols not in metadata")
  }

  meta <- meta[, meta.cols, with=FALSE]
  dt <- cbind(meta, data.table(assay(x, assay.idx)))

  if (isTRUE(do.group)) {
    if (!all(group.cols %in% names(meta))) {
      stop("group.cols not found in metadata")
    }
    if (na.rm.keys) {
      ## Remove any rows that have NA's in the keys
      keep <- Reduce("&", lapply(group.cols, function(kcol) !is.na(dt[[kcol]])))
      n.axe <- sum(!keep)
      if (n.axe) {
        warning("Removing ", n.axe, " rows due to NA in key column", .immediate=TRUE)
        dt <- dt[keep,]
      }
    }
    setkeyv(dt, group.cols)
    more.cols <- setdiff(meta.cols, group.cols)
    dt <- dt[, {
      cnts <- sapply(expt.names, function(wut) sum(.SD[[wut]]), simplify=FALSE)
      more <- sapply(more.cols, function(wut) .SD[[wut]][1], simplify=FALSE)
      c(more, cnts)
    }, by=group.cols]
  }

  tpmf <- getMethod("TPM", c("numeric", "numeric"))

  for (wut in expt.names) {
    dt[[wut]] <- tpmf(dt[[wut]], normer[[wut]])
  }

  as.data.frame(dt)
})

setMethod("TPM", c("RangedSummarizedExperiment", "missing"),
function(x, normer, assay.idx=1, ...) {
  TPM(x, colSums(assay(x, assay.idx)), assay.idx=assay.idx, ...)
})


setMethod("TPM", c("SummarizedExperiment", "numeric"),
function(x, normer, do.group=!missing(group.cols),
         group.cols=c('seqnames', 'strand', 'entrez.id'),
         meta.cols=c(group.cols, "symbol"),
         assay.idx=1, na.rm.keys=TRUE, ...) {
  expt.names <- colnames(x)
  if (is.null(expt.names) || length(expt.names) != ncol(x)) {
    stop("Experiment names are required for the SummarizedExperiment")
  }
  if (!all(expt.names %in% names(normer))) {
    stop("Library sizes are not present for all expt.names")
  }

  meta <- data.table( as(rowData(x), 'data.frame') )
  if (!all(meta.cols %in% names(meta))) {
    stop("Requested meta.cols not in metadata")
  }

  meta <- meta[, meta.cols, with=FALSE]
  dt <- cbind(meta, data.table(assay(x, assay.idx)))

  if (isTRUE(do.group)) {
    if (!all(group.cols %in% names(meta))) {
      stop("group.cols not found in metadata")
    }
    if (na.rm.keys) {
      ## Remove any rows that have NA's in the keys
      keep <- Reduce("&", lapply(group.cols, function(kcol) !is.na(dt[[kcol]])))
      n.axe <- sum(!keep)
      if (n.axe) {
        warning("Removing ", n.axe, " rows due to NA in key column", .immediate=TRUE)
        dt <- dt[keep,]
      }
    }
    setkeyv(dt, group.cols)
    more.cols <- setdiff(meta.cols, group.cols)
    dt <- dt[, {
      cnts <- sapply(expt.names, function(wut) sum(.SD[[wut]]), simplify=FALSE)
      more <- sapply(more.cols, function(wut) .SD[[wut]][1], simplify=FALSE)
      c(more, cnts)
    }, by=group.cols]
  }

  tpmf <- getMethod("TPM", c("numeric", "numeric"))

  for (wut in expt.names) {
    dt[[wut]] <- tpmf(dt[[wut]], normer[[wut]])
  }

  as.data.frame(dt)
})

setMethod("TPM", c("SummarizedExperiment", "missing"),
function(x, normer, assay.idx=1, ...) {
  TPM(x, colSums(assay(x, assay.idx)), assay.idx=assay.idx, ...)
})

##' Calculates the number of read counts needed to get a TPM value of \code{x}
##' from a sample with \code{y}-many reads.
setGeneric("iTPM", signature=c('x', 'normer'),
function(x, normer, ...) {
  standardGeneric("iTPM")
})

setMethod("iTPM", c(x="numeric", normer='numeric'),
function(x, normer, ...) {
  (normer * x) / 1e6
})

setMethod("iTPM", c(x="numeric", normer="TagSeq"),
function(x, normer, ...) {
  iTPM(x, getReadCount(normer, ...))
})



