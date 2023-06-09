##' Calculates percent overlaps of query against subject.
setGeneric("quantifyOverlaps",
function(query, subject, maxgap=0L, minoverlap=1L, ...) {
  standardGeneric("quantifyOverlaps")
})

setMethod("quantifyOverlaps", c(query="Ranges", subject="Ranges"),
function(query, subject, maxgap, minoverlap, ...) {
  o <- findOverlaps(query, subject, maxgap, minoverlap)
  overlapped <- ranges(o, query, subject)
  p.overlap <- width(overlapped) / width(query)[queryHits(o)]
  p.overlap <- ifelse(p.overlap > 1, 1, p.overlap)
  cbind(as.matrix(as.data.frame(o)), p.overlap)
})

setMethod("quantifyOverlaps", c(query="GRanges", subject="GRanges"),
function(query, subject, maxgap, minoverlap, ...) {

  ##############################################################################
  ## TODO: Fix quantifyOverlaps to work on circular chromosomes
  ##       Had a hell of a time chasing a bug down from a read that went
  ##       "circular" around chrM ... fix ends of query to be no longer
  ##       than the end of the chromosome
  ##
  ##       I think this error only ever happened because the reads object
  ##       didn't have any proper values for seqlengths, so the "overhang"
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
  cbind(as.matrix(as.data.frame(o)), p.overlap)
})

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

##' Counts the number of ranges in \code{subject} that land in \code{query}.
##'
##' This function differs from the normal \code{countOverlaps} function by
##' assigning each read in subject to a unique range in query.
##'
##' The default is to assign each read in \code{subject} uniquely by calculating
##' which range in \code{query} it overlaps the most with
##'
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param query \code{IRanges} object
##' @param subject \code{IRanges} object
##' @param assign.by Determine how ranges in subject are assigned to a unique
##' query range
##' @param fix Where to "fix" the resizing of ranges in \code{subject} when
##' \code{assign.by == 'fix'}
##'
##' @return An integer vector indicating how many ranges in \code{subject}
##' overlap with the ranges defined in \code{query}.
countUniqueOverlaps <- function(query, subject, assign.by=c('quantify', 'fix'),
                                maxgap=0L, minoverlap=1L,
                                fix=c('start', 'end', 'center'), ...) {
  assign.by <- match.arg(assign.by)
  if (assign.by == 'fix') {
    fix <- match.arg(fix)
    counts <- countOverlaps(query, resize(subject, width=1, fix=fix))
  } else {
    counts <- integer(length(query))
    qo <- data.table(quantifyOverlaps(query, subject), key='subjectHits')
    if (nrow(qo) == 0L) {
      return(counts)
    }
    qo.unique <- qo[, list(queryHits=queryHits[which.max(p.overlap)]),
                    by='subjectHits']
    key(qo.unique) <- 'queryHits'
    tally <- qo.unique[, list(count=length(query)), by='queryHits']
    counts[tally$query] <- tally$count
  }
  counts
}

###############################################################################
## Utility functions, courtesy of Tim Triche, Jr.
## Example usage:
##   (hESC.H3K4ME1 %s% hESC.P300) %u% (hESC.H3K4ME3 %s% hESC.H3K27ME3)

# intersection
setGeneric('%i%', function(x, y) standardGeneric('%i%'))
setMethod('%i%', c('ANY','ANY'), function(x, y) {
  intersect(x, y)
})

# union
setGeneric('%u%', function(x, y) standardGeneric('%u%'))
setMethod('%u%', c('ANY','ANY'), function(x, y) {
  union(x, y)
})

# setdiff
setGeneric('%d%', function(x, y) standardGeneric('%d%'))
setMethod('%d%', c('ANY','ANY'), function(x, y) {
  setdiff(x, y)
})

# subsetByOverlaps
setGeneric('%s%', function(x, y) standardGeneric('%s%'))
setMethod('%s%', c('GRanges','GRanges'), function(x, y) {
  subsetByOverlaps(x, y)
})

# intersection
setGeneric('%i%', function(x, y) standardGeneric('%i%'))
setMethod('%i%', c('ANY','ANY'), function(x, y) {
  intersect(x, y)
})

# union
setGeneric('%u%', function(x, y) standardGeneric('%u%'))
setMethod('%u%', c('ANY','ANY'), function(x, y) {
  union(x, y)
})

# setdiff
setGeneric('%d%', function(x, y) standardGeneric('%d%'))
setMethod('%d%', c('ANY','ANY'), function(x, y) {
  setdiff(x, y)
})

# subsetByOverlaps
setGeneric('%s%', function(x, y) standardGeneric('%s%'))
setMethod('%s%', c('GRanges','GRanges'), function(x, y) {
  subsetByOverlaps(x, y)
})


test.quantifyOverlaps <- function() {
  ir.1 <- IRanges(c(1, 30, 50), width=10)
  ir.2 <- IRanges(c(2, 25, 80), width=10)

  gr.1 <- GRanges(c('chr1'), ir.1, strand='*')
  gr.2 <- GRanges(c("chr1", 'chr2', 'chr2'), ir.2, strand='*')
}
