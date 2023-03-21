setGeneric("transform", function(`_data`, ...) standardGeneric("transform"),
           useAsDefault=function(`_data`, ...) base::transform)
##' transform function for DataFrame objects that behaves like the one for
##' data.frames
setMethod("transform", "DataFrame", function(`_data`, ...) {
  e <- eval(substitute(list(...)), `_data`, parent.frame())
  tags <- names(e)
  inx <- match(tags, names(`_data`))
  matched <- !is.na(inx)
  if (any(matched)) {
    for (cname in tags[matched]) {
      `_data`[[cname]] <- e[[cname]]
    }
  }
  if (!all(matched)) {
    `_data` <- cbind(`_data`, as(e[!matched], 'DataFrame'))
  }
  `_data`
})

##' Checks that a directory exists and will create it if not.
##'
##' If the directory does not exist, and the caller does not want to create it
##' an error will be thrown
##'
##' @export
##' @author Steve Lianoglou \email{slianoglou@@gmail.com}
##'
##' @param path The path to the directory to check.
##' @param create A logical indicating whether or not the directory should be
##' created if it doesn't exist
##' @param verbose Let us know what's going on
##'
##' @return \code{TRUE} if everything is kosher, otherwise an error is thrown.
checkOrCreateDirectory <- function(path, create=FALSE, verbose=TRUE,
                                   raise.error=TRUE) {
  if (!dir.exists(path)) {
    if (!create) {
      if (raise.error) {
        stop("Directory", path, "does not exist", sep=" ")
      } else {
        return(FALSE)
      }
    } else {
      if (verbose) cat("Creating directory", path, "...\n")
      if (!dir.create(path)) {
        if (raise.error) {
          stop("Error! Check permissions? Parent directory exists?")
        } else {
          return(FALSE)
        }
      }
    }
  }

  path
}

checkVerbose <- function(...) {
  verbose <- list(...)$verbose
  if (is.null(verbose)) verbose <- options()$verbose
  verbose
}

dir.exists <- function(path) {
  path <- as.character(path)
  sapply(file.info(path)$isdir, isTRUE)
}

