## A specialized version in peaks.R was built for fenceposting events in a
## smooth signal using Rcpp as in R/peaks.R. It makes some assumptions but
## ends up being much faster.

setOldClass('turnpoints')

turns <- function(x, eps=1e-8) {

  n <- length(x)
  diffs <- abs(c(x[1] - 1, x[-n]) - x) > eps
  uniques <- x[diffs]
  n.unique <- length(uniques)

  poss <- (1:n)[diffs] ## position of unique points
  exaequos <- c(poss[2:n.unique], n + 1) - poss - 1

  if (n.unique < 3) {
    ## Less than 3 values can't give us anything goold
    ## TODO: return dummy value
    return(NULL)
  }

  ex <- embed(uniques, 3)
  e1 <- ex[,1]
  e2 <- ex[,2]
  e3 <- ex[,3]
  peaks <- c(FALSE, pmax(e1, e2, e3) == e2, FALSE)
  pits <- c(FALSE, pmin(e1, e2, e3) == e2, FALSE)

  ret <- list()

  if (sum(peaks) > 0L) {
    ret$peaks.which <- (poss + exaequos)[peaks]
  } else {
    ret$peaks.which <- integer()
  }


  if (sum(pits) > 0L) {
    ret$pits.which <- (poss + exaequos)[pits]
  } else {
    ret$pits.which <- integer()
  }

  ret
}


if (FALSE) {
  ## This is rewritten for our particular usecase in peaks.R
  ## TODO: Optimze along with `findPeakBoundsInSignal` by assuming that `x`
  ##       is always a 'decent' coverage vector. This will add minima to
  ##       the first & last non-zero elements in x, and will put a maxima
  ##       in between minima fenceposts if non is found (ie. from a flat coverage
  ##       vector)
  turnfun.rcpp <- cxxfunction(signature(x="numeric", eps="numeric"), '
    Rcpp::NumericVector xx(x);
    Rcpp::NumericVector epsilon(eps);
    std::vector<int> maxima;
    std::vector<int> minima;
    double pslope = 0;
    double pvalue = 0;
    double slope = 0;
    double value = 0;
    double precision = epsilon[0];
    int run_start = 0;
    int run_is_peak = 0;
    int nxx = xx.size();
    int i;
    int no_maxima;

    for (i = 0; i < nxx; i++) {
      value = xx[i];
      slope = value - pvalue;

      if ((-precision < slope) && (slope < precision)) {
        // TODO: Slope is essentially 0 -- handle a plateu case here
        pslope = 0;
        pvalue = value;
      } else {
        if (slope < 0 && pslope >= 0) {
          // previous point was a maxima
          maxima.push_back(i);
        } else if (slope > 0 && pslope <= 0) {
          // previous point was a minima or this if first non-zero element
          minima.push_back(i);
        }

        pvalue = value;
        pslope = slope;
      }
    }

    // Add a ghost maxima to deal with flat coverage vector?
    no_maxima = (maxima.size() == 0) ? 1 : 0;
    if (no_maxima) {
      maxima.push_back((int) nxx / 2);
    }

    if (minima.size() > 0) {
      if (minima[minima.size() - 1] < maxima[maxima.size() - 1]) {
        // add a ghost end
        minima.push_back(nxx);
      }
      if (minima[0] == 0) {
        minima[0] = 1; // R indexing
      }
    }

    return Rcpp::List::create(Rcpp::Named("maxima") = maxima,
                              Rcpp::Named("minima") = minima,
                              Rcpp::Named("ghost.maxima") = no_maxima);
  ', plugin="Rcpp")
}

turnfun <- function(x, eps=1e-6) {
  turnfun.rcpp(x, eps)
}

##' Locates the turnpoints in a "smooth" series.
##'
##' This function was taken/modified from the turnpoints function in the
##' pastecs package.
##'
##' @param x The object to extract turnpoints from.
##' @param eps How close are two numbers before they are considered "equal"?
##' @return A \code{turnpoints} object
turnpoints <- function(x, eps=1e-8) {
  if (!is.null(dim(x))) {
    stop("Can only work on vectors, not matrices")
  }
  data <- deparse(substitute(x))
  x <- as.numeric(x)
  n <- length(x)
  ## diffs <- c(x[1] - 1, x[-n]) != x
  diffs <- abs(c(x[1] - 1, x[-n]) - x) > eps
  uniques <- x[diffs]
  n2 <- length(uniques)
  poss <- (1:n)[diffs] ## position of unique points
  exaequos <- c(poss[2:n2], n+1) - poss - 1
  if (n2 < 3) {
    warning("Less than 3 unique values, no calculation!")
    nturns <- NA
    firstispeak <- FALSE
    peaks <- rep(FALSE, n2)
    pits <- rep(FALSE, n2)
    tppos <- NA
    proba <- NA
    info <- NA
  } else {
    ## The original function had two versions of the code -- I think the note
    ## about it not working all the time using the `embed` + `max.col` route
    ## might have been because of how ties are broken with max.col -- I'm using
    ## a modification of this version in a way I think should work.
    ##
    ## Original (2) versions of the code
    ## The following code is faster in R, but do not work all the time!
    ## ex <- embed(uniques, 3) # Works only in R!
    ## peaks <- c(FALSE, max.col(ex) == 2, FALSE)
    ## pits <- c(FALSE, max.col(-ex) == 2, FALSE)
    ##
    ## Not sure why the above won't work all the time, but this was the
    ## the substituted code -- why they use apply I dunno
    m <- n2 - 2
    ex <- matrix(uniques[1:m + rep(3:1, rep(m, 3)) - 1], m)
    peaks <- c(FALSE, apply(ex, 1, max, na.rm=TRUE) == ex[, 2], FALSE)
    pits <- c(FALSE, apply(ex, 1, min, na.rm=TRUE) == ex[, 2], FALSE)
    tpts <- peaks | pits
    if (sum(tpts) == 0) { # No turning point
      nturns <- 0
      firstispeak <- FALSE
      peaks <- rep(FALSE, n2)
      pits <- rep(FALSE, n2)
      tppos <- NA
      proba <- NA
      info <- NA
    } else {
      ## Consider the last element of duplicates
      tppos <- (poss + exaequos)[tpts]
      tptspos <- (1:n2)[tpts]
      firstispeak <- tptspos[1] == (1:n2)[peaks][1]
      nturns <- length(tptspos)
      if (nturns < 2) {
        inter <- n2 + 1
        posinter1 <- tptspos[1]
      } else {
        inter <- c(tptspos[2:nturns], n2) - c(1, tptspos[1:(nturns-1)]) + 1
        posinter1 <- tptspos - c(1, tptspos[1:(nturns-1)])
      }
      posinter2 <- inter - posinter1
      posinter <- pmax(posinter1, posinter2)
      proba <- 2 / (inter * gamma(posinter) * gamma(inter - posinter + 1))
      info <- -log(proba, base=2)
    }
  }

  res <- list(data=data, n=n, points=uniques, pos=(poss + exaequos),
              exaequos=exaequos, nturns=nturns, firstispeak=firstispeak,
              peaks=peaks, pits=pits, tppos=tppos, proba=proba, info=info)
  class(res) <- "turnpoints"
  res
}


##' Extracts part of a dataset according to an analysis
##'
##' @param x The object on which extraction is performed
##' @param n The number of elements from the original dataset to extract
##' @param ... Additional parameters affecting the extraction
setGeneric("extract", function(x, n, ...) standardGeneric("extract"))

##' @nord
setMethod("extract", c(x="turnpoints"),
function(x, n, no.tp=0L, peak=1L, pit=-1L, ...) {
  if (missing(n)) {
    n <- -1L
  }
  res <- rep(no.tp, length.out=x$n)
  res[x$pos[x$peaks]] <- peak
  res[x$pos[x$pits]] <- pit

  # Keep only the first n points
  if (n < length(res) && n > 0) {
    res <- res[1:n]
  }

  res
})

##' Locates maxima in a continues signal
##'
##' @param x A numeric, Rle, or turnpoints object containing the signal
##' @param as If \code{index}, the index with the maxima into \code{x},
##' is returned. If \code{value}, the value of the signal at each maxmima.
##'
##' @return The indices (or values) of the maxima in \code{x}.
setGeneric("maxima", function(x, as=c('index', 'value'), ...) {
  standardGeneric("maxima")
})

##' @nord
setMethod("maxima", c(x="turnpoints"),
function(x, as, ...) {
  as <- match.arg(as)
  tp <- extract(x, ...)
  switch(as, index=which(tp == 1), value=x$points[x$peaks])
})

##' @nord
setMethod("maxima", c(x="numeric"),
function(x, as, eps=1e-8, ...) {
  maxima(turnpoints(x, eps=eps), as, ...)
})

##' @nord
setMethod("maxima", c(x="Rle"),
function(x, as, eps=1e-8, ...) {
  maxima(turnpoints(as.numeric(x), eps=eps), as, ...)
})


##' Locates minima in a continues signal
##'
##' @param x A numeric, Rle, or turnpoints object containing the signal
##' @param as If \code{index}, the index with the maxima into \code{x},
##' is returned. If \code{value}, the value of the signal at each maxmima.
##'
##' @return The indices (or values) of the minima in \code{x}.
setGeneric("minima", function(x, as=c('index', 'value'), ...) {
  standardGeneric("minima")
})

##' @nord
setMethod("minima", c(x="turnpoints"),
function(x, as=c('index', 'value'), ...) {
  as <- match.arg(as)
  tp <- extract(x, ...)
  switch(as, index=which(tp == -1), value=x$points[x$pits])
})

##' @nord
setMethod("minima", c(x="numeric"),
function(x, as, eps=1e-8, ...) {
  minima(turnpoints(x, eps=eps), as, ...)
})

##' @nord
setMethod("minima", c(x="Rle"),
function(x, as, eps=1e-8, ...) {
  minima(turnpoints(as.numeric(x), eps=eps), as, ...)
})

