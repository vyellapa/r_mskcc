setGeneric("convolve1d", function(x, ...) standardGeneric("convolve1d"))

setMethod("convolve1d", c(x="numeric"),
function(x, kernel='normal', rescale=TRUE, bandwidth=20, ...) {
  if (is.character(kernel)) {
    kernel <- generateKernel(kernel, bandwidth=bandwidth, ...)
  }
  stopifnot(is.numeric(kernel))
  stopifnot(is.logical(rescale) && length(rescale) == 1L)

  # kernel <- rev(kernel)

  ret <- .Call("Rconvolve_1d", x, kernel, rescale, PACKAGE="biosignals")
  ret
})

setMethod("convolve1d", c(x="Rle"),
function(x, kernel='normal', rescale=TRUE, bandwidth=20, lower=0,
         eps=1e-6, ...) {
  views <- slice(x, lower=lower, includeLower=lower != 0)
  ret <- convolve1d(views, kernel, rescale, bandwidth, eps=eps, ...)
  subject(ret)
})

# June 2016: I believe convolve1d-RleViews was never called at a higher level
# For now it should be fine if we just remove it and see if Rconvolve_rle can be fixed in the future 
# Yuheng 
#
# setMethod("convolve1d", c(x="RleViews"),
# function(x, kernel='normal', rescale=TRUE, bandwidth=20, eps=1e-6, ...) {
  # if (is.character(kernel)) {
    # kernel <- generateKernel(kernel, bandwidth=bandwidth, ...)
  # }
  # stopifnot(is.numeric(kernel))
# 
  # x.rle <- subject(x)
  # islands <- ranges(x)
  # starts <- start(islands)
  # if (is.unsorted(starts)) {
    # stop("Starts in convolve1d,Rle must be sorted non-decreasing order")
  # }
# 
  # if (length(islands) == 0L) {
    # x.rle <- Rle(values=0, lengths=length(x.rle))
  # } else {
    # x.rle <- .Call("Rconvolve_rle", x.rle, kernel, starts, width(islands),
                   # rescale, eps, PACKAGE="biosignals")
  # }
# 
  # Views(x.rle, islands)
# })
