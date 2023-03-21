restrictionSiteDistribution <- function(x, width=4) {
  leadingKmerDisribution(x, width)
}

leadingKmerDistribution <- function(x, k=4) {
  if (inherits(x, 'ShortRead')) {
    x <- sread(x)
  }
  stopifnot(inherits(x, 'XStringSet'))

  x <- narrow(x, start=1, width=k)
  counts <- table(as.character(x))
  ## counts <- counts[!is.na(counts)]
  ## counts <- counts[order(counts, decreasing=TRUE)]
  counts
}

##' @param counts named character vector. names() holds the kmer, and the value
##' is its count
plotDiscreteDistribution <- function(counts, top.n=5,
                                     title="Discrete Distribution",
                                     clean=TRUE, float.text=TRUE) {
  if (clean) {
    counts <- counts[!is.na(counts)]
  }
  
  o <- order(counts, decreasing=TRUE)
  is.top <- head(o, top.n)
  top.colors <- rainbow(top.n)
  
  bar.colors <- rep('grey', length(counts))
  bar.colors[is.top] <- top.colors
  bar.labels <- names(counts)[is.top]
  
  bp <- barplot(counts, main=title, xaxt='n', col=bar.colors)
  if (float.text) {
    text(bp[is.top],
         jitter(rep(.8 * max(counts), top.n), amount=.2*max(counts)),
         bar.labels, col=top.colors, cex=0.8)
  }
  legend('topleft', legend=bar.labels, text.col=top.colors)
  
  invisible(bp)
}

plotKmerDistribution <- function(...) {
  args <- list(...)
  if (is.null(args$title)) {
    args$title <- 'Kmer Distribution'
  }
  do.call(plotDiscreteDistribution, args)
}
