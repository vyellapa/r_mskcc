##' Creates a quality metrics report on a TagSeq dataset
##'
##' This delegates to the methods the sub-QC methods, which are responsible
##' for generating their own portion of the report.
##'
##' Each sub-QC method should have the same signature:
##'   x          : TagSeq experimetn
##'   gcache     : The GenomicCache object to use as a reference
##'   atags=NULL : The annotated tags (fetched if null)
##'   fh         : The file handle to write to
##'   path       : The directory the report is being written to
##'   ...        : Other vars to pass through to internal function calls
##'
##' @param x The \code{TagSeq} object
##' @param metrics A list of function names to use for report generation.
##' If \code{NULL}, a default set of qc stats will be presented.
##' @param path The directory to dump the report into
setGeneric("qcReport",
function(x, gcache, atags=NULL, metrics=NULL, path='.', ...) {
  standardGeneric("qcReport")
})

setMethod("qcReport", c(x="TagSeq"),
function(x, gcache, atags=NULL, metrics=NULL, path='.', ...) {
  stopifnot(inherits(x, 'TagSeq'))
  stopifnot(inherits(gcache, 'GenomicCache'))
  if (is.null(metrics)) {
    metrics <- c('annotatedTagDistribution', 'topGenes')
  }
  if (is.null(atags)) {
    atags <- getAnnotatedReads(x, gcache, ...)
  }

  fh <- file(file.path(path, 'index.html'), 'w')
  on.exit(close(fh))

  ## HTML Header
  cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"\n', file=fh)
  cat('  "http://www.w3.org/TR/html4/strict.dtd">\n', file=fh)
  cat('<html>\n\n', file=fh)
  cat('<head>\n', file=fh)
  cat('<title>TagSeq QC Report</title>\n', file=fh)
  cat('<link href="tagseq-qc.css" rel="stylesheet" type="text/css">\n', file=fh)
  cat('</head>\n\n', file=fh)
  cat('<body>\n', file=fh)

  ## Print TOC
  cat('<div id="toc">\n', file=fh)
  cat('<ul>\n', file=fh)
  for (metric in metrics) {
    cat(sprintf('<li><a href="#%s">%s</li>\n', metric, metric), file=fh)
  }
  cat('</ul>\n', file=fh)
  cat('</div>\n', file=fh)

  ## Delegate
  err <- function(e) {
    cat("Internal QC block failed, moving on ...\n")
  }
  for (metric in metrics) {
    cat("...", metric, "...\n")
    cat(sprintf('<a name="%s">&nbsp;</a>', metric), file=fh)
    cat(sprintf('<div id="%s" class="qc-section">\n', metric), file=fh)

    tryCatch({
      getQcMethod(metric)(x, gcache, atags, fh, path=path, ...)
    }, error=function(e) cat("QC Block failed, moving on ...\n"))

    cat(sprintf('</div> <!-- %s -->\n\n', metric), file=fh)
  }

  cat('</body>\n', file=fh)
  cat('</html>\n', file=fh)
})

getQcMethod <- function(what) {
  if (substring(what, 1, 3) != 'qc.') {
    what <- paste('qc', what, sep='.')
  }
  fn <- getFunction(what)
  if (!is.function(fn)) {
    stop("Unknown function `", what, "`")
  }
  fn
}

qc.annotatedTagDistribution <- function(x, gcache, atags=NULL, fh="", path='.',
                                        min.cover=iTPM(3, x), ...) {
  if (is.null(atags)) atags <- getAnnotatedReads(x, gcache, ...)

  ##################################################################### Raw Tags
  png(file.path(path, "expression.tag.distribution.png"), 800, 800)
  etags <- annotatedTagDistribution(atags, expression.weighted=TRUE, do.plot=TRUE)
  dev.off()

  png(file.path(path, "genomic.tag.distribution.png"), 800, 800)
  gtags <- annotatedTagDistribution(atags, expression.weighted=FALSE,
                                    do.plot=TRUE)
  dev.off()

  cat('<h2>Annotated Tag Distributions</h2>\n', file=fh)
  cat('<table border="0" cellpadding="2" cellspacing="0">\n', file=fh)
  cat('<tr><th>Tag Distribution</th><th>Genomic Positions</th></tr>\n', file=fh)
  cat('<tr>\n', file=fh)
  cat(sprintf('<td><a href="%s"><img src="%s" width="400" /></a></td>\n',
              "expression.tag.distribution.png",
              "expression.tag.distribution.png"), file=fh)
  cat(sprintf('<td><a href="%s"><img src="%s" width="400" /></a></td>\n',
              'genomic.tag.distribution.png', 'genomic.tag.distribution.png'),
      file=fh)
  cat('</tr>\n</table>\n', file=fh)

  ################################################################# Cleaned tags
  png(file.path(path, "clean.expression.tag.distribution.png"), 800, 800)
  etags <- annotatedTagDistribution(atags, expression.weighted=TRUE,
                                    do.plot=TRUE, min.count=min.cover)
  dev.off()

  png(file.path(path, "clean.genomic.tag.distribution.png"), 800, 800)
  gtags <- annotatedTagDistribution(atags, expression.weighted=FALSE,
                                    do.plot=TRUE, min.count=min.cover)
  dev.off()

  cat(sprintf('<h2>Cleaned tag distribution (> %d support)</h2>\n',
              floor(min.cover)), file=fh)
  cat('<table border="0" cellpadding="2" cellspacing="0">\n', file=fh)
  cat('<tr><th>Tag Distribution</th><th>Genomic Positions</th></tr>\n', file=fh)
  cat('<tr>\n', file=fh)
  cat(sprintf('<td><a href="%s"><img src="%s" width="400" /></a></td>\n',
              "clean.expression.tag.distribution.png",
              "clean.expression.tag.distribution.png"), file=fh)
  cat(sprintf('<td><a href="%s"><img src="%s" width="400" /></a></td>\n',
              'clean.genomic.tag.distribution.png',
              'clean.genomic.tag.distribution.png'),
      file=fh)
  cat('</tr>\n</table>\n', file=fh)
  cat("</table>", file=fh)
}

qc.topGenes <- function(x, gcache, atags=NULL, fh="", path='.', n.genes=10,
                        tag.type=c('cds', 'utr3', 'utr3*'), ...) {
  cat(sprintf('<h2>Reads over the top %d most highly expressed genes</h2>\n',
              n.genes), file=fh)
  if (is.null(atags)) atags <- getAnnotatedReads(x, ...)
  tag.table <- tabularize(atags)
  ## tag.table <- subset(tag.table, exon.anno %in% tag.type)
  expr <- tag.table[, list(count=sum(count)), by='symbol']
  expr <- expr[!is.na(expr$symbol)]
  take <- order(expr$count, decreasing=TRUE)[1:n.genes]

  .processGene <- function(name) {
    gene <- GFGene(name, gcache)
    cat("   ...", symbol(gene), "\n")
    html <- paste(
      '<table border="0" cellpading="2" cellspacing="0">',
      sprintf('<tr class="gene-name"><td colspan="2">%s</td></tr>', symbol(gene)),
      '<tr><th>Both strands</th><th>Gene strand only</th></tr>\n', sep="\n")

    same <- paste(symbol(gene), 'same-strand', 'png', sep='.')
    all <- paste(symbol(gene), 'png', sep='.')

    png(file.path(path, same), 800, 1000)
    plotReads(gene, list(tags=x), reads.strand=strand(gene), ...)
    dev.off()

    png(file.path(path, all), 800, 1000)
    plotReads(gene, list(tags=x), ...)
    dev.off()

    html <- paste(html,
      '<tr class="plot-reads">\n',
      sprintf('<td><a href="%s"><img src="%s" width="400" /></a></td>\n', all, all),
      sprintf('<td><a href="%s"><img src="%s" width="400" /></a></td>\n', same, same),
      '</tr>', sep="\n")

    tryCatch({
      ttg <- basename(ttview(x, gene, path=path))
    }, error=function(e) NULL)

    if (!is.null(ttg)) {
      html <- paste(html,
        '<tr class="text-view"><td colspan="2">',
        sprintf('<a href="%s">Text View</a>', ttg),
        '</td></tr>', sep="\n")
    } else {
      cat("ttview error for", symbol(gene), "\n")
    }

    html <- paste(html, '</table>', sep="\n")
    cat(html, file=fh)
    gene
  }
  genes <- lapply(as.character(expr[take]$symbol), function(name) {
    gene <- tryCatch(.processGene(name), error=function(e) NULL)
    if (is.null(gene)) {
      cat("failed processing gene", name, "... moving on\n")
    }
    gene
  })


  invisible(genes)
}


## Counts per read
countsPerRead <- function(annotated.reads, as.log=TRUE) {
  counts <- values(annotated.reads)$count
  if (as.log) {
    counts <- log10(counts)
  }
  plot(ecdf(counts), verticals=TRUE, do.points=FALSE)
}

readComplexity <- function(annotated.reads, bsg, k=3, as.fraction=TRUE,
                           do.square=FALSE, with.sequences=FALSE) {
  chrs <- as.character(unique(seqnames(annotated.reads)))
  cplx <- ldply(chrs, function(chr) {
    cat("===", chr, "===\n")
    bsg.chr <- unmasked(bsg[[chr]])
    anno <- annotated.reads[seqnames(annotated.reads) == chr]
    .ranges <- ranges(anno)
    max.scores <- width(.ranges) - k
    max.squared <- max.scores^2

    seqs <- Views(bsg.chr, .ranges)
    tri <- trinucleotideFrequency(seqs) - 1L
    nkmers <- rowSums(tri >= 0)
    tri[tri < 0] <- 0L
    scores.squared <- rowSums(tri * tri)
    scores <- rowSums(tri)

    df <- data.frame(chr=chr, exon.anno=values(anno)$exon.anno,
                     count=values(anno)$count,
                     score=scores, score.sq=scores.squared,
                     norm.score=scores / max.scores,
                     norm.score.sq=scores.squared / max.squared,
                     nkmers=nkmers)
    if (with.sequences) {
      df$sequences <- as.character(seqs)
    }

    df
  })
  cplx
}


