## These functions are deprecated and use "the old" way of annotating each read
## by iterating over a gene-by-gene basis. Look in summaries.R for newer methods
## that use annotatedChromosomes.
##
## Functions to summarise tags into an excel-like report format
##
## These were written to get summary stats for Christine's grant during the
## week of September 8, 2010. They rely on a cached set of idealized gene models
## that were also built at the same time.
##
## There is mysterious problem with islandCoverage smoothing when looping
## so I'm saving after each iteration.

##' @param expt The SeqStore object to summarize tags over
##' @param gene.models A list of list of gene models
generateSummaries <- function(expt, gene.models, unique.only=TRUE, min.cover=2,
                              chrs=paste('chr', c(1:22,'X','Y','M'), sep=""),
                              smooth.by=c('islandCoverage', 'repackTags'),
                              sample.name=NULL, path='.') {
  sums <- vector('list', length(chrs))
  names(sums) <- chrs
  
  for (chr in chrs) {
    cat(chr, "\n")
    fname <- 'summaries'
    if (!is.null(sample.name)) {
      fname <- paste(fname, sample.name, sep=".")
    }
    fname <- sprintf('%s.%s.unique-%s.cover%d.rda', fname, chr, unique.only,
                     min.cover)
    
    reads <- getReadsOnChromosome(expt, chr, smooth.by=NULL,
                                  unique.only=unique.only)
    
    sm <- tryCatch({
      smoothReads(reads, smooth.by=smooth.by, min.cover=min.cover)
    }, error=function(e) NULL)
    
    if (is.null(sm)) {
      cat("!!! ERROR -- saving b0rked reads!!!\n")
      save(reads, file=file.path(path, paste('b0rked', fname, sep=".")))
    } else {
      cat("... summarizing ...\n")
      models <- gene.models[[chr]]
      cs <- summarizeTags(sm, models)
      if (sum(cs$nreads == -1) > nrow(cs) / 3) {
        cat("////////", chr, "summary is suspicious. ///////\n")
        fname <- paste('suspicious', fname, sep=".")
        save(cs, file=file.path(path, fname))
      } else {
        cat("... writing file to:", file.path(path, fname), "\n")
        save(cs, file=file.path(path, fname))
      }
      
      sums[[chr]] <- cs
    }
  }
  
  sums <- do.call(rbind, sums)
  fname <- 'summaries'
  if (!is.null(sample.name)) {
    fname <- paste(fname, sample.name, sep=".")
  }
  fname <- sprintf('%s.unique-%s.cover%d.rda', fname, unique.only, min.cover)
  save(sums, file=file.path(path, fname))
  sums
}

##' Assumes that the reads and gene.models list are on the same chromosome
##' Using cleaned reads, return a table with:
##'     symbol   chr   strand   nreads   npiles   utr3.piles   utr3.count
##'     intron.piles   intron.count   cds.piles   cds.count   utr5.piles
##'     utr5.count
summarizeTags <- function(reads, gene.models, stranded=TRUE) {
  fwd.reads <- ranges(reads[strand(reads) == '+'])
  rev.reads <- ranges(reads[strand(reads) == '-'])
  all.reads <- ranges(reads)

  .empty <- data.frame(chr='', strand='', nreads=0L, npiles=0L, utr3.piles=0L,
                       utr3.count=0L, intron.piles=0L, intron.count=0L,
                       cds.piles=0L, cds.count=0L, utr5.piles=0L, utr5.count=0L)
  
  summaries <- ldply(gene.models, function(gm) {
    g.strand <- as.character(strand(gm))
    if (length(g.strand) == 0) {
      return(.empty)
    }
    g.strand <- g.strand[1]
    if (stranded) {
      the.reads <- if (g.strand == '+') fwd.reads else rev.reads
    } else {
      the.reads <- all.reads
    }
    
    s <- tryCatch(summarizeTagsOverIdealizedGene(gm, the.reads),
                  error=function(e) NULL)
    if (is.null(s)) {
      cat(geterrmessage())
      s <- .empty
      s$nreads <- -1L
    }
    s
  }, .progress='text')

  summaries
}

.summarize <- function(.range, .reads) {
  if (length(.range) == 0L) {
    return(c(npiles=0L, nreads=0L))
  }
  
  .reads <- subsetByOverlaps(.reads, .range)
  if (length(.reads) > 0L) {
    c(npiles=length(reduce(.reads)), nreads=length(.reads))
  } else {
    c(npiles=0L, nreads=0L)
  }
}

## Expect reads to be iranges object (all strand info out the window)
summarizeTagsOverIdealizedGene <- function(igene, reads) {
  if (is.null(igene) || length(igene) == 0) {
    df <- data.frame(chr='', strand='', nreads=0L, npiles=0L, utr3.piles=0L,
                     utr3.count=0L, intron.piles=0L, intron.count=0L,
                     cds.piles=0L, cds.count=0L, utr5.piles=0L, utr5.count=0L)
    return(df)
  }

  g.strand <- as.character(strand(igene[1]))
  g.chr <- as.character(seqnames(igene[1]))

  g.ranges <- ranges(igene)
  g.anno <- values(igene)$exon.anno
  these.reads <- subsetByOverlaps(reads, range(g.ranges))
  
  nreads <- length(these.reads)
  piles <- reduce(these.reads)
  npiles <- length(piles)

  ## Narrow width of reads to be 1, so that they only are counted in
  ## one category
  narrowed <- resize(these.reads, width=1,
                     fix=if (g.strand == '+') 'start' else 'end')
  ## utr3
  is.utr3 <- g.anno == 'utr3' | g.anno == 'utr3*'
  utr3.summary <- .summarize(g.ranges[is.utr3], narrowed)
  
  ## intron
  intron.summary <- .summarize(gaps(g.ranges), narrowed)
  
  ## cds
  cds.summary <- .summarize(g.ranges[g.anno == 'cds'], narrowed)
  
  ## utr5
  is.utr5 <- g.anno == 'utr5' | g.anno == 'utr5*'
  utr5.summary <- .summarize(g.ranges[is.utr5], narrowed)

  data.frame(chr=g.chr, strand=g.strand, nreads=nreads, npiles=npiles,
             utr3.piles=utr3.summary[1], utr3.count=utr3.summary[2],
             intron.piles=intron.summary[1], intron.count=intron.summary[2],
             cds.piles=cds.summary[1], cds.count=cds.summary[2],
             utr5.piles=utr5.summary[1], utr5.count=utr5.summary[2])
}

## Load the summary rdas generated from the experiment from generateSummaries
loadSummaries <- function(min.cover=2, unique.only=TRUE, sample.name=NULL,
                          path='.') {
  fname <- 'summaries'
  if (!is.null(sample.name)) {
    fname <- paste(fname, sample.name, sep=".")
  }
  regex <- sprintf('^%s.*unique-%s.cover%d', fname, unique.only, min.cover)
  files <- list.files(path, pattern=regex, full.names=TRUE)
  sums <- ldply(files, function(fname) {
    name <- load(fname)
    get(name)
  })
  names(sums)[1] <- 'symbol'
  sums
}
