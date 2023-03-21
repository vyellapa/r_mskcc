#!/usr/bin/env Rscript

details <- paste(
  "Generates an annotated tag distribution, indicating the proportion of reads",
  "that land in different (annotated) regions of the genome\n",
  "Packages required:",
  "  + optparse",
  "  + Rsamtools",
  "  + GenomicRanges",
  sep="\n")

###############################################################################
## Utilitfy functions
## ----------------------------------------------------------------------------

##' Bomb out of the script with an error message.
##'
##' @param ... A character vector of "message pieces" to string together
##' @param sep The separator to use to concat the pieces together
err <- function(..., sep=" ") {
  msg <- paste(unlist(list(...)), sep=sep)
  cat("ERROR:", msg, "\n")
  quit(save='no', status=1)
}

###############################################################################
## Load up `optparse` to configure and parse command line arguments.
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  if (!require(optparse)) {
    err('The optparse library needed, please install via:\n',
        'install.packages("optparse", repos="http://R-Forge.R-project.org")')
  }
})

usage <- paste("%prog [OPTIONS] BAMFILE GENOMEANNO\n\n", details)
option.list <- 
  list(make_option(c('-n', '--no.strand'), default=FALSE, action="store_true",
                   help="Set to ignore the strand of alignments."),
       make_option(c('-v', '--verbose'), default=FALSE, action="store_true",
                   help="Keeps you informed of what's happening."))

parser <- OptionParser(usage=usage, option_list=option.list)
parsed <- parse_args(parser, positional_arguments=TRUE)
opts <- parsed$options
args <- parsed$args

verbose <- opts$verbose
if (opts$no.strand) {
  err("Ignoring strand is not implemented yet -- left as exercise for reader.")
}

###############################################################################
## Check option/arguments
## ----------------------------------------------------------------------------
if (length(args) != 2) {
  err("2 arguments required")
}

bam.file <- args[1]
if (!file.exists(bam.file)) {
  err("BAM file not found:", bam.file)
}

if (!file.exists(args[2])) {
  err("Genome annotation not found:", args[2])
}

if (verbose) cat("Loading genome annotation object ...\n")
tmpenv <- new.env()
g.name <- load(args[2], tmpenv)
genome.anno <- get(g.name[1], tmpenv, inherits=FALSE)
if (!inherits(genome.anno, 'GRanges')) {
  err("Genome Annotation object looks wrong:", class(genome.anno))
}

###############################################################################
## Everything looks kosher, let's do the deed.
## ----------------------------------------------------------------------------
if (verbose) {
  cat("Loading required packages ...\n")
}
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
})

if (verbose) {
  cat("Loading aligned reads ...\n")
}
what <- c('rname', 'strand', 'pos', 'qwidth')
flags <- scanBamFlag(isUnmappedQuery=FALSE)
aln <- scanBam(bam.file, param=ScanBamParam(what=what, flag=flags))[[1]]
gr <- GRanges(aln$rname, IRanges(aln$pos, width=aln$qwidth), aln$strand)

if (verbose) {
  cat("Tallying distribution of reads ...\n")
}

## Note that you might be introducing some sort of bias by only selecting the
## "first" match, but findOverlaps,GRanges doesn't support select="arbitrary"
## so maybe you can live that for now.
idx <- findOverlaps(gr, genome.anno, select='first')
nomatch <- is.na(idx)
idx <- idx[!nomatch]

distro <- table(values(genome.anno)$exon.anno[idx])
distro <- data.frame(anno=names(distro), count=as.integer(distro))

if (any(nomatch)) {
  distro <- rbind(distro, data.frame(anno="unknown", count=sum(nomatch)))
}

show(distro)
