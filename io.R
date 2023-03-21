##' Expand a frequency count file to one that lists each repitition of a tag on
##' its own line
frequencyCountToTags <- function(freq.file, outfile, prepend=NULL,
                                 append=NULL,...) {
  tags <- read.table(freq.file, ..., stringsAsFactors=FALSE)
  if (!is.null(prepend)) {
    tags[,1] <- paste(prepend, tags[,1], sep="")
  }
  if (!is.null(append)) {
    tags[,1] <- paste(tags[,1], append, sep="")
  }
  
  tags.ex <- rep(tags[,1], times=tags[,2])
  writeLines(tags.ex, outfile)
}

tagsToFrequencyCount <- function(tag.file, freq.file) {

}
