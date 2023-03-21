##' The main class representing a TailSeq experiment. It contains the aligned
##' reads along with other metadata for the experiment.
TailSeq <- function(directory, ...) {
  ## Check for bam files and ignore any from restriction-info (lumiDGE)
  ## TODO Use the .bamAndIndexPath function to cut boilerplate
  ss <- SeqStore(directory, 'TailSeq', read.length=c(20:22), is.stranded=TRUE,
                 is.paired=FALSE, ...)
  ss
}
