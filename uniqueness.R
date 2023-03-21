##' Create genomic chromosome FASTA files by masking out non-transcribed space
##'
##' Use these chromosome files to make a uniqueness map for aligning
##' tag-sequenced experiments.
createTranscribedGenome <- function(gcache, flank.up=1000L,
                                    flank.down=flank.up, path='.') {
  xcripts <- transcripts(gcache)
  if (flank.up > 0L) {
    xcripts <- resize(xcripts, width=width(xcripts) + flank.up, fix="end")
  }
  if (flank.down > 0L) {
    xcripts <- resize(xcripts, width=width(xcripts) + flank.down, fix="start")
  }

  xcripts <- reduce(xcripts)
  bsg <- getBsGenome(gcache)
  for (.xcripts in split(xcripts, seqnames(xcripts))) {
    chr <- as.character(seqnames(.xcripts)[1])
    cat(chr, "...\n")
    bsg.chr <- unmasked(bsg[[chr]])
    genomic <- gaps(.xcripts) ## returns gaps for all other chromosomes too!
    genomic <- genomic[seqnames(genomic) == chr]
    genomic <- ranges(genomic[strand(genomic) != '*'])
    for (i in head(1:length(genomic))) {
      cat(i, "...\n")
      repl <- DNAString(paste(rep("N", width(genomic[i])), collapse=""))
      subseq(bsg.chr, start(genomic[i]), end(genomic[i])) <- repl
    }
    ## Write the fasta file
    dset <- DNAStringSet(bsg.chr)
    names(dset) <- chr
    write.XStringSet(dset, paste(chr, 'fa', sep='.'))
  }
}
