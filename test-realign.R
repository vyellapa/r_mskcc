## The reads on the 3'utr of DICER and the multimapped locations of the same
## should act as a simple test case for remapping:
## 
## samtools view -H Mcf10A-1.preprocessed-sorted.bam > dicer.sam
## samtools view Mcf10A-1.preprocessed-sorted.bam chr5:67155938-67156111 >> dicer.sam
## samtools view Mcf10A-1.preprocessed-sorted.bam chr14:95552175-95557401 >> dicer.sam

context("Religning tags to unique lawn")

test_that("parsing BWA multimap (XA) tags works ...", {
  mm.string <- c("chr14,+19782316,35M,0;",
                 "chr5,+154194375,28M,0;chr5,+173363818,28M,2;chr2,+191537794,28M,2;",
                 NA,
                 "chr14,+19758953,24M,0;")

  ## Expected values
  sl <- c(chr14=NA_integer_, chr2=NA_integer_, chr5=NA_integer_)

  e1 <- GRanges('chr14', IRanges(19782316, width=35), strand='+', seqlengths=sl)
  values(e1) <- DataFrame(edit.distance=0L, cigar='35M', reference=1L)

  e2 <- GRanges(c('chr5', 'chr5', 'chr2'),
                IRanges(c(154194375, 173363818, 191537794), width=28),
                strand='+', seqlengths=sl)
  values(e2) <- DataFrame(edit.distance=c(0L, 2L, 2L), cigar=rep('28M', 3),
                          reference=c(2L, 2L, 2L))

  e3 <- GRanges('chr14', IRanges(19758953, width=24), strand='+', seqlengths=sl)
  values(e3) <- DataFrame(edit.distance=0L, cigar='24M', reference=4L)

  grl <- GRangesList(e1, e2, e3)
  result <- convertBwaMultimap(mm.string)
  result <- split(result, values(result)$reference)

  expect_identical(length(result), length(grl))
  for (idx in 1:length(grl)) {
    expect_identical(result[[idx]], grl[[idx]])
  }
})
