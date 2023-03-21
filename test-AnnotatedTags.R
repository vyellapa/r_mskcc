context("Compressing Tags")
################################################################################
## Test Data
################################################################################
tag.chrs <- c("chr21", "chr21", "chr22", "chr21", 'chr21')
tag.times <- c(5,       5,      8,       1,       3)
tag.start <- c(1000,    1000,   1000,    100,     105)
tag.width <- c(21,      21,     21,      21,      5)
tag.strand <- c("+",    "-",    "+",     "-",     "-")
ir <- IRanges(start=rep(tag.start, tag.times), width=rep(tag.width, tag.times))
test.tags <- GRanges(
  seqnames=rep(tag.chrs, tag.times), ranges=ir,
  strand=rep(tag.strand, tag.times))
## compressed <- compressReads(test.tags)

################################################################################
## GO
################################################################################
.check <- function(compressed, unrolled) {
  compressed <- as.data.frame(compressed)
  unrolled <- as.data.frame(unrolled)

  for (idx in 1:nrow(compressed)) {
    g <- got[idx,]
    ref <- subset(original, seqnames == g$seqnames & start == g$start & end == g$end
                  & strand == g$strand)
    expect_equal(g$count, nrow(ref))
  }
}
test_that("Compressed reads ranges returns correct counts ...", {
  .check(compressReads(test.tags), test.tags)
})

test_that("Compressing reads with counts adds to existing counts ...", {
  g1 <- compressReads(test.tags[1:13])
  g2 <- compressReads(test.tags[-(1:13)])
  final <- compressReads(c(g1, g2))
  .check(final, test.tags)
  expect_identical(final, compressReads(test.tags))
})
