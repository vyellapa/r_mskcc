context("Smoothing Reads")

## An object represeting "ragged" reads rep
ragged <- GRanges(seqnames='chr1', strand='+',
                  ranges=IRanges(start=c(sample(18:22, 10, replace=TRUE),
                                         sample(100:108, 10, replace=TRUE)),
                                 width=sample(18:20, 20, replace=TRUE)))

test_that("repacking tags evens-out ends", {
  
  group.1 <- ragged[1:10]
  group.2 <- ragged[11:20]
  expected.starts <- c(min(start(group.1)), min(start(group.2)))
  expected.ends <- c(max(end(group.1)), max(end(group.2)))
  
  clean <- smoother.repackTags(ragged)
  
  expect_that(unique(start(clean)), equals(expected.starts),
              info="Starts aren't clean")
  expect_that(unique(end(clean)), equals(expected.ends),
              info="Ends aren't clean")
})
