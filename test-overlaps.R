context("More Overlaps")

test_that("Assigning overlaps by 'majority' works", {
  ir1 <- IRanges(c(1, 10, 20), width=9)
  ir2 <- IRanges(c(1, 5, 9, 11, 15, 18, 22), width=5)

  counts <- countUniqueOverlaps(ir1, ir2, 'quantify')
  expect_equal(sum(counts), length(ir2))
  expect_equal(counts, c(2, 3, 2))
})

test_that("assign unique overlaps by majority works", {
  ir1 <- IRanges(c(1, 3, 10, 15, 30), c(4, 12, 15, 20, 50))
  ir2 <- IRanges(c(1, 5, 10, 15), width=4)
  u <- assignUniqueOverlaps(ir1, ir2)
  expected <- c(1, 2, 3, 4, NA)
  expect_equal(u, expected)
})
