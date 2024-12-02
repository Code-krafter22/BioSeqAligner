library(testthat)
library(stringi)

test_that("reverse_complement returns correct results", {
  expect_equal(reverse_complement("ATCGCTCA"), "TGAGCGAT")
  expect_equal(reverse_complement("AAAA"), "TTTT")
  expect_error(reverse_complement("ATCX"), "Input is not a valid DNA sequence.")
  expect_equal(reverse_complement(""), "")
})
