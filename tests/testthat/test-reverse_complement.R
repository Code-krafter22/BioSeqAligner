library(testthat)
library(stringi)

test_that("reverse_complement returns correct results", {

  # Test basic function
  expect_equal(reverse_complement("ATCGCTCA"), "TGAGCGAT")

  # Test single base
  expect_equal(reverse_complement("AAAA"), "TTTT")

  #Test invalid DNA sequence
  expect_error(reverse_complement("ATCX"), "Input is not a valid DNA sequence.")
})
