library(testthat)
library(stringi)

test_that("transcribe_dna returns correct RNA sequence", {

  # Test basic function
  expect_equal(transcribe_dna("ATGC"), "AUGC")
  expect_equal(transcribe_dna("TTAA"), "UUAA")

  #Test invalid DNA sequence
  expect_error(transcribe_dna("ATGX"), "Input is not a valid DNA sequence.")
})
