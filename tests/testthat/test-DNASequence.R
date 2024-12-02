library(testthat)

test_that("DNASequence constructor works correctly", {
  dna <- DNASequence("ATGCGC")
  expect_s3_class(dna, "DNASequence")
  expect_equal(attr(dna, "sequence"), "ATGCGC")

  #Test invalid DNA sequence
  expect_error(DNASequence("ATGX"), "Input is not a valid DNA sequence.")
})
