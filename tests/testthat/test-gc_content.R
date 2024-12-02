library(testthat)
library(stringr)

test_that("gc_content.S3 calculates GC content correctly", {

  # Test basic function
  dna <- DNASequence("ATGCGC")
  expect_equal(gc_content.S3(dna), 66.666667)

  #Test sequence with no GC base
  dna <- DNASequence("AAAA")
  expect_equal(gc_content.S3(dna), 0)

  #Test sequence with only GC base
  dna <- DNASequence("GGCC")
  expect_equal(gc_content.S3(dna), 100)

  #Test invalid DNA sequence
  expect_error(gc_content.S3("ATGC"), "Input must be of class 'DNASequence'.")
})
