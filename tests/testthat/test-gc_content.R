library(testthat)
library(stringr)

test_that("gc_content.S3 calculates GC content correctly", {
  dna <- DNASequence("ATGCGC")
  expect_equal(gc_content.S3(dna), 66.67, tolerance = 0.01)
  dna <- DNASequence("AAAA")
  expect_equal(gc_content.S3(dna), 0)
  dna <- DNASequence("GGCC")
  expect_equal(gc_content.S3(dna), 100)
  expect_error(gc_content.S3("ATGC"), "Input must be of class 'DNASequence'.")
})
