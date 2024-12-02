library(testthat)
library(ggplot2)

test_that("generate_dot_plot function works correctly", {

  # Test basic functionality
  seq1 <- "ATCGATCGATGC"
  seq2 <- "ATCCGATCGTAT"
  plot <- generate_dot_plot(seq1, seq2)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 1)
  expect_s3_class(plot$layers[[1]]$geom, "GeomPoint")
})
