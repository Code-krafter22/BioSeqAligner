library(testthat)

test_that("find_motif function works correctly", {
  #Test basic functions
  expect_equal(find_motif("ATCGATCGATGC", "ATC"), c(1, 5))
  expect_equal(find_motif("GATCAGATCGTAT", "GAT"), c(1, 6))

  #Test no motif match
  expect_equal(find_motif("AAAAAA", "TT"), NULL)

  #Test exact match
  expect_equal(find_motif("ATCG", "ATCG"), 1)

  # Test overlapping matches
  expect_equal(find_motif("AAAAA", "AA"), c(1, 2, 3, 4))

})
