test_that("msa_align returns equal-length rows in input order", {
  seqs <- c(s1 = "ACGTGGAA", s2 = "ACGTGCAA", s3 = "ACGTGGATA")
  m <- msa_align(seqs)
  expect_s3_class(m, "msa")
  expect_length(unique(nchar(m$aligned)), 1)
  expect_equal(names(m$aligned), names(seqs))
  expect_equal(nchar(m$consensus), unname(nchar(m$aligned[1])))
})

test_that("conservation scores lie in [0, 1] with one per column", {
  m <- msa_align(c(a = "ACGT", b = "ACGA", c = "ACGT"))
  cs <- conservation_scores(m)
  expect_length(cs, nchar(m$aligned[1]))
  expect_true(all(cs >= 0 & cs <= 1))
})

test_that("msa_distance is symmetric with a zero diagonal", {
  m <- msa_align(c(a = "ACGTGGAA", b = "ACGTGCAA", c = "TTGTGGATA"))
  d <- msa_distance(m)
  expect_true(all(diag(d) == 0))
  expect_equal(d, t(d))
})
