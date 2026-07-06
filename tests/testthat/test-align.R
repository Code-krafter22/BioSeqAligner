test_that("global alignment returns equal-length gapped strings", {
  a <- align("ACGTGGA", "ACGTA", method = "global")
  expect_s3_class(a, "alignment")
  expect_equal(nchar(a$aligned1), nchar(a$aligned2))
})

test_that("identical sequences give 100% identity", {
  a <- align("ACGTACGT", "ACGTACGT")
  expect_equal(alignment_identity(a), 100)
})

test_that("local alignment locates a subregion with correct coordinates", {
  a <- align("TTACGTGGATT", "ACGTGGA", method = "local")
  expect_equal(a$method, "local")
  expect_equal(a$start1, 3)
  expect_equal(a$end1, 9)
  expect_equal(gsub("-", "", a$aligned2), "ACGTGGA")
})

test_that("protein alignment with BLOSUM62 scores the diagonal", {
  mp <- scoring_matrix("BLOSUM62")
  a <- align("MKVLA", "MKVLA", method = "global", submat = mp)
  expect_equal(a$score, sum(diag(mp[c("M", "K", "V", "L", "A"),
                                    c("M", "K", "V", "L", "A")])))
})

test_that("alignment_stats returns one row with expected coverage", {
  a <- align("TTACGTGGATT", "ACGTGGA", method = "local")
  st <- alignment_stats(a)
  expect_equal(nrow(st), 1)
  expect_equal(st$query_coverage_pct, round(100 * 7 / 11, 2))
})
