test_that("find_cpg_islands detects a CpG-dense region flanked by AT-rich sequence", {
  cpg_core <- paste(rep("CG", 150), collapse = "")   # 300bp, GC-dense, high obs/exp
  at_flank <- paste(rep("AT", 100), collapse = "")   # 200bp AT-rich flank each side
  seq <- paste0(at_flank, cpg_core, at_flank)

  r <- find_cpg_islands(seq, window = 200, step = 5)
  expect_true(nrow(r) >= 1)
  # The detected island should overlap the true CpG-rich core (201-500),
  # allowing for the natural boundary padding of a sliding-window scan.
  expect_true(r$start[1] <= 300 && r$end[1] >= 400)
  expect_true(r$gc_percent[1] > 50)
  expect_true(r$obs_exp_ratio[1] > 0.6)
})

test_that("find_cpg_islands reports no islands in a uniformly AT-rich sequence", {
  at_only <- paste(rep("ATAT", 200), collapse = "")
  r <- find_cpg_islands(at_only, window = 200, step = 5)
  expect_equal(nrow(r), 0)
})

test_that("find_cpg_islands returns an empty data frame for sequences shorter than the window", {
  r <- find_cpg_islands("ACGTACGT", window = 200)
  expect_equal(nrow(r), 0)
  expect_equal(colnames(r), c("start", "end", "length", "gc_percent", "obs_exp_ratio"))
})

test_that("find_cpg_islands respects min_length", {
  cpg_core <- paste(rep("CG", 150), collapse = "")
  at_flank <- paste(rep("AT", 100), collapse = "")
  seq <- paste0(at_flank, cpg_core, at_flank)
  r <- find_cpg_islands(seq, window = 200, step = 5, min_length = 10000)
  expect_equal(nrow(r), 0)
})

test_that("cpg_disease_lookup returns a known curated entry", {
  r <- cpg_disease_lookup("MLH1")
  expect_equal(nrow(r), 1)
  expect_equal(r$gene, "MLH1")
  expect_true(nchar(r$disease) > 0)
  expect_true(nchar(r$mechanism) > 0)
})

test_that("cpg_disease_lookup is case-insensitive", {
  expect_equal(nrow(cpg_disease_lookup("mlh1")), 1)
  expect_equal(nrow(cpg_disease_lookup("Mlh1")), 1)
})

test_that("cpg_disease_lookup returns zero rows (with a message, not an error) for unknown genes", {
  expect_message(res <- cpg_disease_lookup("NOT_A_REAL_GENE_XYZ"), "not in this package's")
  expect_equal(nrow(res), 0)
})
