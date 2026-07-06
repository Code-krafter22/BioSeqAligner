skip_if_not(nzchar(Sys.which("primer3_core")), "primer3_core not installed on this system")

target <- paste0(
  "ATGGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGC",
  "ATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGC"
)

test_that("design_primers returns well-formed primer pairs", {
  p <- design_primers(target, product_size_range = c(70, 120), num_return = 3)
  expect_s3_class(p, "data.frame")
  expect_equal(nrow(p), 3)
  expect_true(all(c("left_seq", "right_seq", "left_tm", "right_tm", "penalty") %in% colnames(p)))
  expect_true(all(nchar(p$left_seq) > 0))
})

test_that("design_and_validate_primers flags a truly distinct off-target as specific", {
  distinct <- strrep("GGGGGGCCCCCC", 10)
  offs <- c(GENE_A = target, UNRELATED = distinct)
  r <- design_and_validate_primers(target, offs,
    intended_target = "GENE_A", product_size_range = c(70, 120), num_return = 2
  )
  expect_true(all(r$specific))
  expect_true(all(r$max_offtarget_identity == 0))
})

test_that("design_and_validate_primers flags an identical paralog as not specific", {
  offs <- c(GENE_A = target, PARALOG_B = target)
  r <- design_and_validate_primers(target, offs,
    intended_target = "GENE_A", product_size_range = c(70, 120), num_return = 2
  )
  expect_true(all(!r$specific))
  expect_true(all(r$max_offtarget_identity == 100))
})

test_that("a bad primer3_path errors clearly", {
  expect_error(
    design_primers(target, primer3_path = "/definitely/not/here"),
    "does not exist"
  )
})
