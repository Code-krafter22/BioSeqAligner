test_that("bisulfite_convert matches a hand-computed example", {
  expect_equal(bisulfite_convert("ACGTCCGACG"), "ATGTTTGATG")
  expect_equal(bisulfite_convert("ACGTCCGACG", assume_methylated_cpg = TRUE), "ACGTTCGACG")
})

test_that("fully unmethylated conversion leaves zero cytosines", {
  seq <- "CGCGATCGATCGCCGGATCGATCGCCCGATCG"
  converted <- bisulfite_convert(seq)
  expect_equal(nchar(gsub("[^C]", "", converted)), 0)
  expect_equal(nchar(converted), nchar(seq))
})

test_that("CpG-protected conversion retains exactly the original CpG dinucleotide count", {
  seq <- "CGCGATCGATCGCCGGATCGATCGCCCGATCG"
  converted <- bisulfite_convert(seq, assume_methylated_cpg = TRUE)
  cg_count <- lengths(regmatches(seq, gregexpr("CG", seq)))
  c_remaining <- nchar(gsub("[^C]", "", converted))
  expect_equal(c_remaining, cg_count)
  expect_equal(nchar(converted), nchar(seq))
})

test_that("non-C bases are never altered", {
  seq <- "AGGTAGGTAGGT"  # no cytosines at all
  expect_equal(bisulfite_convert(seq), seq)
  expect_equal(bisulfite_convert(seq, assume_methylated_cpg = TRUE), seq)
})

test_that("design_bisulfite_primers produces different primers for the two states", {
  skip_if_not(nzchar(Sys.which("primer3_core")), "primer3_core not installed on this system")
  seq <- paste0(
    "CCACGTGAGGTGCATCCATTCACCATTCTCGCTGACAACGTTACTCCGCGTTTGAGGAAT",
    "GCGTCAACGAAAATAGATTACGCACTGCCGCACTGCCGCACTGCCGCACTGCCGCACTGC",
    "CGCACTGCCGCACTGCCGCACTGCCGCACTGCCGCACTGCCGCACTGCCGCACTGCTATG",
    "GTAAATGAACCGTTGGGAATCCGGTAGCGTTTATGCTTTTGTCCAGCGGCCTCAGGAATG",
    "GCACAAGTGTGGAAAG"
  )
  r_u <- design_bisulfite_primers(seq, "unmethylated", product_size_range = c(70, 180), num_return = 2)
  r_m <- design_bisulfite_primers(seq, "methylated", product_size_range = c(70, 180), num_return = 2)

  expect_true(nrow(r_u) > 0)
  expect_true(nrow(r_m) > 0)
  expect_true(all(r_u$methylation_state == "unmethylated"))
  expect_true(all(r_m$methylation_state == "methylated"))
  expect_false(identical(r_u$left_seq, r_m$left_seq))
})
