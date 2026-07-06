test_that("FASTA read handles multi-line records and round-trips", {
  fa <- tempfile(fileext = ".fasta")
  writeLines(c(">seqA", "ACGT", "GGCA", ">seqB", "TTTT"), fa)
  r <- read_fasta(fa)
  expect_equal(unname(r["seqA"]), "ACGTGGCA")
  expect_equal(unname(r["seqB"]), "TTTT")

  out <- tempfile(fileext = ".fasta")
  write_fasta(r, out)
  expect_equal(read_fasta(out), r)
})

test_that("FASTQ read parses id, sequence, and quality", {
  fq <- tempfile(fileext = ".fastq")
  writeLines(c("@r1", "ACGT", "+", "IIII", "@r2", "TTGG", "+", "JJJJ"), fq)
  q <- read_fastq(fq)
  expect_equal(nrow(q), 2)
  expect_equal(q$sequence[1], "ACGT")
  expect_equal(q$quality[2], "JJJJ")
})

test_that("sequence utilities behave as expected", {
  expect_true(validate_sequence("ACGTN", "iupac"))
  expect_error(validate_sequence("ACGTZ", "dna"))
  expect_equal(unname(kmer_count("ATATAT", 2)["AT"]), 3L)
  expect_equal(sum(codon_usage("ATGAAATAG")), 3L)
  expect_equal(round(seq_composition("ATGCGC")$gc_percent), 67)
})

test_that("batch_align ranks the best reference first", {
  refs <- c(refA = "TTACGTGGATT", refB = "GGGGCCCC", refC = "AAACGTGGAAA")
  b <- batch_align("ACGTGGA", refs, method = "local")
  expect_equal(b$reference[1], "refA")
  expect_true(all(diff(b$score) <= 0))
})
