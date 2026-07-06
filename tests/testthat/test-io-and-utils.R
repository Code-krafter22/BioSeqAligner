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

test_that("parse_fasta never touches the filesystem, even if text matches a real path", {
  real_file <- tempfile(fileext = ".fasta")
  writeLines(c(">from_disk", "AAAA"), real_file)

  # A single-line "path-like" string with no FASTA header must error, not
  # silently read the file at that path (this was the path-traversal bug:
  # read_fasta() used to fall back to treating such input as a file path).
  expect_error(parse_fasta(real_file), "No FASTA header")

  # Confirm it truly never opened the file: a path-like string appearing
  # INSIDE FASTA text is treated as literal sequence text, never resolved
  # against the filesystem (if it had been opened, this would contain
  # "AAAA" from disk instead of the path string itself).
  pasted <- paste0(">pasted_header\n", real_file, "\nACGT")
  result <- parse_fasta(pasted)
  expect_equal(unname(result["pasted_header"]), toupper(paste0(real_file, "ACGT")))
})

test_that("read_fasta rejects non-existent paths instead of treating them as content", {
  expect_error(read_fasta(">seqA\nACGT"), "existing file")
  expect_error(read_fasta("/definitely/not/a/real/path.fasta"), "existing file")
})

test_that("read_fastq rejects non-existent paths instead of treating them as content", {
  expect_error(read_fastq("@r1\nACGT\n+\nIIII"), "existing file")
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
