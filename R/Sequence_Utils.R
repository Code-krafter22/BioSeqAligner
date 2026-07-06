#' @title Validate a Biological Sequence
#'
#' @description Checks that a sequence contains only characters allowed for the
#'   given alphabet. Unlike the strict DNA-only helpers, this supports DNA, RNA,
#'   protein, and IUPAC ambiguity codes.
#'
#' @param sequence A character string.
#' @param type One of \code{"dna"}, \code{"rna"}, \code{"protein"}, or
#'   \code{"iupac"} (nucleotide ambiguity codes).
#'
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#'
#' @examples
#' validate_sequence("ACGTN", "iupac")
#' validate_sequence("MKV", "protein")
#'
#' @export
validate_sequence <- function(sequence, type = c("dna", "rna", "protein", "iupac")) {
  type <- match.arg(type)
  patt <- switch(type,
    dna = "^[ACGT]*$",
    rna = "^[ACGU]*$",
    protein = "^[ACDEFGHIKLMNPQRSTVWYBZX*]*$",
    iupac = "^[ACGTURYSWKMBDHVN]*$"
  )
  if (!grepl(patt, toupper(sequence))) {
    stop(sprintf("Input is not a valid %s sequence.", toupper(type)))
  }
  invisible(TRUE)
}

#' @title Count k-mers in a Sequence
#'
#' @description Counts all overlapping substrings of length \code{k}.
#'
#' @param sequence A character string.
#' @param k k-mer length (positive integer).
#'
#' @return A named integer vector of k-mer counts, sorted in decreasing order.
#'
#' @examples
#' kmer_count("ATATAT", 2)  # AT = 3, TA = 2
#'
#' @export
kmer_count <- function(sequence, k = 3) {
  sequence <- toupper(as.character(sequence))
  n <- nchar(sequence)
  if (k < 1 || k > n) stop("k must be between 1 and the sequence length.")
  starts <- seq_len(n - k + 1)
  kmers <- substring(sequence, starts, starts + k - 1)
  sort(table(kmers), decreasing = TRUE)
}

#' @title Codon Usage of a Coding Sequence
#'
#' @description Splits a sequence into non-overlapping codons (triplets) and
#'   tabulates their frequencies. Trailing bases that do not form a full codon
#'   are ignored.
#'
#' @param sequence A character string (a coding DNA/RNA sequence).
#'
#' @return A named integer vector of codon counts, sorted decreasing.
#'
#' @examples
#' codon_usage("ATGAAATAG")
#'
#' @export
codon_usage <- function(sequence) {
  sequence <- toupper(as.character(sequence))
  n <- nchar(sequence)
  ncod <- n %/% 3
  if (ncod == 0) stop("Sequence is shorter than one codon.")
  starts <- seq(1, by = 3, length.out = ncod)
  codons <- substring(sequence, starts, starts + 2)
  sort(table(codons), decreasing = TRUE)
}

#' @title Base / Residue Composition
#'
#' @description Tabulates the frequency of each residue in a sequence and, for
#'   nucleotides, reports GC content.
#'
#' @param sequence A character string.
#'
#' @return A list with element \code{counts} (named integer vector) and, when
#'   the sequence looks like a nucleotide sequence, \code{gc_percent}.
#'
#' @examples
#' seq_composition("ATGCGC")
#'
#' @export
seq_composition <- function(sequence) {
  sequence <- toupper(as.character(sequence))
  chars <- strsplit(sequence, "")[[1]]
  counts <- table(chars)
  out <- list(counts = counts)
  if (all(names(counts) %in% c("A", "C", "G", "T", "U", "N"))) {
    gc <- sum(counts[names(counts) %in% c("G", "C")])
    out$gc_percent <- 100 * gc / nchar(sequence)
  }
  out
}
