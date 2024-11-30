#' @title Reverse Complement of a DNA Sequence
#'
#' @description This function computes the reverse complement of a given DNA sequence by reversing it and substituting complementary bases.
#'
#'
#' @param sequence A character string representing a DNA sequence. Must contain only A, C, G, and T.
#'
#' @return A character string representing the reverse complement of the input DNA sequence.
#'
#' @examples
#' sequence <- ("ATCGCTCA") # Returns "TGAGCGAT"
#'
#' @import stringi
#' @export
reverse_complement <- function(sequence) {
  # Check if the sequence is a DNA sequence
  if (!grepl("^[ACGT]*$", sequence)) {
    stop("Input is not a valid DNA sequence.")
  }

  # Reverse complement the sequence
  complement <- chartr("ACGT", "TGCA", sequence)
  reverse_complement <- stri_reverse(complement)

  return(reverse_complement)
}
#'
#' @title Transcribe DNA to RNA
#'
#' @description This function converts a DNA sequence into an RNA sequence by replacing thymine (T) with uracil (U).
#'
#' @param sequence A character string representing a DNA sequence. Must contain only A, C, G, and T.
#' @return A character string representing the transcribed RNA sequence.
#'
#' @examples
#' transcribe_dna("ATGC") # Returns "AUGC"
#' transcribe_dna("TTAA") # Returns "UUAA"
#'
#' @export
transcribe_dna <- function(sequence) {
  if (!grepl("^[ACGT]*$", sequence)) {
    stop("Input is not a valid DNA sequence.")
  }

  # Replace Thymine (T) with Uracil (U)
  rna_sequence <- chartr("T", "U", sequence)

  return(rna_sequence)
}
#'
#' @title Calculate GC Content of a DNA Sequence
#'
#' @description This function calculates the GC content (percentage of guanine and cytosine) in a given DNA sequence.
#'
#' @param sequence A character string representing a DNA sequence. Must contain only A, C, G, and T.
#' @return A numeric value representing the GC content percentage of the DNA sequence.
#'
#' @examples
#' gc_content("ATGC") # Returns 50
#' gc_content("GCGC") # Returns 100
#' gc_content("ATAT") # Returns 0
#'
#' @import stringr
#'
#' @export
gc_content <- function(sequence) {
  if (!grepl("^[ACGT]*$", sequence)) {
    stop("Input is not a valid DNA sequence.")
  }

  # Calculate GC content
  gc_count <- sum(str_count(sequence, "[GC]"))
  gc_percent <- (gc_count / nchar(sequence)) * 100

  return(gc_percent)
}
