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
#' @title Create a DNASequence Object
#'
#' @description Constructs an S3 object of class \code{DNASequence} representing a DNA sequence.
#' This constructor validates the input to ensure it contains only valid DNA bases (A, C, G, T).
#'
#' @param sequence A character string representing a DNA sequence. Must contain only A, C, G, and T.
#'
#' @return An object of class \code{DNASequence} with the sequence attribute.
#'
#' @examples
#' dna_seq <- DNASequence("ATGCGC")
#' print(dna_seq)
#'
#' @export
DNASequence <- function(sequence) {
  if (!grepl("^[ACGT]*$", sequence)) {
    stop("Input is not a valid DNA sequence.")
  }
  structure(sequence, class = "DNASequence", sequence = sequence)
}

#' @title Calculate GC Content of a DNA Sequence (S3 Method)
#'
#' @description Calculates the GC content (percentage of guanine and cytosine) in a DNA sequence represented by a \code{DNASequence} object.
#'
#' @param sequence An object of class \code{DNASequence}.
#'
#' @return A numeric value representing the GC content percentage of the DNA sequence.
#'
#' @examples
#' dna_seq <- DNASequence("ATGCGC")
#' gc_content.S3(dna_seq)  # Returns 66.67
#'
#' @import stringr
#'
#' @export
gc_content.S3 <- function(sequence) {
  if (!inherits(sequence, "DNASequence")) {
    stop("Input must be of class 'DNASequence'.")
  }
  sequence_value <- attr(sequence, "sequence")
  gc_count <- sum(str_count(sequence_value, "[GC]"))
  gc_percent <- (gc_count / nchar(sequence_value)) * 100
  return(gc_percent)
}

