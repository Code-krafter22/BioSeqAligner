#' @title Locates all positions of a motif in a sequence.
#'
#' @description
#' This function identifies all starting positions of a specified motif within a given DNA sequence.
#' It performs a linear scan and returns the positions where the motif matches the sequence.
#'
#' @param seq A character string representing the DNA sequence to search.
#' @param motif A character string representing the motif to search for within the sequence.
#'
#' @return #A vector of starting positions where the motif occurs
#' If no matches are found, it returns an empty vector.
#'
#' #' @examples
#' find_motif("ATCGATCGATGC", "ATC")  # Returns positions 1 and 5
#' find_motif("GATCGATCGTAT", "GAT")  # Returns positions 1 and 6
#' find_motif("AAAAAA", "TT")         # Returns an empty vector
#'
#' @export
#'
find_motif <- function(seq, motif) {
  seq_len <- nchar(seq)
  motif_len <- nchar(motif)
  positions <- c()
  for (i in 1:(seq_len - motif_len + 1)) {
    if (substr(seq, i, i + motif_len - 1) == motif) positions <- c(positions, i)
  }
  positions
}
