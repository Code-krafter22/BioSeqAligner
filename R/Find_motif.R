#' Locates all positions of a motif in a sequence.
#'
#' @param seq DNA sequence (character string)
#' @param motif Motif to search for (character string)
#'
#' @return #A vector of starting positions where the motif occurs
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
