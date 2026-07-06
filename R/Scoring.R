#' @title Build a Substitution (Scoring) Matrix
#'
#' @description
#' Constructs a named symmetric substitution matrix used to score aligned
#' residues. Two ready-made schemes are provided: a simple nucleotide matrix
#' (match / mismatch) for DNA/RNA, and the BLOSUM62 matrix for proteins.
#' You can also supply your own \code{match} / \code{mismatch} values.
#'
#' @param type One of \code{"nucleotide"} (default) or \code{"BLOSUM62"}.
#' @param match Score for identical residues (nucleotide scheme only).
#' @param mismatch Score for differing residues (nucleotide scheme only).
#'
#' @return A numeric matrix with row/column names equal to the residue
#'   alphabet, suitable for passing to \code{\link{align}}.
#'
#' @examples
#' m <- scoring_matrix("nucleotide", match = 2, mismatch = -1)
#' m["A", "A"]  # 2
#' m["A", "G"]  # -1
#' p <- scoring_matrix("BLOSUM62")
#' p["W", "W"]  # 11
#'
#' @export
scoring_matrix <- function(type = c("nucleotide", "BLOSUM62"),
                           match = 2, mismatch = -1) {
  type <- match.arg(type)

  if (type == "nucleotide") {
    alphabet <- c("A", "C", "G", "T", "U", "N")
    m <- matrix(mismatch,
      nrow = length(alphabet), ncol = length(alphabet),
      dimnames = list(alphabet, alphabet)
    )
    diag(m) <- match
    # Treat T and U as equivalent so DNA can be scored against RNA
    m["T", "U"] <- match
    m["U", "T"] <- match
    # N (unknown) is neutral against anything
    m["N", ] <- 0
    m[, "N"] <- 0
    return(m)
  }

  # BLOSUM62 (standard 20 amino acids + B, Z, X, *)
  aa <- c(
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
    "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*"
  )
  vals <- c(
     4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4,
    -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4,
    -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4,
    -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4,
     0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
    -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4,
    -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
     0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4,
    -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4,
    -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4,
    -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4,
    -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4,
    -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4,
    -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4,
    -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
     1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4,
     0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4,
    -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4,
    -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4,
     0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4,
    -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4,
    -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
     0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1
  )
  m <- matrix(vals, nrow = length(aa), byrow = TRUE, dimnames = list(aa, aa))
  m
}
