#' @title Distance Matrix from a Multiple Sequence Alignment
#'
#' @description Computes a pairwise distance matrix from an MSA using
#'   \emph{p-distance} (proportion of differing, non-gap columns between each
#'   pair of sequences). This is the standard hand-off point to phylogenetics.
#'
#' @param x An object of class \code{msa}.
#'
#' @return A symmetric \code{dist}-compatible numeric matrix of pairwise
#'   distances, with sequence names as row/column labels.
#'
#' @examples
#' m <- msa_align(c(a = "ACGTGGAA", b = "ACGTGCAA", c = "TTGTGGATA"))
#' msa_distance(m)
#'
#' @export
msa_distance <- function(x) {
  stopifnot(inherits(x, "msa"))
  mat <- do.call(rbind, strsplit(x$aligned, ""))
  n <- nrow(mat)
  d <- matrix(0, n, n, dimnames = list(rownames(mat), rownames(mat)))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      keep <- mat[i, ] != "-" & mat[j, ] != "-"
      if (!any(keep)) {
        dist_ij <- 1
      } else {
        dist_ij <- mean(mat[i, keep] != mat[j, keep])
      }
      d[i, j] <- dist_ij
      d[j, i] <- dist_ij
    }
  }
  d
}

#' @title Build and Plot a Neighbour-Joining Tree from an MSA
#'
#' @description Computes a distance matrix from an MSA and builds a
#'   neighbour-joining phylogenetic tree using the \pkg{ape} package. This
#'   bridges alignment straight into phylogenetics in one call.
#'
#' @param x An object of class \code{msa}.
#' @param plot Logical; if \code{TRUE} (default) the tree is plotted.
#'
#' @return An object of class \code{phylo} (from \pkg{ape}), invisibly.
#'
#' @examples
#' \donttest{
#' m <- msa_align(c(a = "ACGTGGAA", b = "ACGTGCAA", c = "TTGTGGATA", d = "ACGTGGTA"))
#' nj_tree(m)
#' }
#'
#' @export
nj_tree <- function(x, plot = TRUE) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for nj_tree(). Install it with install.packages('ape').")
  }
  d <- stats::as.dist(msa_distance(x))
  tree <- ape::nj(d)
  if (plot) {
    ape::plot.phylo(tree, main = "Neighbour-Joining Tree", cex = 0.9)
  }
  invisible(tree)
}
