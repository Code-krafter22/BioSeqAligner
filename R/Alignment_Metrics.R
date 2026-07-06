#' @title Percent Identity of a Pairwise Alignment
#'
#' @description Fraction of aligned columns (excluding gap columns) at which the
#'   two sequences carry the same residue, expressed as a percentage.
#'
#' @param x An object of class \code{alignment} from \code{\link{align}}.
#'
#' @return A numeric percentage in \code{[0, 100]}.
#'
#' @examples
#' a <- align("ACGTGGA", "ACGTA")
#' alignment_identity(a)
#'
#' @export
alignment_identity <- function(x) {
  stopifnot(inherits(x, "alignment"))
  a1 <- strsplit(x$aligned1, "")[[1]]
  a2 <- strsplit(x$aligned2, "")[[1]]
  keep <- a1 != "-" & a2 != "-"
  if (!any(keep)) {
    return(0)
  }
  100 * sum(a1[keep] == a2[keep]) / sum(keep)
}

#' @title Alignment Summary Statistics
#'
#' @description Returns a one-row data frame of common alignment-quality
#'   metrics: length, number of matches / mismatches / gaps, percent identity,
#'   percent gaps, and query coverage.
#'
#' @param x An object of class \code{alignment}.
#'
#' @return A data frame with one row of metrics.
#'
#' @examples
#' alignment_stats(align("TTACGTGGATT", "ACGTGGA", method = "local"))
#'
#' @export
alignment_stats <- function(x) {
  stopifnot(inherits(x, "alignment"))
  a1 <- strsplit(x$aligned1, "")[[1]]
  a2 <- strsplit(x$aligned2, "")[[1]]
  aln_len <- length(a1)
  gap <- a1 == "-" | a2 == "-"
  matches <- sum(a1 == a2 & !gap)
  mismatches <- sum(a1 != a2 & !gap)
  gaps <- sum(gap)
  coverage <- 100 * (x$end1 - x$start1 + 1) / nchar(x$seq1)

  data.frame(
    method = x$method,
    score = x$score,
    length = aln_len,
    matches = matches,
    mismatches = mismatches,
    gaps = gaps,
    identity_pct = round(alignment_identity(x), 2),
    gap_pct = round(100 * gaps / aln_len, 2),
    query_coverage_pct = round(coverage, 2),
    stringsAsFactors = FALSE
  )
}
