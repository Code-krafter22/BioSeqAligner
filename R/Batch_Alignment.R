#' @title Align One Query Against Many References (mini-BLAST)
#'
#' @description
#' Aligns a single query sequence against a set of reference sequences and
#' returns a ranked results table (best score first), giving a lightweight
#' BLAST-like search over a small, in-memory database.
#'
#' @param query A character string, the query sequence.
#' @param refs A named character vector of reference sequences (e.g. from
#'   \code{\link{read_fasta}}).
#' @param method \code{"local"} (default, good for finding subregions) or
#'   \code{"global"}.
#' @param submat Substitution matrix (see \code{\link{scoring_matrix}}).
#' @param gap_open,gap_extend Affine gap penalties.
#' @param top If not \code{NULL}, keep only the top \code{n} hits.
#'
#' @return A data frame ranked by score, with reference name, score, percent
#'   identity, alignment length, and the reference match coordinates.
#'
#' @examples
#' refs <- c(refA = "TTACGTGGATT", refB = "GGGGCCCC", refC = "AAACGTGGAAA")
#' batch_align("ACGTGGA", refs, method = "local")
#'
#' @export
batch_align <- function(query, refs,
                        method = c("local", "global"),
                        submat = scoring_matrix("nucleotide"),
                        gap_open = 10, gap_extend = 1, top = NULL) {
  method <- match.arg(method)
  if (is.null(names(refs))) names(refs) <- paste0("ref", seq_along(refs))

  rows <- lapply(names(refs), function(nm) {
    a <- align(query, refs[nm],
      method = method, submat = submat,
      gap_open = gap_open, gap_extend = gap_extend
    )
    st <- alignment_stats(a)
    data.frame(
      reference = nm,
      score = a$score,
      identity_pct = st$identity_pct,
      aln_length = st$length,
      ref_start = a$start2,
      ref_end = a$end2,
      stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, rows)
  res <- res[order(-res$score), ]
  rownames(res) <- NULL
  if (!is.null(top)) res <- utils::head(res, top)
  res
}
