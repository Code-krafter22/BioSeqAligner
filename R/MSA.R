#' @title Progressive Multiple Sequence Alignment
#'
#' @description
#' Aligns three or more sequences using a simple progressive strategy: the two
#' most similar sequences are aligned first, then remaining sequences are added
#' one at a time to the growing profile (guided by pairwise alignment scores).
#' This is a lightweight, dependency-free MSA suitable for teaching and
#' small-to-moderate sequence sets; for very large or divergent sets, wrap a
#' dedicated tool such as MUSCLE or Clustal.
#'
#' @param seqs A named character vector of sequences (e.g. from
#'   \code{\link{read_fasta}}). Names become row labels in the result.
#' @param submat A substitution matrix (see \code{\link{scoring_matrix}}).
#' @param gap_open,gap_extend Affine gap penalties passed to \code{\link{align}}.
#'
#' @return An object of class \code{msa}: a list with \code{aligned} (a named
#'   character vector of equal-length, gapped sequences) and \code{consensus}.
#'   Has \code{print} and \code{plot} methods.
#'
#' @examples
#' seqs <- c(s1 = "ACGTGGAA", s2 = "ACGTGCAA", s3 = "ACGTGGATA")
#' msa <- msa_align(seqs)
#' print(msa)
#'
#' @export
msa_align <- function(seqs,
                      submat = scoring_matrix("nucleotide"),
                      gap_open = 10, gap_extend = 1) {
  if (length(seqs) < 2) stop("Provide at least two sequences.")
  if (is.null(names(seqs))) names(seqs) <- paste0("seq", seq_along(seqs))
  seqs <- toupper(seqs)

  # Pairwise score matrix to pick the best starting pair / merge order
  n <- length(seqs)
  score <- matrix(-Inf, n, n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      score[i, j] <- align(seqs[i], seqs[j],
        method = "global",
        submat = submat, gap_open = gap_open, gap_extend = gap_extend
      )$score
      score[j, i] <- score[i, j]
    }
  }

  # Merge a new sequence into an aligned block by aligning it to the block's
  # consensus, then propagating any newly introduced gaps to all rows.
  insert_gaps <- function(block, gapped_consensus) {
    cons_chars <- strsplit(gapped_consensus, "")[[1]]
    gap_pos <- which(cons_chars == "-")
    if (length(gap_pos) == 0) {
      return(block)
    }
    vapply(block, function(s) {
      chars <- strsplit(s, "")[[1]]
      for (p in gap_pos) {
        chars <- append(chars, "-", after = p - 1)
      }
      paste(chars, collapse = "")
    }, character(1))
  }

  order_idx <- which(score == max(score), arr.ind = TRUE)[1, ]
  first_pair <- align(seqs[order_idx[1]], seqs[order_idx[2]],
    method = "global", submat = submat,
    gap_open = gap_open, gap_extend = gap_extend
  )
  block <- c(first_pair$aligned1, first_pair$aligned2)
  names(block) <- names(seqs)[order_idx]
  remaining <- setdiff(seq_len(n), order_idx)

  while (length(remaining) > 0) {
    cons <- .msa_consensus(block)
    # Add whichever remaining sequence aligns best to the current consensus
    cand_scores <- vapply(remaining, function(k) {
      align(cons, seqs[k],
        method = "global", submat = submat,
        gap_open = gap_open, gap_extend = gap_extend
      )$score
    }, numeric(1))
    pick <- remaining[which.max(cand_scores)]

    al <- align(cons, seqs[pick],
      method = "global", submat = submat,
      gap_open = gap_open, gap_extend = gap_extend
    )
    block <- insert_gaps(block, al$aligned1)
    new_row <- al$aligned2
    names(new_row) <- names(seqs)[pick]
    # Pad any length mismatch (defensive) so all rows stay equal length
    L <- max(nchar(block), nchar(new_row))
    block <- vapply(block, function(s) formatC(s, width = -L, flag = "-"), character(1))
    block <- gsub(" ", "-", block)
    new_row <- gsub(" ", "-", formatC(new_row, width = -L, flag = "-"))
    block <- c(block, new_row)

    remaining <- setdiff(remaining, pick)
  }

  block <- block[names(seqs)] # restore input order
  result <- list(aligned = block, consensus = .msa_consensus(block))
  class(result) <- "msa"
  result
}

# Column-wise majority consensus of an aligned block (internal helper)
.msa_consensus <- function(block) {
  mat <- do.call(rbind, strsplit(block, ""))
  cons <- apply(mat, 2, function(col) {
    col <- col[col != "-"]
    if (length(col) == 0) {
      return("-")
    }
    names(sort(table(col), decreasing = TRUE))[1]
  })
  paste(cons, collapse = "")
}

#' @title Per-column Conservation Scores of an MSA
#'
#' @description Computes, for each alignment column, the fraction of rows that
#'   carry the most common residue in that column (gaps excluded). A value of
#'   1 means fully conserved.
#'
#' @param x An object of class \code{msa}.
#'
#' @return A numeric vector, one conservation score per column.
#'
#' @examples
#' conservation_scores(msa_align(c(a = "ACGT", b = "ACGA", c = "ACGT")))
#'
#' @export
conservation_scores <- function(x) {
  stopifnot(inherits(x, "msa"))
  mat <- do.call(rbind, strsplit(x$aligned, ""))
  apply(mat, 2, function(col) {
    col <- col[col != "-"]
    if (length(col) == 0) {
      return(0)
    }
    max(table(col)) / length(col)
  })
}

#' @title Print an MSA
#'
#' @description Prints the aligned rows plus a consensus line.
#'
#' @param x An object of class \code{msa}.
#' @param width Columns per printed block.
#' @param ... Ignored.
#'
#' @return The msa object, invisibly.
#' @export
print.msa <- function(x, width = 60, ...) {
  aln <- x$aligned
  label_w <- max(nchar(names(aln)), nchar("consensus"))
  n <- nchar(aln[1])
  cat(sprintf("Multiple sequence alignment: %d sequences x %d columns\n\n",
              length(aln), n))
  for (start in seq(1, n, by = width)) {
    end <- min(start + width - 1, n)
    for (nm in names(aln)) {
      cat(formatC(nm, width = -label_w), " ",
          substring(aln[nm], start, end), "\n", sep = "")
    }
    cat(formatC("consensus", width = -label_w), " ",
        substring(x$consensus, start, end), "\n\n", sep = "")
  }
  invisible(x)
}
