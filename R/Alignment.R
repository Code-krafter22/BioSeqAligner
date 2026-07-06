#' @title Pairwise Sequence Alignment (Global or Local, Affine Gaps)
#'
#' @description
#' Aligns two sequences using dynamic programming. Supports global alignment
#' (Needleman-Wunsch) and local alignment (Smith-Waterman), both with affine
#' gap penalties (Gotoh algorithm). Works for DNA, RNA, and protein sequences
#' depending on the supplied scoring matrix.
#'
#' @param seq1 A character string, the first (query) sequence.
#' @param seq2 A character string, the second (reference) sequence.
#' @param method \code{"global"} (default) for Needleman-Wunsch or
#'   \code{"local"} for Smith-Waterman.
#' @param submat A substitution matrix from \code{\link{scoring_matrix}} or a
#'   compatible named numeric matrix. Defaults to the nucleotide matrix.
#' @param gap_open Penalty (a positive number) for opening a gap.
#' @param gap_extend Penalty (a positive number) for extending an existing gap.
#'
#' @return An object of class \code{alignment}: a list with the aligned
#'   strings (\code{aligned1}, \code{aligned2}), the alignment \code{score},
#'   \code{method}, and the start/end coordinates of the aligned region in
#'   each sequence. Has \code{print} and \code{plot} methods.
#'
#' @examples
#' a <- align("ACGTGGA", "ACGTA", method = "global")
#' print(a)
#'
#' # Local alignment of a short motif inside a longer sequence
#' align("TTACGTGGATT", "ACGTGGA", method = "local")
#'
#' @export
align <- function(seq1, seq2,
                  method = c("global", "local"),
                  submat = scoring_matrix("nucleotide"),
                  gap_open = 10, gap_extend = 1) {
  method <- match.arg(method)
  seq1 <- toupper(as.character(seq1))
  seq2 <- toupper(as.character(seq2))
  a <- strsplit(seq1, "")[[1]]
  b <- strsplit(seq2, "")[[1]]
  n <- length(a)
  m <- length(b)

  s <- function(x, y) {
    if (!(x %in% rownames(submat)) || !(y %in% colnames(submat))) {
      stop(sprintf("Residue '%s' or '%s' not found in scoring matrix.", x, y))
    }
    submat[x, y]
  }

  NEG <- -Inf
  # M = ends in a (mis)match; X = gap in seq2 (consumes seq1); Y = gap in seq1
  M <- matrix(NEG, n + 1, m + 1)
  X <- matrix(NEG, n + 1, m + 1)
  Y <- matrix(NEG, n + 1, m + 1)
  # Pointers: 1 = M, 2 = X, 3 = Y, 0 = stop (local, or origin)
  pM <- matrix(0L, n + 1, m + 1)
  pX <- matrix(0L, n + 1, m + 1)
  pY <- matrix(0L, n + 1, m + 1)

  if (method == "global") {
    M[1, 1] <- 0
    for (i in seq_len(n)) X[i + 1, 1] <- -(gap_open + (i - 1) * gap_extend)
    for (j in seq_len(m)) Y[1, j + 1] <- -(gap_open + (j - 1) * gap_extend)
  } else {
    # local: all borders start at 0 in M
    M[, 1] <- 0
    M[1, ] <- 0
  }

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      sc <- s(a[i], b[j])

      # X: gap in seq2 (vertical move, consume a[i])
      x_open <- M[i, j + 1] - gap_open
      x_ext <- X[i, j + 1] - gap_extend
      if (x_ext >= x_open) {
        X[i + 1, j + 1] <- x_ext
        pX[i + 1, j + 1] <- 2L
      } else {
        X[i + 1, j + 1] <- x_open
        pX[i + 1, j + 1] <- 1L
      }

      # Y: gap in seq1 (horizontal move, consume b[j])
      y_open <- M[i + 1, j] - gap_open
      y_ext <- Y[i + 1, j] - gap_extend
      if (y_ext >= y_open) {
        Y[i + 1, j + 1] <- y_ext
        pY[i + 1, j + 1] <- 3L
      } else {
        Y[i + 1, j + 1] <- y_open
        pY[i + 1, j + 1] <- 1L
      }

      # M: diagonal move from best of the three
      diag_prev <- c(M[i, j], X[i, j], Y[i, j])
      best <- which.max(diag_prev)
      mval <- diag_prev[best] + sc
      if (method == "local" && mval < 0) {
        M[i + 1, j + 1] <- 0
        pM[i + 1, j + 1] <- 0L
      } else {
        M[i + 1, j + 1] <- mval
        pM[i + 1, j + 1] <- best
      }
    }
  }

  # --- Find start cell for traceback ---
  if (method == "global") {
    finals <- c(M[n + 1, m + 1], X[n + 1, m + 1], Y[n + 1, m + 1])
    layer <- which.max(finals)
    score <- finals[layer]
    ci <- n + 1
    cj <- m + 1
  } else {
    score <- max(M)
    idx <- which(M == score, arr.ind = TRUE)[1, ]
    ci <- unname(idx[1])
    cj <- unname(idx[2])
    layer <- 1L
  }

  al1 <- character(0)
  al2 <- character(0)
  i <- ci
  j <- cj
  cur <- layer

  repeat {
    if (method == "local" && cur == 1L && M[i, j] == 0) break
    if (i == 1 && j == 1) break

    if (cur == 1L) { # match/mismatch -> diagonal
      al1 <- c(a[i - 1], al1)
      al2 <- c(b[j - 1], al2)
      nxt <- pM[i, j]
      i <- i - 1
      j <- j - 1
      cur <- if (nxt == 0L) 1L else nxt
    } else if (cur == 2L) { # gap in seq2 -> vertical
      al1 <- c(a[i - 1], al1)
      al2 <- c("-", al2)
      nxt <- pX[i, j]
      i <- i - 1
      cur <- nxt
    } else { # cur == 3L, gap in seq1 -> horizontal
      al1 <- c("-", al1)
      al2 <- c(b[j - 1], al2)
      nxt <- pY[i, j]
      j <- j - 1
      cur <- nxt
    }
  }

  # Remaining leading gaps for global alignment
  if (method == "global") {
    while (i > 1) {
      al1 <- c(a[i - 1], al1)
      al2 <- c("-", al2)
      i <- i - 1
    }
    while (j > 1) {
      al1 <- c("-", al1)
      al2 <- c(b[j - 1], al2)
      j <- j - 1
    }
  }

  start1 <- i
  start2 <- j

  result <- list(
    aligned1 = paste(al1, collapse = ""),
    aligned2 = paste(al2, collapse = ""),
    score = score,
    method = method,
    start1 = start1, end1 = ci - 1,
    start2 = start2, end2 = cj - 1,
    seq1 = seq1, seq2 = seq2
  )
  class(result) <- "alignment"
  result
}

#' @title Print an Alignment
#'
#' @description Pretty-prints an \code{alignment} object as a side-by-side
#'   block with a match line, plus score and percent identity.
#'
#' @param x An object of class \code{alignment}.
#' @param width Number of alignment columns per printed block.
#' @param ... Ignored.
#'
#' @return The alignment object, invisibly.
#' @export
print.alignment <- function(x, width = 60, ...) {
  a1 <- strsplit(x$aligned1, "")[[1]]
  a2 <- strsplit(x$aligned2, "")[[1]]
  midline <- ifelse(a1 == a2 & a1 != "-", "|", " ")
  ident <- alignment_identity(x)

  cat(sprintf(
    "%s alignment  |  score = %g  |  identity = %.1f%%\n\n",
    tools::toTitleCase(x$method), x$score, ident
  ))
  n <- length(a1)
  for (start in seq(1, max(1, n), by = width)) {
    end <- min(start + width - 1, n)
    cat("seq1  ", paste(a1[start:end], collapse = ""), "\n")
    cat("      ", paste(midline[start:end], collapse = ""), "\n")
    cat("seq2  ", paste(a2[start:end], collapse = ""), "\n\n")
  }
  invisible(x)
}
