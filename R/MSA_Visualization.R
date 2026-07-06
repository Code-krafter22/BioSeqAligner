#' @title Plot a Pairwise Alignment
#'
#' @description ggplot2 tile view of a pairwise alignment, colouring matches,
#'   mismatches, and gaps so regions of similarity stand out.
#'
#' @param x An object of class \code{alignment}.
#' @param ... Ignored.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' plot(align("ACGTGGA", "ACGTGCA"))
#' }
#'
#' @import ggplot2
#' @export
plot.alignment <- function(x, ...) {
  a1 <- strsplit(x$aligned1, "")[[1]]
  a2 <- strsplit(x$aligned2, "")[[1]]
  n <- length(a1)

  status <- ifelse(a1 == "-" | a2 == "-", "Gap",
    ifelse(a1 == a2, "Match", "Mismatch")
  )
  df <- data.frame(
    pos = rep(seq_len(n), 2),
    row = rep(c("seq1", "seq2"), each = n),
    residue = c(a1, a2),
    status = rep(status, 2),
    stringsAsFactors = FALSE
  )
  df$row <- factor(df$row, levels = c("seq2", "seq1"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pos, y = .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$status), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$residue), size = 3) +
    ggplot2::scale_fill_manual(values = c(
      Match = "#8ecae6", Mismatch = "#ffb703", Gap = "#e5e5e5"
    )) +
    ggplot2::labs(
      title = sprintf("%s alignment (identity %.1f%%)",
                      tools::toTitleCase(x$method), alignment_identity(x)),
      x = "Alignment position", y = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' @title Plot a Multiple Sequence Alignment
#'
#' @description ggplot2 tile view of an MSA, colouring each residue and
#'   optionally shading columns by conservation.
#'
#' @param x An object of class \code{msa}.
#' @param color_by \code{"residue"} (default) colours each nucleotide/amino
#'   acid; \code{"conservation"} shades columns by how conserved they are.
#' @param ... Ignored.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' m <- msa_align(c(a = "ACGTGGAA", b = "ACGTGCAA", c = "ACGTGGATA"))
#' plot(m)
#' }
#'
#' @import ggplot2
#' @export
plot.msa <- function(x, color_by = c("residue", "conservation"), ...) {
  color_by <- match.arg(color_by)
  aln <- x$aligned
  mat <- do.call(rbind, strsplit(aln, ""))
  rownames(mat) <- names(aln)
  nrow_a <- nrow(mat)
  ncol_a <- ncol(mat)

  # as.vector(mat) fills column-major, so the sequence label varies fastest
  df <- data.frame(
    row = factor(rep(names(aln), times = ncol_a), levels = rev(names(aln))),
    pos = rep(seq_len(ncol_a), each = nrow_a),
    residue = as.vector(mat),
    stringsAsFactors = FALSE
  )

  base <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pos, y = .data$row))

  if (color_by == "residue") {
    p <- base +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$residue), colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = .data$residue), size = 2.6) +
      ggplot2::scale_fill_manual(values = .residue_palette(unique(df$residue)))
  } else {
    cons <- conservation_scores(x)
    df$conservation <- cons[df$pos]
    p <- base +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$conservation), colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = .data$residue), size = 2.6) +
      ggplot2::scale_fill_gradient(low = "#f1faee", high = "#e63946",
                                   limits = c(0, 1), name = "Conservation")
  }

  p +
    ggplot2::labs(
      title = "Multiple Sequence Alignment",
      x = "Alignment column", y = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# Colour palette for residues (nucleotides + gap; falls back to grey)
.residue_palette <- function(residues) {
  pal <- c(
    A = "#8ecae6", C = "#ffb703", G = "#90be6d", T = "#f28482",
    U = "#f28482", N = "#cccccc", "-" = "#eeeeee"
  )
  missing <- setdiff(residues, names(pal))
  if (length(missing) > 0) {
    extra <- grDevices::rainbow(length(missing))
    names(extra) <- missing
    pal <- c(pal, extra)
  }
  pal[as.character(residues)]
}
