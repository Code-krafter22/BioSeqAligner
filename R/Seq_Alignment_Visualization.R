#' @title Generate a Dot Plot for Sequence Alignment
#'
#' @description
#' This function creates a dot plot comparing two DNA sequences, highlighting matching positions.
#' The dot plot visualizes sequence similarities and can be useful for identifying regions of alignment.
#
#' @param seq1 A character string representing the first DNA sequence.
#' @param seq2  A character string representing the second DNA sequence.
#'
#' @return A ggplot object displaying the dot plot for sequence alignment.
#' @import ggplot2
#'
#' @examples
#' seq1 <- "ATCGATCGATGC"
#' seq2 <- "ATCCGATCGTAT"
#' generate_dot_plot(seq1, seq2)
#' @export
generate_dot_plot <- function(seq1, seq2) {
  # Create matrix to store matches
  dot_plot_matrix <- matrix(0, nrow = nchar(seq1), ncol = nchar(seq2))

  # Fill matrix with match indicators
  for (i in 1:nchar(seq1)) {
    for (j in 1:nchar(seq2)) {
      if (substr(seq1, i, i) == substr(seq2, j, j)) {
        dot_plot_matrix[i, j] <- 1
      }
    }
  }

  # Convert matrix into data frame
  dot_plot_data <- expand.grid(Seq1_Pos = 1:nchar(seq1), Seq2_Pos = 1:nchar(seq2))
  dot_plot_data$Match <- as.vector(dot_plot_matrix)

  # Plot with ggplot2 using namespace-safe syntax
  ggplot2::ggplot(dot_plot_data, ggplot2::aes(x = .data$Seq1_Pos, y = .data$Seq2_Pos)) +
    ggplot2::geom_point(ggplot2::aes(color = factor(.data$Match)), size = 3) +
    ggplot2::scale_color_manual(values = c("0" = "lightblue", "1" = "darkred")) +
    ggplot2::labs(
      title = "Dot Plot for Sequence Alignment",
      x = "Seq1 Position",
      y = "Seq2 Position"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}
