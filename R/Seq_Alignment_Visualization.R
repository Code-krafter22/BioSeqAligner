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
#' seq1<- seq1 <- "ATCGATCGATGC"
#' seq2 <- "ATCCGATCGTAT"
#' generate_dot_plot(seq1, seq2)
#' @export
generate_dot_plot <- function(seq1, seq2) {
# Create a matrix to store dot plot data
dot_plot_matrix <- matrix(0, nrow = nchar(seq1), ncol = nchar(seq2))

# Compare characters of seq1 and seq2
for (i in 1:nchar(seq1)) {
  for (j in 1:nchar(seq2)) {
    if (substr(seq1, i, i) == substr(seq2, j, j)) {
      dot_plot_matrix[i, j] <- 1  # Mark a match with a 1 (dot)
    }
  }
}

# Plot the dot plot using ggplot2
library(ggplot2)

# Convert the matrix into a data frame for ggplot
dot_plot_data <- expand.grid(Seq1_Pos = 1:nchar(seq1), Seq2_Pos = 1:nchar(seq2))
dot_plot_data$Match <- as.vector(dot_plot_matrix)

# Plot the dot plot
ggplot(dot_plot_data, aes(x = Seq1_Pos, y = Seq2_Pos)) +
  geom_point(aes(color = factor(Match)), size = 3) +
  scale_color_manual(values = c("0" = "lightblue", "1" = "darkred")) +  # White for no match, blue for match
  labs(title = "Dot Plot for Sequence Alignment", x = "Seq1 Position", y = "Seq2 Position") +
  theme_minimal() +
  theme(legend.position = "none")
}
