% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seq_Alignment_Visualization.R
\name{generate_dot_plot}
\alias{generate_dot_plot}
\title{Generate a Dot Plot for Sequence Alignment}
\usage{
generate_dot_plot(seq1, seq2)
}
\arguments{
\item{seq1}{A character string representing the first DNA sequence.}

\item{seq2}{A character string representing the second DNA sequence.}
}
\value{
A ggplot object displaying the dot plot for sequence alignment.
}
\description{
This function creates a dot plot comparing two DNA sequences, highlighting matching positions.
The dot plot visualizes sequence similarities and can be useful for identifying regions of alignment.
}
\examples{
seq1 <- "ATCGATCGATGC"
seq2 <- "ATCCGATCGTAT"
generate_dot_plot(seq1, seq2)
}
