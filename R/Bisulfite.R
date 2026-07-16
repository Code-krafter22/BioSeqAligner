#' @title Simulate Bisulfite Conversion of a DNA Sequence
#'
#' @description
#' Simulates sodium bisulfite treatment of genomic DNA, the standard
#' first step for methylation-specific PCR (MSP) and bisulfite sequencing
#' PCR (BSP) primer design. Bisulfite treatment converts unmethylated
#' cytosine to uracil (read as thymine after PCR); 5-methylcytosine is
#' chemically protected and remains cytosine.
#'
#' Since non-CpG cytosine methylation is negligible in mammalian genomic
#' DNA, this function uses the standard simplifying assumption: every
#' non-CpG cytosine always converts, and CpG-context cytosines convert only
#' if \code{assume_methylated_cpg = FALSE} (the default, simulating a fully
#' unmethylated template). Set \code{assume_methylated_cpg = TRUE} to
#' simulate a template where every CpG site is methylated (protected).
#' Real samples are usually a mix of methylated and unmethylated CpGs
#' across cells/alleles -- these two calls give you the two extreme
#' reference states MSP primer design is built around, not a prediction of
#' your sample's actual methylation level.
#'
#' CpG context is always determined from the original, unconverted
#' sequence, so conversion is a single deterministic pass with no
#' order-dependent effects.
#'
#' @param sequence A character string, the genomic DNA sequence (top/sense
#'   strand) to convert.
#' @param assume_methylated_cpg If \code{FALSE} (default), every cytosine
#'   converts to thymine (fully unmethylated template). If \code{TRUE},
#'   cytosines in a CpG context (\code{C} immediately followed by \code{G})
#'   are left as \code{C}; all other cytosines still convert (fully
#'   CpG-methylated template).
#'
#' @return A character string, the bisulfite-converted sequence.
#'
#' @examples
#' bisulfite_convert("ACGTCCGACG")                          # "ATGTTTGATG"
#' bisulfite_convert("ACGTCCGACG", assume_methylated_cpg = TRUE)  # "ACGTTCGACG"
#'
#' @export
bisulfite_convert <- function(sequence, assume_methylated_cpg = FALSE) {
  sequence <- toupper(as.character(sequence))
  chars <- strsplit(sequence, "")[[1]]
  n <- length(chars)
  is_c <- chars == "C"

  if (assume_methylated_cpg) {
    is_cpg_c <- c(is_c[-n] & chars[-1] == "G", FALSE)
    convert <- is_c & !is_cpg_c
  } else {
    convert <- is_c
  }

  chars[convert] <- "T"
  paste(chars, collapse = "")
}

#' @title Design Primers Against a Bisulfite-Converted Template
#'
#' @description
#' Runs \code{\link{bisulfite_convert}} on \code{sequence} to simulate one
#' of the two reference bisulfite-treated states, then runs
#' \code{\link{design_primers}} on the converted template -- the standard
#' workflow for designing methylation-specific PCR (MSP) or bisulfite
#' sequencing PCR (BSP) primers.
#'
#' This function does \strong{not} verify that a returned primer pair is
#' actually methylation-discriminating (i.e. that its binding site overlaps
#' a CpG position and would fail to anneal to the other converted state).
#' That is a property of exactly where in the template the primer lands,
#' which you should confirm yourself -- e.g. by checking the primer's
#' position against \code{\link{find_cpg_islands}} output, or by comparing
#' the primer sequence directly against \code{bisulfite_convert(sequence,
#' assume_methylated_cpg = !assume_methylated_cpg)} to see whether it still
#' matches. It only automates the conversion + primer-design steps.
#'
#' @param sequence A character string, the original (unconverted) genomic
#'   DNA sequence.
#' @param methylation_state Which reference state to design primers
#'   against: \code{"unmethylated"} (default) or \code{"methylated"}.
#' @param ... Additional arguments passed to \code{\link{design_primers}}
#'   (e.g. \code{product_size_range}, \code{num_return}).
#'
#' @return A data frame like \code{\link{design_primers}}'s output, with an
#'   added \code{methylation_state} column.
#'
#' @examples
#' \dontrun{
#' design_bisulfite_primers(promoter_seq,
#'   methylation_state = "unmethylated", product_size_range = c(100, 200)
#' )
#' }
#'
#' @export
design_bisulfite_primers <- function(sequence,
                                     methylation_state = c("unmethylated", "methylated"),
                                     ...) {
  methylation_state <- match.arg(methylation_state)
  converted <- bisulfite_convert(
    sequence,
    assume_methylated_cpg = (methylation_state == "methylated")
  )
  pairs <- design_primers(converted, ...)
  if (nrow(pairs) > 0) {
    pairs$methylation_state <- methylation_state
  }
  pairs
}
