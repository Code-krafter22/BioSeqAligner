#' @title Read Sequences from a FASTA File
#'
#' @description Parses a FASTA file (or a character vector of FASTA lines) into
#'   a named character vector, where names are the header lines (without the
#'   leading \code{>}) and values are the concatenated sequences.
#'
#' @param path Path to a \code{.fasta}/\code{.fa} file, or a character vector
#'   of already-read lines.
#'
#' @return A named character vector of sequences.
#'
#' @examples
#' tmp <- tempfile(fileext = ".fasta")
#' writeLines(c(">seqA", "ACGT", "GGCA", ">seqB", "TTTT"), tmp)
#' read_fasta(tmp)
#'
#' @export
read_fasta <- function(path) {
  lines <- if (length(path) == 1 && file.exists(path)) readLines(path) else path
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  if (length(lines) == 0) {
    return(character(0))
  }
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) stop("No FASTA header (>) found.")

  headers <- sub("^>", "", lines[header_idx])
  # Sequence blocks run from each header to the line before the next
  starts <- header_idx + 1
  ends <- c(header_idx[-1] - 1, length(lines))
  seqs <- mapply(function(s, e) {
    if (s > e) {
      return("")
    }
    toupper(paste(lines[s:e], collapse = ""))
  }, starts, ends)

  stats::setNames(seqs, headers)
}

#' @title Write Sequences to a FASTA File
#'
#' @description Writes a named character vector of sequences to a FASTA file,
#'   wrapping long sequences at a fixed line width.
#'
#' @param seqs A named character vector of sequences.
#' @param path Output file path.
#' @param width Line-wrap width for sequence lines (default 60).
#'
#' @return The output path, invisibly.
#'
#' @examples
#' s <- c(gene1 = "ACGTACGT", gene2 = "TTGGCCAA")
#' out <- tempfile(fileext = ".fasta")
#' write_fasta(s, out)
#'
#' @export
write_fasta <- function(seqs, path, width = 60) {
  if (is.null(names(seqs))) names(seqs) <- paste0("seq", seq_along(seqs))
  con <- file(path, open = "wt")
  on.exit(close(con))
  for (i in seq_along(seqs)) {
    writeLines(paste0(">", names(seqs)[i]), con)
    s <- seqs[i]
    starts <- seq(1, nchar(s), by = width)
    writeLines(substring(s, starts, starts + width - 1), con)
  }
  invisible(path)
}

#' @title Read Sequences from a FASTQ File
#'
#' @description Parses a FASTQ file into a data frame of read IDs, sequences,
#'   and quality strings. Handles the standard 4-lines-per-record layout.
#'
#' @param path Path to a \code{.fastq}/\code{.fq} file, or a character vector
#'   of lines.
#'
#' @return A data frame with columns \code{id}, \code{sequence}, and
#'   \code{quality}.
#'
#' @examples
#' tmp <- tempfile(fileext = ".fastq")
#' writeLines(c("@r1", "ACGT", "+", "IIII"), tmp)
#' read_fastq(tmp)
#'
#' @export
read_fastq <- function(path) {
  lines <- if (length(path) == 1 && file.exists(path)) readLines(path) else path
  lines <- lines[lines != ""]
  if (length(lines) %% 4 != 0) {
    stop("FASTQ file does not contain a multiple of 4 lines.")
  }
  idx <- seq(1, length(lines), by = 4)
  data.frame(
    id = sub("^@", "", lines[idx]),
    sequence = toupper(lines[idx + 1]),
    quality = lines[idx + 3],
    stringsAsFactors = FALSE
  )
}
