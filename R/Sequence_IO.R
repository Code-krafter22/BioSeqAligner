#' @title Parse FASTA-Formatted Text
#'
#' @description Parses FASTA-formatted content already in memory (e.g. pasted
#'   text, or lines already read some other way) into a named character
#'   vector, where names are the header lines (without the leading
#'   \code{>}) and values are the concatenated sequences. Unlike
#'   \code{\link{read_fasta}}, this never touches the filesystem — pass it
#'   raw text directly instead of letting it guess whether a string is a file
#'   path or literal content.
#'
#' @param text A single character string of FASTA content (with embedded
#'   newlines), or a character vector already split into lines.
#'
#' @return A named character vector of sequences.
#'
#' @examples
#' parse_fasta(">seqA\nACGT\nGGCA\n>seqB\nTTTT")
#' parse_fasta(c(">seqA", "ACGT", "GGCA", ">seqB", "TTTT"))
#'
#' @export
parse_fasta <- function(text) {
  lines <- unlist(strsplit(text, "\n", fixed = TRUE))
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

#' @title Read Sequences from a FASTA File
#'
#' @description Reads a FASTA file from disk and parses it into a named
#'   character vector, where names are the header lines (without the leading
#'   \code{>}) and values are the concatenated sequences. \code{path} is
#'   always treated as a filesystem path — to parse FASTA content you
#'   already have in memory (e.g. user-pasted text), use
#'   \code{\link{parse_fasta}} instead, which never touches the filesystem.
#'
#' @param path Path to a \code{.fasta}/\code{.fa} file.
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
  if (length(path) != 1 || !file.exists(path)) {
    stop(
      "'path' must be the path to an existing file. To parse FASTA text you ",
      "already have in memory, use parse_fasta() instead."
    )
  }
  parse_fasta(readLines(path))
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

#' @title Parse FASTQ-Formatted Text
#'
#' @description Parses FASTQ-formatted content already in memory into a data
#'   frame of read IDs, sequences, and quality strings. Handles the standard
#'   4-lines-per-record layout. Unlike \code{\link{read_fastq}}, this never
#'   touches the filesystem.
#'
#' @param text A single character string of FASTQ content (with embedded
#'   newlines), or a character vector already split into lines.
#'
#' @return A data frame with columns \code{id}, \code{sequence}, and
#'   \code{quality}.
#'
#' @examples
#' parse_fastq("@r1\nACGT\n+\nIIII")
#'
#' @export
parse_fastq <- function(text) {
  lines <- unlist(strsplit(text, "\n", fixed = TRUE))
  lines <- lines[lines != ""]
  if (length(lines) %% 4 != 0) {
    stop("FASTQ content does not contain a multiple of 4 lines.")
  }
  idx <- seq(1, length(lines), by = 4)
  data.frame(
    id = sub("^@", "", lines[idx]),
    sequence = toupper(lines[idx + 1]),
    quality = lines[idx + 3],
    stringsAsFactors = FALSE
  )
}

#' @title Read Sequences from a FASTQ File
#'
#' @description Reads a FASTQ file from disk and parses it into a data frame
#'   of read IDs, sequences, and quality strings. \code{path} is always
#'   treated as a filesystem path — to parse FASTQ content you already have
#'   in memory, use \code{\link{parse_fastq}} instead, which never touches
#'   the filesystem.
#'
#' @param path Path to a \code{.fastq}/\code{.fq} file.
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
  if (length(path) != 1 || !file.exists(path)) {
    stop(
      "'path' must be the path to an existing file. To parse FASTQ text you ",
      "already have in memory, use parse_fastq() instead."
    )
  }
  parse_fastq(readLines(path))
}
