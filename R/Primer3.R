# Locate the primer3_core executable and its thermodynamic parameter files.
# Both are required for primer3_core to run at all; neither ships inside this
# R package, so we look for a system installation (e.g. `brew install primer3`
# on macOS, `apt install primer3` on Linux) and fail with install instructions
# if it isn't found.
.primer3_exe <- function(primer3_path = NULL) {
  if (!is.null(primer3_path)) {
    if (!file.exists(primer3_path)) {
      stop("primer3_path does not exist: ", primer3_path)
    }
    return(primer3_path)
  }
  exe <- Sys.which("primer3_core")
  if (nzchar(exe)) {
    return(unname(exe))
  }
  stop(
    "Could not find 'primer3_core' on your PATH. Install Primer3 first:\n",
    "  macOS:  brew install primer3\n",
    "  Linux:  apt-get install primer3 (or conda install -c bioconda primer3)\n",
    "Then either make sure primer3_core is on your PATH, or pass its full ",
    "path via the primer3_path argument."
  )
}

# Locate the primer3_config/ directory of thermodynamic parameter files,
# required by PRIMER_THERMODYNAMIC_PARAMETERS_PATH. We search a few
# conventional install locations relative to the primer3_core binary.
.primer3_config_dir <- function(exe, config_path = NULL) {
  if (!is.null(config_path)) {
    return(config_path)
  }
  candidates <- c(
    file.path(dirname(exe), "primer3_config"),
    file.path(dirname(dirname(exe)), "share", "primer3", "primer3_config"),
    file.path(dirname(dirname(exe)), "share", "primer3_config")
  )
  found <- candidates[dir.exists(candidates)]
  if (length(found) == 0) {
    stop(
      "Could not locate the primer3_config/ directory (thermodynamic ",
      "parameter files) near '", exe, "'. Pass its path explicitly via ",
      "the config_path argument."
    )
  }
  paste0(found[1], "/")
}

#' @title Design PCR Primers with Primer3
#'
#' @description
#' Runs the Primer3 engine (\url{https://primer3.org}) on a target sequence to
#' design forward/reverse primer pairs, using Primer3's nearest-neighbour
#' thermodynamic model for melting temperature and secondary-structure
#' checks. Requires a system installation of Primer3 (\code{primer3_core} on
#' the \code{PATH}); this function does not reimplement Primer3, it drives
#' the real engine and parses its output into a data frame.
#'
#' @param sequence A character string, the DNA template to design primers
#'   from (e.g. a gene/transcript sequence, or the region around a locus of
#'   interest).
#' @param target_start,target_length Optional 1-based start position and
#'   length of a region primers must flank (e.g. an exon-exon junction or a
#'   SNP). If either is \code{NULL}, primers are picked anywhere in
#'   \code{sequence}.
#' @param product_size_range A length-2 integer vector, the acceptable PCR
#'   product size range.
#' @param num_return Maximum number of primer pairs to return.
#' @param primer3_path Optional full path to the \code{primer3_core}
#'   executable, if it is not on your \code{PATH}.
#' @param config_path Optional full path to Primer3's \code{primer3_config/}
#'   directory, if it cannot be auto-detected next to the executable.
#'
#' @return A data frame, one row per primer pair, with columns
#'   \code{pair_id}, \code{left_seq}, \code{right_seq}, \code{left_tm},
#'   \code{right_tm}, \code{left_gc_percent}, \code{right_gc_percent},
#'   \code{product_size}, and \code{penalty} (lower is better).
#'
#' @examples
#' \dontrun{
#' design_primers(
#'   "ATGGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGC",
#'   product_size_range = c(70, 120)
#' )
#' }
#'
#' @export
design_primers <- function(sequence,
                           target_start = NULL, target_length = NULL,
                           product_size_range = c(100, 300),
                           num_return = 5,
                           primer3_path = NULL, config_path = NULL) {
  sequence <- toupper(as.character(sequence))
  exe <- .primer3_exe(primer3_path)
  cfg <- .primer3_config_dir(exe, config_path)

  lines <- c(
    "SEQUENCE_ID=biosaligner",
    paste0("SEQUENCE_TEMPLATE=", sequence),
    "PRIMER_TASK=generic",
    "PRIMER_PICK_LEFT_PRIMER=1",
    "PRIMER_PICK_RIGHT_PRIMER=1",
    paste0("PRIMER_NUM_RETURN=", as.integer(num_return)),
    paste0("PRIMER_PRODUCT_SIZE_RANGE=", product_size_range[1], "-", product_size_range[2]),
    paste0("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=", cfg)
  )
  if (!is.null(target_start) && !is.null(target_length)) {
    lines <- c(lines, paste0("SEQUENCE_TARGET=", target_start, ",", target_length))
  }
  lines <- c(lines, "=")

  input_file <- tempfile()
  writeLines(lines, input_file)
  on.exit(unlink(input_file))

  out <- suppressWarnings(system2(exe, args = character(0), stdin = input_file, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")

  kv <- strsplit(out[grepl("^[A-Z0-9_]+=", out)], "=", fixed = FALSE)
  kv <- lapply(kv, function(x) c(x[1], paste(x[-1], collapse = "=")))
  vals <- stats::setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))

  if (!is.null(status) && status != 0 || "PRIMER_ERROR" %in% names(vals)) {
    msg <- if ("PRIMER_ERROR" %in% names(vals)) vals[["PRIMER_ERROR"]] else paste(out, collapse = "\n")
    stop("Primer3 reported an error: ", msg)
  }

  n_pairs <- as.integer(vals[["PRIMER_PAIR_NUM_RETURNED"]])
  if (is.na(n_pairs) || n_pairs == 0) {
    warning("Primer3 did not find any primer pairs for this sequence/constraints.")
    return(data.frame(
      pair_id = integer(0), left_seq = character(0), right_seq = character(0),
      left_tm = numeric(0), right_tm = numeric(0),
      left_gc_percent = numeric(0), right_gc_percent = numeric(0),
      product_size = integer(0), penalty = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(0:(n_pairs - 1), function(i) {
    g <- function(key) unname(vals[[sprintf(key, i)]])
    data.frame(
      pair_id = i,
      left_seq = g("PRIMER_LEFT_%d_SEQUENCE"),
      right_seq = g("PRIMER_RIGHT_%d_SEQUENCE"),
      left_tm = as.numeric(g("PRIMER_LEFT_%d_TM")),
      right_tm = as.numeric(g("PRIMER_RIGHT_%d_TM")),
      left_gc_percent = as.numeric(g("PRIMER_LEFT_%d_GC_PERCENT")),
      right_gc_percent = as.numeric(g("PRIMER_RIGHT_%d_GC_PERCENT")),
      product_size = as.integer(g("PRIMER_PAIR_%d_PRODUCT_SIZE")),
      penalty = as.numeric(g("PRIMER_PAIR_%d_PENALTY")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @title Design Primers and Check Their Specificity in One Step
#'
#' @description
#' Runs \code{\link{design_primers}} to generate candidate primer pairs, then
#' immediately checks each primer's specificity with \code{\link{batch_align}}
#' against a set of sequences you supply (e.g. the intended target plus known
#' paralogs/isoforms/off-targets). Ties primer design and specificity
#' screening together in a single call, instead of designing primers in one
#' tool and separately checking them in another.
#'
#' @param sequence The DNA template to design primers from (passed to
#'   \code{\link{design_primers}}).
#' @param off_targets A named character vector of sequences to screen each
#'   candidate primer against — typically the intended target plus anything
#'   it might cross-react with (paralogs, other isoforms, related species).
#' @param identity_threshold Percent identity (0-100) above which a hit
#'   against something other than the intended target is flagged as a
#'   specificity concern. Default 80.
#' @param coverage_threshold Minimum percent of the primer's length that must
#'   actually be aligned (query coverage) for a hit to count as a
#'   specificity concern. Default 70. This guards against local alignment
#'   reporting a misleadingly high identity on a short coincidental match
#'   (e.g. 100% identity over just 3 of 20 primer bases) — a real off-target
#'   binding site needs most of the primer to align, not a few lucky bases.
#' @param intended_target Optional name (matching one entry in
#'   \code{off_targets}) of the sequence the primers are meant to bind. Hits
#'   against this name are not counted against specificity.
#' @param ... Additional arguments passed to \code{\link{design_primers}}
#'   (e.g. \code{product_size_range}, \code{num_return}, \code{target_start},
#'   \code{target_length}).
#'
#' @return A data frame like \code{\link{design_primers}}'s output, with two
#'   extra columns: \code{max_offtarget_identity} (the highest percent
#'   identity either primer scored against anything other than
#'   \code{intended_target}, considering only hits that meet
#'   \code{coverage_threshold}) and \code{specific} (\code{TRUE} if that
#'   value is below \code{identity_threshold} for both primers).
#'
#' @examples
#' \dontrun{
#' targets <- c(
#'   GENE_A = "ATGGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGTAGCTAGC",
#'   PARALOG_B = "ATGGCTAGCATCGATCGTTTTTAGCATCGATCGTTTTTAGCATCGATCGTAGCTAGC"
#' )
#' design_and_validate_primers(targets["GENE_A"], targets,
#'   intended_target = "GENE_A", product_size_range = c(70, 120)
#' )
#' }
#'
#' @export
design_and_validate_primers <- function(sequence, off_targets,
                                        identity_threshold = 80,
                                        coverage_threshold = 70,
                                        intended_target = NULL, ...) {
  pairs <- design_primers(sequence, ...)
  if (nrow(pairs) == 0) {
    return(pairs)
  }

  screen_targets <- off_targets
  if (!is.null(intended_target) && intended_target %in% names(off_targets)) {
    screen_targets <- off_targets[names(off_targets) != intended_target]
  }

  max_identity <- function(primer_seq) {
    if (length(screen_targets) == 0) {
      return(0)
    }
    hits <- batch_align(primer_seq, screen_targets, method = "local")
    hits <- hits[hits$query_coverage_pct >= coverage_threshold, ]
    if (nrow(hits) == 0) {
      return(0)
    }
    max(hits$identity_pct)
  }

  pairs$max_offtarget_identity <- mapply(function(l, r) {
    max(max_identity(l), max_identity(r))
  }, pairs$left_seq, pairs$right_seq)
  pairs$specific <- pairs$max_offtarget_identity < identity_threshold

  pairs[order(pairs$penalty), ]
}
