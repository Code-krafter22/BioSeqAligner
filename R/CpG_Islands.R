#' @title Detect CpG Islands in a DNA Sequence
#'
#' @description
#' Scans a DNA sequence with a sliding window and flags regions matching the
#' standard sequence-based definition of a CpG island (Gardiner-Garden &
#' Frommer 1987, refined by Takai & Jones 2002): GC content above a
#' threshold, an observed/expected CpG dinucleotide ratio above a threshold,
#' sustained over a minimum length. Adjacent/overlapping windows that pass
#' are merged into a single island.
#'
#' This is a purely sequence-based, structural computation: it identifies
#' CpG-dense regions in the DNA sequence you supply (typically a gene
#' promoter). It does \strong{not} tell you whether that CpG island is
#' actually methylated in any real sample — methylation is a chemical mark
#' layered on top of the DNA sequence that varies by tissue, individual, and
#' disease state, and cannot be determined from sequence alone. Measuring
#' actual methylation status requires real bisulfite-sequencing or
#' methylation-array data (e.g. via Bioconductor's \code{methylKit} or
#' \code{bsseq}), which is out of scope for this function. See
#' \code{\link{cpg_disease_lookup}} for a curated (not computed) list of
#' genes whose promoter CpG islands are well-established in the literature
#' as recurrently hyper/hypomethylated in specific diseases.
#'
#' @param sequence A character string, the DNA sequence to scan (typically a
#'   gene promoter region).
#' @param window Sliding window size in bp. Default 200 (Gardiner-Garden &
#'   Frommer); use 500 for the stricter Takai & Jones convention.
#' @param step Step size in bp between windows. Default 1 (exhaustive scan).
#'   Larger steps run faster on long sequences at the cost of coordinate
#'   precision.
#' @param gc_threshold Minimum GC percent (0-100) for a window to pass.
#'   Default 50.
#' @param oe_threshold Minimum observed/expected CpG ratio for a window to
#'   pass. Default 0.6. Expected CpG count in a window is
#'   \code{(count(C) * count(G)) / window_length}, the standard formula.
#' @param min_length Minimum length (bp) for a merged region to be reported
#'   as an island. Default equal to \code{window}.
#'
#' @return A data frame, one row per detected CpG island, with columns
#'   \code{start}, \code{end}, \code{length}, \code{gc_percent}, and
#'   \code{obs_exp_ratio} (values are the mean across the island's
#'   constituent windows). Empty (0-row) if no island is found.
#'
#' @examples
#' # A CpG-dense region flanked by AT-rich sequence, e.g. a synthetic promoter
#' cpg_core <- paste(rep("CG", 150), collapse = "")
#' at_flank <- paste(rep("AT", 100), collapse = "")
#' find_cpg_islands(paste0(at_flank, cpg_core, at_flank), step = 5)
#'
#' # In real use, scan an actual promoter sequence:
#' # promoter <- read_fasta("my_gene_promoter.fasta")
#' # find_cpg_islands(promoter)
#'
#' @export
find_cpg_islands <- function(sequence, window = 200, step = 1,
                             gc_threshold = 50, oe_threshold = 0.6,
                             min_length = window) {
  sequence <- toupper(as.character(sequence))
  n <- nchar(sequence)
  no_islands <- function() data.frame(
    start = integer(0), end = integer(0), length = integer(0),
    gc_percent = numeric(0), obs_exp_ratio = numeric(0),
    stringsAsFactors = FALSE
  )
  if (n < window) {
    return(no_islands())
  }

  chars <- strsplit(sequence, "")[[1]]
  is_c <- as.integer(chars == "C")
  is_g <- as.integer(chars == "G")
  # "CG" dinucleotide occurs at position i when chars[i]=="C" & chars[i+1]=="G"
  is_cg <- as.integer(chars[-n] == "C" & chars[-1] == "G")

  cum_c <- c(0, cumsum(is_c))
  cum_g <- c(0, cumsum(is_g))
  cum_cg <- c(0, cumsum(is_cg))

  starts <- seq(1, n - window + 1, by = step)
  ends <- starts + window - 1

  c_count <- cum_c[ends + 1] - cum_c[starts]
  g_count <- cum_g[ends + 1] - cum_g[starts]
  # CpG dinucleotides starting within [start, end-1]
  cg_count <- cum_cg[pmin(ends, n - 1) + 1] - cum_cg[starts]

  gc_pct <- 100 * (c_count + g_count) / window
  expected_cg <- (c_count * g_count) / window
  oe_ratio <- ifelse(expected_cg > 0, cg_count / expected_cg, 0)

  passes <- gc_pct >= gc_threshold & oe_ratio >= oe_threshold
  if (!any(passes)) {
    return(no_islands())
  }

  pass_starts <- starts[passes]
  pass_ends <- ends[passes]
  pass_gc <- gc_pct[passes]
  pass_oe <- oe_ratio[passes]

  # Merge overlapping/adjacent passing windows into islands
  order_idx <- order(pass_starts)
  pass_starts <- pass_starts[order_idx]
  pass_ends <- pass_ends[order_idx]
  pass_gc <- pass_gc[order_idx]
  pass_oe <- pass_oe[order_idx]

  island_start <- pass_starts[1]
  island_end <- pass_ends[1]
  island_gc <- c(pass_gc[1])
  island_oe <- c(pass_oe[1])
  islands <- list()

  flush <- function(s, e, gc, oe) {
    data.frame(
      start = s, end = e, length = e - s + 1,
      gc_percent = mean(gc), obs_exp_ratio = mean(oe),
      stringsAsFactors = FALSE
    )
  }

  for (i in seq_along(pass_starts)[-1]) {
    if (pass_starts[i] <= island_end + 1) {
      island_end <- max(island_end, pass_ends[i])
      island_gc <- c(island_gc, pass_gc[i])
      island_oe <- c(island_oe, pass_oe[i])
    } else {
      islands[[length(islands) + 1]] <- flush(island_start, island_end, island_gc, island_oe)
      island_start <- pass_starts[i]
      island_end <- pass_ends[i]
      island_gc <- c(pass_gc[i])
      island_oe <- c(pass_oe[i])
    }
  }
  islands[[length(islands) + 1]] <- flush(island_start, island_end, island_gc, island_oe)

  result <- do.call(rbind, islands)
  result <- result[result$length >= min_length, ]
  rownames(result) <- NULL
  result
}

#' @title Curated Disease Associations for Well-Known CpG Island Genes
#'
#' @description
#' Looks up a gene symbol in a small, hand-curated reference table of genes
#' whose promoter CpG islands are well-established in the published
#' literature as recurrently hyper- or hypomethylated in specific diseases
#' (e.g. \code{MLH1} promoter hypermethylation silencing DNA mismatch repair
#' in sporadic colorectal cancer).
#'
#' This is a \strong{lookup against a fixed, non-exhaustive list}, not a
#' computation on any sequence you provide — it cannot tell you whether a
#' given gene's CpG island is actually methylated in your sample, and
#' finding no entry for a gene does not mean no disease association exists
#' in the wider literature, only that it isn't in this small curated table.
#' Use \code{\link{find_cpg_islands}} separately to locate CpG islands from
#' sequence; the two are intentionally not combined into one function so
#' structural sequence facts are never conflated with literature curation.
#'
#' @param gene A gene symbol, e.g. \code{"MLH1"}. Matching is case-insensitive.
#'
#' @return A one-row data frame (\code{gene}, \code{disease},
#'   \code{mechanism}) if found, or a zero-row data frame with a message
#'   printed if the gene isn't in the table.
#'
#' @examples
#' cpg_disease_lookup("MLH1")
#' cpg_disease_lookup("MGMT")
#'
#' @export
cpg_disease_lookup <- function(gene) {
  ref <- .cpg_disease_reference()
  hit <- ref[toupper(ref$gene) == toupper(gene), ]
  if (nrow(hit) == 0) {
    message(
      "'", gene, "' is not in this package's small curated reference table. ",
      "This does not mean no disease association exists in the literature -- ",
      "only that it isn't in this list. See ?cpg_disease_lookup for scope."
    )
  }
  rownames(hit) <- NULL
  hit
}

# A small, hand-curated table of well-established (textbook-level, widely
# published) gene-disease CpG island methylation associations. Not
# exhaustive, not computed -- see cpg_disease_lookup() documentation.
.cpg_disease_reference <- function() {
  data.frame(
    gene = c(
      "MLH1", "BRCA1", "MGMT", "CDKN2A", "VHL", "RASSF1A",
      "GSTP1", "APC", "FMR1", "SNRPN"
    ),
    disease = c(
      "Sporadic colorectal cancer (Lynch-like)",
      "Breast and ovarian cancer",
      "Glioblastoma (predicts temozolomide response)",
      "Multiple cancers (melanoma, lung, pancreatic, others)",
      "Clear cell renal cell carcinoma",
      "Lung, breast, and other cancers",
      "Prostate cancer",
      "Colorectal cancer",
      "Fragile X syndrome",
      "Prader-Willi syndrome"
    ),
    mechanism = c(
      "Promoter hypermethylation silences DNA mismatch repair",
      "Promoter hypermethylation silences DNA damage repair",
      "Promoter hypermethylation silences DNA repair enzyme",
      "Promoter hypermethylation silences cell-cycle regulator (p16)",
      "Promoter hypermethylation silences tumor suppressor",
      "Promoter hypermethylation silences Ras-association tumor suppressor",
      "Promoter hypermethylation silences detoxification enzyme (most common epigenetic alteration in prostate cancer)",
      "Promoter hypermethylation silences tumor suppressor (Wnt pathway)",
      "CGG-repeat expansion triggers CpG island hypermethylation, silencing the gene",
      "Loss of normal paternal-allele methylation pattern at the imprinting center"
    ),
    stringsAsFactors = FALSE
  )
}
