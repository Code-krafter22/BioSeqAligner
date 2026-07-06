# **Introduction**
*Genomic analysis made elegant using R*

## **Overview**

BioSeqAligner is an R package for biological sequence analysis and alignment.
Whether you're a bioinformatician, researcher, or student, this package
provides intuitive functions for sequence manipulation (motif detection,
reverse complement, DNA-to-RNA transcription, GC content) alongside a full
alignment engine (global and local pairwise alignment, multiple sequence
alignment, alignment metrics, FASTA/FASTQ I/O, batch search, phylogenetics,
and rich visualization) — so the common tasks that would otherwise require
juggling several bioinformatics tools live in one package.

## **Key Features**

-   🔍 Motif Finder: Identify all occurrences of a motif in a DNA sequence.
-   🧬 Reverse Complement: Generate the reverse complement of a sequence.
-   🧪 DNA to RNA Transcription: Convert DNA sequences to RNA by replacing thymine (T) with uracil (U).
-   📊 GC Content Calculator: Calculate the percentage of guanine (G) and cytosine (C) in your DNA sequence using S3 objects.
-   📈 Dot Plot Visualization: Create dot plots to visualize sequence alignments.
-   🧮 Pairwise Alignment: True global (Needleman-Wunsch) and local (Smith-Waterman) alignment with affine gap penalties, for DNA, RNA, and protein.
-   🎯 Scoring Schemes: Built-in nucleotide and BLOSUM62 substitution matrices, or supply your own.
-   🧷 Multiple Sequence Alignment: Progressive MSA with consensus and per-column conservation scoring.
-   📐 Alignment Metrics: Percent identity, similarity, coverage, and gap statistics in one call.
-   🗂️ FASTA / FASTQ I/O: Read and write real sequence files.
-   🚀 Batch (mini-BLAST) Search: Rank one query against many references.
-   🌳 Phylogenetics: Distance matrix and neighbour-joining tree straight from an MSA.
-   🎨 Rich Visualization: Colour-coded pairwise and MSA plots (ggplot2).
-   🖥️ Interactive Shiny App: `launch_aligner()` ties it all together in the browser.

## **Installation**

#### **From GitHub**

You can install this package directly from GitHub.

```r
devtools::install_github("Code-krafter22/BioSeqAligner")
library(BioSeqAligner)
```

## **Packages Required**

```r
install.packages("ggplot2")
install.packages("stringi")
install.packages("stringr")
```

Two features are optional and only need their package if you use them:

```r
install.packages("ape")     # for nj_tree()
install.packages("shiny")   # for launch_aligner()
```

## **Load Required Libraries**

```r
library(ggplot2)
library(stringi)
library(stringr)
```

## **Using BioSeqAligner**

### 1. **find_motif()**

Identifies all starting positions of a specified motif within a DNA sequence.
Performs a linear scan and returns the positions where the motif matches the
sequence.

```r
find_motif("GATCGATCGTAT", "GAT")
# Returns [1] 1 5
```

### 2. **generate_dot_plot()**

Creates a dot plot to visualize the alignment between two DNA sequences.
Requires `ggplot2`.

```r
generate_dot_plot("GATCGATCGTAT", "GATATCGTCATC")
# Returns a plot where the X-axis is Sequence 1 and the Y-axis is Sequence 2.
# Dark red dots mark matching positions; blue dots mark mismatches.
```

Overview of how the plot looks:

![](README_files/figure-markdown_strict/dot-plot-1.png)

### 3. **reverse_complement()**

Computes the reverse complement of a given DNA sequence by reversing it and
substituting complementary bases. Requires `stringi`.

```r
reverse_complement("GATCGATCGTAT")
# Returns [1] "ATACGATCGATC"
```

### 4. **transcribe_dna()**

Converts a DNA sequence into an RNA sequence by replacing thymine (T) with
uracil (U). Requires `stringi`.

```r
transcribe_dna("GATCGATCGTAT")
# Returns [1] "GAUCGAUCGUAU"
```

### 5. **gc_content.S3()**

Constructs an S3 object representing a DNA sequence with validation to ensure
it contains only valid DNA bases (A, C, G, T), then calculates the GC content
(percentage of guanine and cytosine). Requires `stringr`.

```r
gc_content.S3(DNASequence("ATGC"))
# Returns [1] 50
```

### 6. **align()**

Performs pairwise sequence alignment via dynamic programming: global
(Needleman-Wunsch) or local (Smith-Waterman), both with affine gap penalties.
Works for DNA, RNA, and protein depending on the scoring matrix supplied.
Returns an `alignment` object with `print()` and `plot()` methods.

```r
a <- align("ACGTGGATCGA", "ACGTGCATCGA", method = "global")
print(a)   # side-by-side view with a match line, score, and % identity
plot(a)    # colour-coded tile plot (matches / mismatches / gaps)

# Local alignment finds a conserved region inside a longer sequence
align("TTACGTGGATT", "ACGTGGA", method = "local")

# Protein alignment with BLOSUM62
align("MKVLA", "MKVLW", method = "global", submat = scoring_matrix("BLOSUM62"))
```

### 7. **scoring_matrix()**

Builds the substitution matrix used to score aligned residues: a nucleotide
match/mismatch scheme, or the BLOSUM62 matrix for proteins.

```r
scoring_matrix("nucleotide", match = 2, mismatch = -1)
scoring_matrix("BLOSUM62")
```

### 8. **alignment_stats()**

Returns a one-row summary of alignment-quality metrics: length, matches,
mismatches, gaps, percent identity, percent gaps, and query coverage.

```r
alignment_stats(a)
# score, length, matches, mismatches, gaps, identity_pct, gap_pct, query_coverage_pct
```

### 9. **msa_align()**

Aligns three or more sequences using a progressive multiple sequence
alignment strategy, returning an `msa` object with the aligned rows and a
consensus sequence. `conservation_scores()` gives a per-column conservation
score, and `plot()` draws a colour-coded MSA.

```r
seqs <- c(s1 = "ACGTGGAA", s2 = "ACGTGCAA", s3 = "ACGTGGATA")
m <- msa_align(seqs)
print(m)                    # aligned rows + consensus
conservation_scores(m)      # per-column conservation (0-1)
plot(m)                     # colour-coded MSA
plot(m, color_by = "conservation")
```

### 10. **batch_align()**

Aligns a single query sequence against a set of reference sequences and
returns a ranked results table — a lightweight, BLAST-like search over a
small, in-memory database.

```r
refs <- c(refA = "TTACGTGGATT", refB = "GGGGCCCC", refC = "AAACGTGGAAA")
batch_align("ACGTGGA", refs, method = "local")   # ranked hit table
```

### 11. **read_fasta() / write_fasta() / read_fastq()**

Reads and writes real sequence files instead of pasted strings.

```r
seqs <- read_fasta("sequences.fasta")
write_fasta(seqs, "out.fasta")
reads <- read_fastq("reads.fastq")
```

### 12. **nj_tree()**

Computes a distance matrix from an MSA and builds a neighbour-joining
phylogenetic tree, bridging alignment straight into phylogenetics. Requires
`ape`.

```r
nj_tree(m)   # neighbour-joining tree built from the MSA distance matrix
```

### 13. **launch_aligner()**

Launches an interactive Shiny app that ties the package together: paste or
upload sequences, run pairwise or multiple sequence alignment, and view
colour-coded plots, summary metrics, and a neighbour-joining tree in the
browser. Requires `shiny`.

```r
launch_aligner()   # paste/upload sequences, align, and view results in the browser
```

## **Error Handling**

**Invalid DNA Sequences**: Functions will stop and display an error if the
input contains characters other than A, C, G, or T (unless the function
explicitly supports RNA, protein, or IUPAC ambiguity codes, e.g.
`validate_sequence()`).

```r
transcribe_dna("ATBX")
# Error: Input is not a valid DNA sequence.
```

## **License**

This project is licensed under the MIT License. See the LICENSE file for
details.

## **Contact**

For questions or feedback, reach out via GitHub Issues.

## **Happy sequencing!**
