# **Introduction**
*Genomic analysis made elegant using R*

## **Overview**

BioSeqAligner is a R package designed for DNA sequence analysis. Whether you're a bioinformatician, 
researcher, or student, this package provides intuitive functions 
for motif detection, sequence alignment
visualization, reverse complement calculation, DNA-to-RNA transcription,
and GC content calculation. It simplifies essential function for manipulating and analyzing DNA sequences.

## **Key Features**


-   🔍 Motif Finder: Identify all occurrences of a motif in a DNA sequence.
-   🧬 Reverse Complement: Generate the reverse complement of a sequence.
-   🧪 DNA to RNA Transcription: Convert DNA sequences to RNA by replacing thymine (T) with uracil (U).
-   📊 GC Content Calculator: Calculate the percentage of guanine (G) and cytosine (C) in your DNA sequence using S3 objects.
-   📈 Dot Plot Visualization: Create dot plots to visualize sequence alignments.
-   🧮 Pairwise Alignment: True global (Needleman–Wunsch) and local (Smith–Waterman) alignment with affine gap penalties, for DNA, RNA, and protein.
-   🎯 Scoring Schemes: Built-in nucleotide and BLOSUM62 substitution matrices, or supply your own.
-   🧷 Multiple Sequence Alignment: Progressive MSA with consensus and per-column conservation scoring.
-   📐 Alignment Metrics: Percent identity, similarity, coverage, and gap statistics in one call.
-   🗂️ FASTA / FASTQ I/O: Read and write real sequence files.
-   🚀 Batch (mini-BLAST) Search: Rank one query against many references.
-   🌳 Phylogenetics: Distance matrix and neighbour-joining tree straight from an MSA.
-   🎨 Rich Visualization: Colour-coded pairwise and MSA plots (ggplot2).
-   🖥️ Interactive Shiny App: `launch_aligner()` ties it all together in the browser.

## **Installation**

#### **For Github**

You can install this package directly from github

-   devtools::install_github(“Code-krafter22/BioSeqAligner”)

-   library(BioSeqAligner)

## **Packages required**

-   install.packages(“ggplot2”)

-   install.packages(“stringi”)

-   install.packages(“stringr”)

## **Load required libraries**

-   library(ggplot2)

-   library(stringi)

-   library(stringr)

## **Using BioSeqAligner**

### 1. **find_motif()**:

This function identifies all starting positions of a specified motif
within a DNA sequence. It performs a linear scan and returns the
positions where the motif matches the sequence.


#### **Example**

-   find_motif(“GATCGATCGTAT”, “GAT”)

-   ***Returns [1] 1 5***

### 2. **generate_dot_plot()**:

This function creates a dot plot to visualize the alignment between two
DNA sequences.

#### **Library required**

-   *library(ggplot2)*

#### **Example**

-   generate_dot_plot(“GATCGATCGTAT”, “GATATCGTCATC”)

***Returns a graph plot for sequence alignment where X-axis
contains the Sequence 1 and Y-axis contains Sequence 2. The dark red dot
signifies the similarities in the sequence and the blue dot signifies
the dissimilarities in the sequence.***

### **Dot plot**

-  Overview of how the plot looks

![](README_files/figure-markdown_strict/dot-plot-1.png)

### 3. **reverse_complement()**:

This function computes the reverse complement of a given DNA sequence by
reversing it and substituting complementary bases.

#### **Library required**

-   *library(stringi)*

#### **Example**

-   reverse_complement(“GATCGATCGTAT”)

-   ***Returns [1] “ATACGATCGATC”***

### 4. **transcribe\_dna()**:

This function converts a DNA sequence into an RNA sequence by replacing
thymine (T) with uracil (U).

#### **Library required**

-   *library(stringi)*

#### **Example**

-    transcribe_dna(“GATCGATCGTAT”)

-    ***Returns [1] “GAUCGAUCGUAU”***

### 5. **gc_content()**:

Constructs an S3 object representing a DNA sequence with validation to
ensure it contains only valid DNA bases (A, C, G, T). Calculates the GC
content (percentage of guanine and cytosine) in a DNA sequence
represented by a DNASequence object.

#### **Library required**

-   *library(stringr)*

#### **Example**

-   gc_content.S3(DNASequence(“ATGC”))

-   ***Returns [1] 50***

## **Sequence Alignment Engine**

Beyond the visualization helpers above, BioSeqAligner now performs *true*
sequence alignment using dynamic programming.

### **align()** — pairwise global or local alignment

Global (Needleman–Wunsch) and local (Smith–Waterman) alignment, both with
affine gap penalties. Works for DNA, RNA, and protein via the chosen scoring
matrix. Returns an `alignment` object with `print()` and `plot()` methods.

```r
a <- align("ACGTGGATCGA", "ACGTGCATCGA", method = "global")
print(a)   # side-by-side view with a match line, score, and % identity
plot(a)    # colour-coded tile plot (matches / mismatches / gaps)

# Local alignment finds a conserved region inside a longer sequence
align("TTACGTGGATT", "ACGTGGA", method = "local")

# Protein alignment with BLOSUM62
align("MKVLA", "MKVLW", method = "global", submat = scoring_matrix("BLOSUM62"))
```

### **alignment_stats()** — quality metrics

```r
alignment_stats(a)
# score, length, matches, mismatches, gaps, identity_pct, gap_pct, query_coverage_pct
```

### **msa_align()** — multiple sequence alignment

```r
seqs <- c(s1 = "ACGTGGAA", s2 = "ACGTGCAA", s3 = "ACGTGGATA")
m <- msa_align(seqs)
print(m)                    # aligned rows + consensus
conservation_scores(m)      # per-column conservation (0–1)
plot(m)                     # colour-coded MSA
plot(m, color_by = "conservation")
```

### **batch_align()** — mini-BLAST search

```r
refs <- c(refA = "TTACGTGGATT", refB = "GGGGCCCC", refC = "AAACGTGGAAA")
batch_align("ACGTGGA", refs, method = "local")   # ranked hit table
```

### **FASTA / FASTQ I/O**

```r
seqs <- read_fasta("sequences.fasta")
write_fasta(seqs, "out.fasta")
reads <- read_fastq("reads.fastq")
```

### **Phylogenetics** (requires the `ape` package)

```r
nj_tree(m)   # neighbour-joining tree built from the MSA distance matrix
```

### **Interactive app** (requires the `shiny` package)

```r
launch_aligner()   # paste/upload sequences, align, and view results in the browser
```

### **Optional packages**

-   install.packages("ape")     # for nj_tree()
-   install.packages("shiny")   # for launch_aligner()

## **Error**:

-    **Invalid DNA Sequences**: Functions will stop and display an error if the input contains characters other than A, C, G, or T.

#### **Example**

-   transcribe_dna("ATBX")
-   ***Error: Input is not a valid DNA sequence***

### **License**:

-   This project is licensed under the MIT License. See the LICENSE file for details.

### **Contact**:

-   For questions or feedback, reach out via GitHub Issues.

## **Happy sequencing!** 
