## R Markdown

## **Introduction**

BioSeqAligner is a R package designed for DNA sequence analysis. It
includes functions for motif detection, sequence alignment
visualization, reverse complement calculation, DNA-to-RNA transcription,
and GC content calculation.

## **Installation**

For local or GitHub installations:
devtools::install\_github(“Code-krafter22/BioSeqAligner”)

library(BioSeqAligner)

## **Packages required**

-   install.packages(“ggplot2”)

-   install.packages(“stringi”)

-   install.packages(“stringr”)

## **Load required libraries**

-   library(ggplot2)

-   library(stringi)

-   library(stringr)

## **Using BioSeqAligner**

### 1. **find\_motif() **:

This function identifies all starting positions of a specified motif
within a DNA sequence. It performs a linear scan and returns the
positions where the motif matches the sequence.

#### A. **Method 1: Define sequence and motif and then call the function**

seq &lt;- (“ATCGATCGATGC”), motif &lt;- (“ATC”)

find\_motif(seq,motif)

-   ***returns [1] 1 5***

#### B. **Method 2: Specify the sequence and motif in the function**

1.  find\_motif(“GATCGATCGTAT”, “GAT”)

-   ***Returns [1] 1 6***

1.  find\_motif(“AAAAAA”, “TT”)

-   ***Returns an empty vector***

### 2. **generate\_dot\_plot()**:

This function creates a dot plot to visualize the alignment between two
DNA sequences.

#### **Library required**

-   *library(ggplot2)*

#### A. **Method 1: Define sequence and then call the function**

-   seq1 &lt;- (“ATCGATCGATGC”), Seq2 &lt;- (“ATCGGCTACGTA”)

-   generate\_dot\_plot(seq1,seq2)

#### B. **Method 2: Specify the sequences in the function**

-   generate\_dot\_plot(“GATCGATCGTAT”, “GATATCGTCATC”)

***Both method returns a graph plot for sequence alignment where X-axis
contains the Sequence 1 and Y-axis contains Sequence 2. The dark red dot
signifies the similarities in the sequence and the blue dot signifies
the dissimilarities in the sequence.***

### **Dot plot**

![](README_files/figure-markdown_strict/dot-plot-1.png)

### 3. **reverse\_complement()**:

This function computes the reverse complement of a given DNA sequence by
reversing it and substituting complementary bases.

#### **Library required**

-   *library(stringi)*

#### A. **Method 1: Define sequence and then call the function**

-   seq &lt;- (“ATCGATCGATGC”)
-   reverse\_complement(seq)
-   ***returns \[1] “GCATCGATCGAT” ***

#### B. **Method 2: Specify the sequences in the function**

-   reverse\_complement(“GATCGATCGTAT”)
-   ***returns \[1] “ATACGATCGATC” ***

### 4. **transcribe\_dna()**:

This function converts a DNA sequence into an RNA sequence by replacing
thymine (T) with uracil (U).

#### **Library required**

-   *library(stringi)*

#### A. **Method 1: Define sequence and then call the function**

-   seq &lt;- (“ATCGATCGATGC”)
-   transcribe\_dna(seq)
-   ***returns \[1] “AUCGAUCGAUGC” ***

#### B. **Method 2: Specify the sequences in the function**

-   transcribe\_dna(“GATCGATCGTAT”)
-   ***returns \[1] “GAUCGAUCGUAU” ***

### 5. **gc\_content()**:

Constructs an S3 object representing a DNA sequence with validation to
ensure it contains only valid DNA bases (A, C, G, T). Calculates the GC
content (percentage of guanine and cytosine) in a DNA sequence
represented by a DNASequence object.

#### **Library required**

-   *library(stringr)*

#### A. **Method 1: Define sequence and then call the function**

-   dna\_seq &lt;- DNASequence(“ATGCGC”)
-   gc\_content.S3(dna\_seq)
-   ***returns \[1] 66.66667 ***

#### B. **Method 2: Specify the sequences in the function**

-   gc\_content.S3(DNASequence(“ATGC”))
-   ***returns \[1] 50 ***
