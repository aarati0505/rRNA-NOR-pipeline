# Methodology

Complete step-by-step methodology for extracting and analyzing 45S rRNA genes
from published whole genome sequences.

---

## Overview

The pipeline extracts individual 45S rRNA gene copies from a whole chromosome
sequence, aligns them, and analyzes sequence diversity. The organism used was
*Solanum lycopersicum* cultivar Micro-Tom, where the NOR region spans approximately
15 Mb on the short arm of chromosome 2, containing ~1500 rRNA gene copies.

The 45S rRNA gene unit contains the following annotated regions (5' to 3'):
```
5'ETS — 18S — ITS1 — 5.8S — ITS2 — 25S — 3'ETS — IGS (promoter)
```

**Note:** All sequences obtained from NCBI were in reverse complement orientation
and were handled accordingly throughout the analysis.

---

## Step 1 — Download chromosome sequence

Download the whole chromosome 2 sequence of *Solanum lycopersicum* Micro-Tom
from NCBI in FASTA format. Build a local BLAST database:

```bash
makeblastdb \
  -in [input_fasta_file] \
  -dbtype nucl \
  -out [output_database_name]
```

**Parameters:**
- `-dbtype nucl` — nucleotide database (use `prot` for protein)
- `-in` — path to the downloaded chromosome FASTA file
- `-out` — name for the output database files

---

## Step 2 — BLAST known sequences against the genome

Use the already-characterized 25S, 18S, and rDNA gene promoter sequences as
queries to locate every gene copy in the chromosome. These three sequences are
highly conserved across copies and serve as anchors for gene boundary detection.

```bash
blastn \
  -query [query_fasta_file] \
  -db [database_name] \
  -perc_identity 95 \
  -qcov_hsp_perc 90 \
  -outfmt 6 \
  -out [output_file.txt]
```

**Parameters:**
- `-perc_identity 95` — only report hits with >=95% sequence identity
- `-qcov_hsp_perc 90` — query must cover >=90% of the alignment
- `-outfmt 6` — tabular output format

**Output columns (outfmt 6):**

| Column | Description |
|--------|-------------|
| 1 | Query accession.version |
| 2 | Subject accession.version |
| 3 | Percentage of identical matches |
| 4 | Alignment length |
| 5 | Number of mismatches |
| 6 | Number of gap openings |
| 7 | Start of alignment in query |
| 8 | End of alignment in query |
| 9 | Start of alignment in subject |
| 10 | End of alignment in subject |
| 11 | Expect value |
| 12 | Bit score |

---

## Step 3 — Convert BLAST output to Excel

Convert the tab-separated BLAST output (.txt) to Excel format. The numbers in
columns 9 and 10 (start/end of alignment in subject) give the chromosomal
positions of each 25S, 18S, and rDNA promoter hit — these are used to define
gene boundaries.

---

## Step 4 — Identify gene patterns using Excel formulas

Two splicing patterns were identified based on the arrangement of query hits:

**Pattern 1** (complete gene with 18S):
```
3' --[25S region]----[18S region]----[rDNA promoter]-- 5'
     ^                                               ^
  Start of alignment                          End of alignment
  in subject (col 9)                       in subject (col 10)
```

**Pattern 2** (gene without 18S in alignment):
```
3' --[25S region]--------------------[rDNA promoter]-- 5'
     ^                                               ^
  Start of alignment                          End of alignment
  in subject (col 9)                       in subject (col 10)
```

Excel formulas were used to:
- Label each row as "Pattern Found" or "No Pattern"
- Extract the start and end positions for each complete gene
- Output two columns: gene start position and gene end position

---

## Step 5 — Extract gene sequences using Biopython

Use the start/end positions from Step 4 to slice individual gene sequences from
the chromosome FASTA file:

```python
from Bio import SeqIO

def slice_fasta(input_fasta, output_fasta, ranges,
                header_format="{id}/{start}_{end}"):
    """
    Slices a FASTA sequence into sub-sequences based on given coordinate ranges.

    Parameters:
        input_fasta  : Path to the input chromosome FASTA file
        output_fasta : Path to the output FASTA file for sliced genes
        ranges       : List of (start, end) tuples from Excel output
        header_format: Format string for output sequence headers
    """
    with open(input_fasta, 'r') as infile, \
         open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            seq = record.seq
            for start, end in ranges:
                # Biopython uses 0-based indexing -- subtract 1 from start
                slice_seq    = seq[start-1:end]
                header       = header_format.format(
                    id=record.id, start=start, end=end
                )
                slice_record = SeqIO.SeqRecord(
                    slice_seq, id=header, description=""
                )
                SeqIO.write(slice_record, outfile, "fasta")

# Example usage:
input_fasta  = "tomato_chr2.fasta"
output_fasta = "sliced_genes.fasta"
ranges       = [(244, 7779), (9725, 17270)]  # from Excel output
slice_fasta(input_fasta, output_fasta, ranges)
```

This produces a FASTA file containing all ~1521 unaligned gene sequences of
varying lengths due to insertions present at different positions in different copies.

---

## Step 6 — Multiple sequence alignment with MAFFT

Align all extracted gene sequences using MAFFT (Multiple Alignment using Fast
Fourier Transform):

```bash
mafft --auto --thread -1 sliced_genes.fasta > aligned_genes.fasta
```

- `--auto` — automatically selects the best alignment strategy
- `--thread -1` — use all available CPU cores

After alignment, sequences were sorted in Jalview software by pairwise identity
to arrange them in chromosomal order (based on start position encoded in IDs).

**Note:** The initial alignment did not include the 3'ETS region.
A second alignment was produced by extending each gene 900 bp upstream of
the 25S start position to capture the 3'ETS region (~790 bp).

---

## Step 7 — Alignment curation and annotation

Three alignment versions were produced:

| Alignment | Sequences | Length | Description |
|-----------|-----------|--------|-------------|
| Alignment 1 | 1521 | 56,734 bp | All sequences including large insertions |
| Alignment 2 | 1506 | 27,018 bp | With 18S insertion variants retained |
| Alignment 3 | 1491 | 19,287 bp | Large-insertion sequences removed |

**Annotation in Jalview:**
- Known regions (25S, 18S, 3'ETS, rDNA promoter) were annotated by BLASTing
  reference sequences against the alignment
- All gaps after the rDNA gene promoter region were manually annotated as
  GAP1, GAP2, GAP3... in order from the promoter
- An Excel table was created recording for each gap:
  - Gap name
  - Length of gap (bp)
  - Number of genes without the gap filled
  - Number of genes with the gap fully filled
  - Number of genes with the gap partially filled

Alignment 3 was produced by identifying and removing genes with large filled
gaps from Alignment 1, then re-running MAFFT. This reduced alignment length
from 56,734 bp to 19,287 bp.

---

## Step 8 — Haplotype analysis in R

Two approaches were used to identify haplotypes from the aligned sequences.

### Method 1 — Standard R libraries (ape + pegas)

```r
library(ape)
library(pegas)

# Load aligned FASTA file
fasta_file    <- file.choose()   # interactive file selection
dna_sequences <- read.dna(fasta_file, format = "fasta")

# Identify haplotypes
haplotypes <- haplotype(dna_sequences)

# Summarize results
hap_summary <- haploNet(haplotypes)
print(haplotypes)
summary(haplotypes)

# Optional: plot haplotype network
plot(hap_summary,
     size        = attr(hap_summary, "freq"),
     scale.ratio = 0.5)

# Save haplotype sequences to FASTA
write.dna(haplotypes,
          file   = "haplotypes_output.fasta",
          format = "fasta")
```

### Method 2 — geneHapR package

```r
library(geneHapR)
library(Biostrings)
library(rtracklayer)

# Load FASTA file
fasta_file <- "C:/path/to/Alignment_3.fasta"
dna_seqs   <- readDNAStringSet(fasta_file)

# Assign names if missing
if (is.null(names(dna_seqs)) ||
    any(is.na(names(dna_seqs))) ||
    any(names(dna_seqs) == "")) {
    names(dna_seqs) <- paste0("seq", seq_along(dna_seqs))
}

# Convert to alignment matrix
seq_matrix <- as.matrix(dna_seqs)

# Remove columns with too many gaps (allow up to 80% gaps per column)
max_gap_freq        <- 0.8
gap_freq            <- colSums(seq_matrix == "-") / nrow(seq_matrix)
keep_cols           <- which(gap_freq <= max_gap_freq)
seq_matrix_filtered <- seq_matrix[, keep_cols]

# Convert filtered matrix back to DNAStringSet
filtered_seqs       <- DNAStringSet(apply(
    seq_matrix_filtered, 1, paste0, collapse = ""
))
names(filtered_seqs) <- names(dna_seqs)

# Identify haplotypes
hap_result <- seqs2hap(filtered_seqs, max.gap.freq = 2)

# View and save
print(hap_result$hap)
print(hap_result$hapInfo)
write.csv(hap_result$hap,
          "haplotypes.csv",      row.names = FALSE)
write.csv(hap_result$hapInfo,
          "hap_assignments.csv", row.names = FALSE)
```

---

## Step 9 —  insertion characterization

**Extraction method:**
1. Identify all gap positions annotated as 141 bp in the Excel gap table
2. Extract the filled sequences from those positions using Biopython
3. Create a separate FASTA file with the extracted sequences
4. Align with MAFFT and visualize in Jalview

**Analysis performed on extracted sequences:**
- Per-position conservation scoring across all copies
- GC content and dinucleotide frequency profiling
- BLAST against NCBI nt database (Viridiplantae organism filter)
- Transposable element signature analysis:
  - Terminal Inverted Repeat (TIR) detection
  - Size classification against known TE families
  - AT/GC composition comparison to known TE signatures
- Regulatory motif scanning:
  - Transcription termination signals (T-rich, Sal box)
  - rRNA processing signals
  - Transcription factor binding sites (TATA box, GC box, CCAAT)

---

## Software and tools used

| Tool | Version | Purpose |
|------|---------|---------|
| NCBI BLAST+ | 2.13+ | Local BLAST database and search |
| Biopython | 1.79+ | Sequence slicing, FASTA I/O, NCBI API |
| MAFFT | 7.x | Multiple sequence alignment |
| Jalview | 2.11+ | Alignment visualization and annotation |
| Geneious | 11+ | Alignment visualization |
| R / RStudio | 4.x | Haplotype analysis |
| ape | CRAN | Phylogenetic and haplotype analysis |
| pegas | CRAN | Haplotype network analysis |
| geneHapR | CRAN | Gene haplotype statistics |
| DnaSP | 6.x | Population genetics and haplotyping |
| scikit-learn | 1.0+ | K-Means, Random Forest, ML analysis |
| umap-learn | 0.5+ | UMAP dimensionality reduction |
| matplotlib | 3.5+ | Figure generation |
| seaborn | 0.11+ | Heatmap visualization |

---

## Data availability

The genome sequence of *Solanum lycopersicum* cultivar Micro-Tom used in this
study is available from NCBI (chromosome 2 accession: AP028936).

Analyzed sequences and results are available from the authors upon reasonable
request pending publication.
