# Methods

Detailed description of computational methods used in the rRNA NOR Pipeline.

---

## Data requirements

### Input format
All input sequences must be in FASTA format, aligned with MAFFT or equivalent
multiple sequence alignment tool. Gap characters (`-`) are preserved and used
as features by the gap pattern classifier.

### Sequence ID format
For the positional analysis (Step 9) to work, sequence IDs should encode
chromosomal start and end positions in the format:
```
ANY_PREFIX/START_END
```
For example: `NOR2_region_(15.8_Mb)/8825_17270`

If IDs do not follow this format, Step 9 is skipped automatically.

---

## Step-by-step methods

### Step 1 — Reference sequence download
Reference 45S rDNA sequences are downloaded from NCBI GenBank using the
Biopython Entrez API. Downloaded sequences are cached locally so repeated
runs do not re-download. The pipeline uses three reference species chosen
to represent different phylogenetic distances from the study organism.

### Step 2 — Alignment loading
FASTA files are parsed using BioPython SeqIO. Sequence length, gap fraction,
and chromosomal position (if encoded in the ID) are extracted for each record.

### Step 3 — Feature extraction

**k-mer frequency vectors**
Each sequence is converted to a normalised k-mer frequency vector of length
4^k (256 for k=4). Gap characters and ambiguous bases (N) are stripped before
counting. Each vector sums to 1.0.

```
sequence → strip gaps → count all 4-mers → normalise → 256-dimensional vector
```

**Gap pattern vectors**
Each aligned sequence is divided into non-overlapping 200 bp windows. Each
window is scored as 1 (contains at least one gap character) or 0 (no gaps).
Vectors are zero-padded to equal length across all sequences.

```
aligned sequence → 200bp windows → binary gap presence → gap vector
```

**Combined feature matrix**
The k-mer and gap vectors are concatenated horizontally to form the combined
feature matrix used for clustering.

### Step 4 — Reference similarity scoring
Cosine similarity is computed between each Micro-Tom sequence's k-mer vector
and each reference species' k-mer vector:

```
similarity(a, b) = dot(a, b) / (||a|| × ||b||)
```

Values range from 0 (completely different composition) to 1 (identical
composition). This metric is alignment-free and robust to length differences.

### Step 5 — K-Means clustering
K-Means clustering (scikit-learn, k=9, n_init=20, random_state=42) is applied
to the StandardScaler-normalised combined feature matrix.

Cross-alignment consistency is measured using Adjusted Rand Index (ARI):
- ARI = 1.0: identical cluster assignments
- ARI = 0.0: cluster assignments no better than random

### Step 6 — Gap pattern Random Forest classifier
A Random Forest classifier (300 trees, random_state=42) is trained on gap
vectors to predict insertion type. Each sequence is labelled by its dominant
gap window (the 200 bp window with the most gap content). Sequences with no
gaps are excluded from training.

5-fold stratified cross-validation reports macro-averaged F1 score.
Feature importance identifies which 200 bp windows (chromosomal positions)
best distinguish between insertion types.

### Step 7 — UMAP visualisation
UMAP (n_neighbors=15, min_dist=0.1, random_state=42) projects k-mer vectors
from all Micro-Tom sequences and reference species into a shared 2D embedding.
Sequences are colored by K-Means cluster assignment.

### Step 8 — Conservation heatmap
For each alignment and reference species pair, 500 bp windows are extracted
from each sequence, converted to k-mer vectors, and compared to the reference
k-mer vector by cosine similarity. Windows are binned into 40 positions along
the gene body and mean similarity is plotted as a heatmap.

### Step 9 — Positional analysis
Chromosomal start positions are extracted from sequence IDs and plotted
against cluster assignment. This tests whether specific haplotypes are
spatially organised within the NOR chromosomal region.

### Step 10 — Gap feature importance
Random Forest feature importances from Step 6 are plotted against approximate
chromosomal position (window index × 200 bp). Known functional region
boundaries (18S, 25S, ITS1, ITS2) are overlaid for interpretation.

### Step 11 — Haplotype similarity
For each haplotype representative sequence, cosine similarity to each
reference species is computed and reported in a CSV table.

### Step 12 — Summary report
A cross-alignment comparison table is generated reporting sequence count,
mean length, cluster count, mean similarity to each reference species, and
percentage of sequences with gap-containing windows.

---

## Insertion analysis methods

### Conservation profiling
Per-position conservation is computed as the fraction of sequences matching
the majority-vote consensus base at each position.

### BLAST search
The consensus sequence is submitted to NCBI blastn against the nt database
with a Viridiplantae organism filter (E-value threshold 0.001, top 50 hits).
Results are cached as XML to avoid repeated submissions.

### Transposable element analysis
- **Size classification**: compared against known TE size ranges
- **TIR detection**: terminal 5-25 bp compared to reverse complement of
  opposite terminal, similarity computed for each length
- **Composition**: GC and AT content compared to known TE family signatures
- **Motif scanning**: regex patterns for known plant TE structural features

### Regulatory motif scanning
The consensus and individual sequences are scanned for known plant rDNA
regulatory sequences using regex patterns including transcription termination
signals, processing signals, transcription factor binding sites, and
repeat-associated elements.

---

## Dependencies

| Package | Version | Use |
|---------|---------|-----|
| biopython | ≥1.79 | Sequence I/O, NCBI API, BLAST |
| scikit-learn | ≥1.0 | K-Means, Random Forest, UMAP preprocessing |
| umap-learn | ≥0.5 | UMAP dimensionality reduction |
| matplotlib | ≥3.5 | All figures |
| seaborn | ≥0.11 | Heatmaps |
| pandas | ≥1.3 | Data tables and CSV output |
| numpy | ≥1.21 | Numerical operations |
| scipy | ≥1.7 | Statistical tests |
