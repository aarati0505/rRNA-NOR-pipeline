# rRNA NOR Pipeline

A computational pipeline for analyzing 45S ribosomal RNA gene sequences from the Nucleolus Organizer Region (NOR) of plant genomes, with cross-species comparative analysis using machine learning.

Developed as part of a Study Oriented Project at BITS Pilani, Hyderabad Campus.

---

## What this pipeline does

Plant genomes contain hundreds of tandemly repeated 45S rRNA gene copies clustered in the NOR region of specific chromosomes. These copies are not identical — they vary in sequence, contain insertions of different sizes, and can be grouped into haplotypes. This pipeline provides tools to:

- Load and process multiple MAFFT-aligned FASTA files from rDNA NOR regions
- Extract k-mer frequency features and gap pattern features from aligned sequences
- Compute sequence similarity to reference species (downloads automatically from NCBI)
- Cluster sequences using K-Means and compare clustering consistency across alignments
- Classify insertion types using a Random Forest classifier
- Visualize sequence diversity using UMAP dimensionality reduction
- Generate conservation heatmaps comparing Micro-Tom to reference species
- Analyze the spatial organisation of haplotypes across the NOR chromosomal region
- Characterize specific insertion sequences: BLAST search, TE analysis, regulatory motif scanning

---

## Pipeline overview

```
Your FASTA files          Reference sequences (NCBI)
      │                           │
      └──────────┬────────────────┘
                 ↓
        Feature extraction
      k-mer vectors + gap vectors
                 ↓
      ┌──────────┴──────────┐
      ↓                     ↓
  K-Means               Random Forest
  clustering            gap classifier
      ↓                     ↓
  UMAP plot         Gap importance plot
      ↓
  Conservation heatmap
      ↓
  Positional NOR map
      ↓
  Summary report (CSV)
```

---

## Repository structure

```
rRNA-NOR-pipeline/
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
│
├── pipeline/
│   ├── rRNA_pipeline.py       # Main ML pipeline (Steps 1-12)
│   ├── analyze_insertion.py   # Insertion sequence analysis (BLAST, TE, regulatory)
│   └── haplotype_analysis.R   # R scripts for haplotyping (ape, pegas, geneHapR)
│
├── demo/
│   ├── generate_demo_data.py  # Generates synthetic FASTA files for testing
│   └── demo_sequences.fasta   # Pre-generated synthetic demo data
│
└── docs/
    ├── methods.md             # Full step-by-step methodology
    └── pipeline_overview.png  # Visual pipeline diagram
```

---

## Installation

**Requirements:** Python 3.8 or higher, Windows/Mac/Linux

**Step 1 — Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/rRNA-NOR-pipeline.git
cd rRNA-NOR-pipeline
```

**Step 2 — Install dependencies**
```bash
pip install -r requirements.txt
```

---

## Quick start with demo data

To verify the pipeline works on your system before using your own data:

```bash
# Generate synthetic demo sequences
python demo/generate_demo_data.py

# Run the full pipeline on demo data
python pipeline/rRNA_pipeline.py --demo
```

This will create a `pipeline_output/` folder with all figures and reports generated from synthetic sequences.

---

## Usage with your own data

**Step 1 — Prepare your files**

Place your MAFFT-aligned FASTA files in a single folder. The pipeline expects:
- One or more alignment files (e.g. `Alignment_1.fasta`, `Alignment_2.fasta`)
- Optionally: haplotype summary files

Sequence IDs should ideally encode chromosomal position in the format:
```
NOR_region/START_END
```
for example: `NOR2_region_(15.8_Mb)/8825_17270`

This allows the pipeline to generate the positional NOR map (Step 9).

**Step 2 — Edit the configuration**

Open `pipeline/rRNA_pipeline.py` and edit the two lines at the top:
```python
ALIGNMENT_DIR = r"C:\path\to\your\fasta\folder"   # your folder
YOUR_EMAIL    = "your_email@example.com"            # for NCBI API
```

**Step 3 — Run**
```bash
python pipeline/rRNA_pipeline.py
```

**Step 4 — Run insertion analysis (optional)**

If you have a FASTA file containing specific insertion sequences to characterize:
```bash
# Edit FASTA_141BP path in analyze_insertion.py first
python pipeline/analyze_insertion.py
```

**Step 5 — Run haplotype analysis in R (optional)**

Open `pipeline/haplotype_analysis.R` in RStudio. Three approaches are available:

```r
# Method 1 — standard ape/pegas haplotyping
hap1 <- method1_haplotyping("path/to/alignment.fasta")

# Method 2 — geneHapR with gap filtering
hap2 <- method2_haplotyping("path/to/alignment.fasta")

# Find optimal number of clusters first
optimal_k <- find_optimal_clusters("path/to/alignment.fasta")

# K-Means sequence clustering
km <- cluster_sequences("path/to/alignment.fasta", n_clusters = optimal_k)
```

---

## Output files

The pipeline generates a `pipeline_output/` folder containing:

| File | Description |
|------|-------------|
| `step7_umap.png` | UMAP of all sequences colored by cluster, reference species as stars |
| `step8_conservation_heatmap.png` | Conservation along gene body vs reference species |
| `step9_positional_map.png` | Cluster distribution across NOR chromosomal region |
| `step10_gap_importance.png` | Random Forest feature importance for insertion classification |
| `step11_haplotype_similarity.csv` | Per-haplotype similarity to each reference species |
| `step12_summary_report.csv` | Cross-alignment comparison table |

For insertion analysis:

| File | Description |
|------|-------------|
| `step1_conservation_profile.png` | Per-position conservation, GC content, dinucleotide heatmap |
| `step2_blast_results.png` | BLAST hits ranked by identity and significance |
| `step2_blast_summary.csv` | Full BLAST results table |
| `step3_te_analysis.png` | Transposable element signature analysis |
| `step4_regulatory_analysis.png` | Regulatory motif positions and prevalence |
| `step5_summary.txt` | Plain text summary with biological interpretation |

---

## Reference species used

The pipeline automatically downloads these sequences from NCBI:

| Species | Accession | Role |
|---------|-----------|------|
| *Arabidopsis thaliana* | X52322 | Well-annotated reference |
| *Solanum lycopersicum* SL4.0 | AY305797 | Closest relative to Micro-Tom |
| *Oryza sativa* | AY373817 | Distant outgroup |

---

## Methods summary

**Feature extraction** — Each sequence is converted to a normalised k-mer frequency vector (k=4, 256 features) capturing sequence composition, and a binary gap pattern vector (200 bp windows) capturing insertion presence/absence.

**Clustering** — K-Means clustering (k=9) is applied to the combined feature matrix. Adjusted Rand Index is computed between alignments to assess consistency.

**Gap classification** — A Random Forest classifier (300 trees) is trained on gap vectors to predict insertion type. Feature importance identifies which chromosomal positions best distinguish insertion variants.

**Reference similarity** — Cosine similarity between k-mer vectors measures sequence divergence from reference species without requiring a new alignment.

**Dimensionality reduction** — UMAP projects sequences into 2D space for visualization, embedding Micro-Tom and reference sequences together for direct comparison.

---

## Citation

If you use this pipeline in your work, please cite:

```
Dhamele, A. (2025). rRNA NOR Pipeline:
GitHub: https://github.com/YOUR_USERNAME/rRNA-NOR-pipeline
BITS Pilani, Hyderabad Campus.
```

---

## License

MIT License — see LICENSE file for details.

---

## Author

**Aarati Dhamele**  
B.E. Computer Science+ M.Sc Biological Sciences  
BITS Pilani, Hyderabad Campus  
ID: 2023B1A70641H  

---

## Acknowledgements

This pipeline was developed under the Study Oriented Project program at BITS Pilani Hyderabad. Reference sequences obtained from NCBI GenBank. Pipeline uses Biopython, scikit-learn, UMAP-learn, matplotlib, and seaborn.
