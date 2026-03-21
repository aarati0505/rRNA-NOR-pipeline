# =============================================================================
# 45S rRNA NOR Region — Full Comparative ML Pipeline
# For: Micro-Tom (3 alignments) vs Arabidopsis, SL4.0, Oryza sativa
#
# HOW TO RUN:
#   1. Put this file in the same folder as your FASTA files
#   2. Change ALIGNMENT_DIR below to your folder path
#   3. Change YOUR_EMAIL to your real email address
#   4. Press the Run button (triangle) in VS Code
#
# OUTPUT:
#   A folder called "pipeline_output" will be created with all results
# =============================================================================

# ── Two lines you MUST edit before running ────────────────────────────────────

ALIGNMENT_DIR = r"C:\Users\YourName\rRNA_project"  # ← change this to your folder path
YOUR_EMAIL    = "your_email@example.com"            # ← change this to your email

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ─────────────────────────────────────────────────────────────────────────────
print("Loading libraries... please wait")

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # works without a display — saves files instead of showing popups
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from itertools import product
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

# Check all libraries are installed before going further
missing = []
try:
    from Bio import Entrez, SeqIO
except ImportError:
    missing.append("biopython")
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    from sklearn.preprocessing import LabelEncoder, StandardScaler
    from sklearn.metrics import adjusted_rand_score
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
except ImportError:
    missing.append("scikit-learn")
try:
    import umap
except ImportError:
    missing.append("umap-learn")

if missing:
    print(f"\nERROR: Missing libraries: {', '.join(missing)}")
    print("Please run in the VS Code terminal:")
    print(f"  pip install {' '.join(missing)}")
    sys.exit(1)

print("All libraries loaded successfully.\n")

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
Entrez.email = YOUR_EMAIL
OUTPUT_DIR   = os.path.join(ALIGNMENT_DIR, "pipeline_output")
KMER_SIZE    = 4
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Your exact file names — edit if yours are named differently
ALIGNMENT_FILES = {
    "MSA1_all_gaps":   "Alignment 1(with all the gaps).fasta",
    "MSA2_18S_gap":    "Alignment 2(with 18s gap).fasta",
    "MSA3_small_gaps": "Alignment 3(with only small gaps).fasta",
}
HAPLOTYPE_FILES = {
    "HAP1": "mafft_aligned_Summary_Filtered_Haplotype_file_without_removing_1.fasta",
    "HAP2": "mafft_aligned_Summary_Filtered_Haplotype_file_without_removing_2.fasta",
    "HAP3": "mafft_aligned_Summary_Filtered_Haplotype_file_without_removing_3.fasta",
}

# Reference sequences from NCBI
REFERENCE_ACCESSIONS = {
    "Arabidopsis_thaliana":       "X52322",
    "Solanum_lycopersicum_SL4.0": "AY305797",
    "Oryza_sativa":               "AY373817",
}

# Colors for each reference species in plots
SPECIES_COLORS = {
    "Arabidopsis_thaliana":       "#E63946",
    "Solanum_lycopersicum_SL4.0": "#FF9800",
    "Oryza_sativa":               "#2196F3",
}

# ─────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

# All possible k-mers (256 for k=4)
ALL_KMERS = [''.join(p) for p in product('ATGC', repeat=KMER_SIZE)]

def kmer_vector(seq, k=KMER_SIZE):
    """
    Convert a DNA sequence to a normalised k-mer frequency vector.
    Gaps (-) and ambiguous bases (N) are removed first.
    Example: "ATCGATCG" → vector of 256 numbers summing to 1.0
    """
    seq = seq.upper().replace('-', '').replace('N', '')
    if len(seq) < k:
        return np.zeros(len(ALL_KMERS))
    counts = Counter(
        seq[i:i+k] for i in range(len(seq) - k + 1)
        if all(b in 'ATGC' for b in seq[i:i+k])
    )
    vec = np.array([counts.get(km, 0) for km in ALL_KMERS], dtype=float)
    total = vec.sum()
    return vec / total if total > 0 else vec

def gap_vector(seq, window=200):
    """
    Divide sequence into windows and mark each as 1 (has gap) or 0 (no gap).
    This turns your insertion data into ML features directly.
    """
    seq   = seq.upper()
    wins  = [seq[i:i+window] for i in range(0, len(seq), window)]
    return np.array([1.0 if '-' in w else 0.0 for w in wins])

def cosine_similarity(a, b):
    """How similar are two vectors? 1.0 = identical, 0.0 = completely different."""
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom > 0 else 0.0

def extract_position(seq_id):
    """
    Get chromosomal start position from IDs like:
    NOR2_region_(15.8_Mb)/8825_17270  →  8825
    """
    try:
        return int(seq_id.split('/')[-1].split('_')[0])
    except Exception:
        return 0

def load_fasta(path, label=None):
    """Load a FASTA file and return a list of sequence records."""
    rows = []
    for rec in SeqIO.parse(path, "fasta"):
        rows.append({
            "id":       rec.id,
            "sequence": str(rec.seq),
            "length":   len(rec.seq),
            "position": extract_position(rec.id),
            "source":   label or os.path.basename(path),
        })
    return rows

def save_figure(filename):
    """Save the current matplotlib figure to the output folder."""
    path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {path}")

def section(title):
    """Print a section header to the terminal."""
    print("\n" + "="*60)
    print(f"  {title}")
    print("="*60)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Download reference sequences from NCBI
# ─────────────────────────────────────────────────────────────────────────────

def step1_download_references():
    section("STEP 1 — Downloading reference sequences from NCBI")
    refs = {}

    for species, accession in REFERENCE_ACCESSIONS.items():
        out_path = os.path.join(OUTPUT_DIR, f"{species}.fasta")

        # Use cached file if already downloaded
        if os.path.exists(out_path):
            rec = next(SeqIO.parse(out_path, "fasta"))
            print(f"  (already downloaded) {species}")
        else:
            print(f"  Downloading {species} ({accession}) from NCBI...")
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=accession,
                    rettype="fasta", retmode="text"
                )
                rec = SeqIO.read(handle, "fasta")
                handle.close()
                SeqIO.write(rec, out_path, "fasta")
                time.sleep(0.4)   # be polite to NCBI servers
            except Exception as e:
                print(f"  WARNING: Could not download {species}: {e}")
                print(f"  Skipping {species} — check your internet connection")
                continue

        refs[species] = {
            "id":       species,
            "sequence": str(rec.seq),
            "length":   len(rec.seq),
        }
        print(f"  {species}: {refs[species]['length']:,} bp")

    if not refs:
        print("\n  WARNING: No reference sequences downloaded.")
        print("  Steps 4 and 8 (similarity analysis) will be skipped.")

    return refs

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Load your Micro-Tom alignment files
# ─────────────────────────────────────────────────────────────────────────────

def step2_load_alignments():
    section("STEP 2 — Loading your Micro-Tom alignment files")

    all_data = {}
    all_hap  = {}

    # Load the three main alignment files
    for label, fname in ALIGNMENT_FILES.items():
        path = os.path.join(ALIGNMENT_DIR, fname)
        if not os.path.exists(path):
            print(f"  WARNING: Could not find '{fname}'")
            print(f"  Make sure the file is in: {ALIGNMENT_DIR}")
            continue
        rows = load_fasta(path, label)
        df   = pd.DataFrame(rows)
        all_data[label] = df
        print(f"  {label}: {len(df):,} sequences | "
              f"avg length {df['length'].mean():,.0f} bp | "
              f"positions found: {(df['position']>0).sum()}")

    # Load the six haplotype files
    for label, fname in HAPLOTYPE_FILES.items():
        path = os.path.join(ALIGNMENT_DIR, fname)
        if not os.path.exists(path):
            print(f"  WARNING: Could not find '{fname}' — skipping")
            continue
        rows = load_fasta(path, label)
        df   = pd.DataFrame(rows)
        all_hap[label] = df
        print(f"  {label}: {len(df)} haplotype sequences")

    if not all_data:
        print("\nERROR: No alignment files were loaded.")
        print(f"Please check that ALIGNMENT_DIR is set correctly: {ALIGNMENT_DIR}")
        sys.exit(1)

    return all_data, all_hap

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — Feature extraction
# ─────────────────────────────────────────────────────────────────────────────

def step3_features(all_data):
    section(f"STEP 3 — Feature extraction (k={KMER_SIZE}, {len(ALL_KMERS)} k-mer features)")

    features = {}
    for label, df in all_data.items():
        print(f"  Processing {label}...")

        # K-mer vectors — biological sequence content
        print(f"    Computing k-mer vectors for {len(df)} sequences...")
        X_kmer = np.array([kmer_vector(s) for s in df["sequence"]])

        # Gap vectors — insertion fingerprint
        print(f"    Computing gap pattern vectors...")
        gap_vecs = [gap_vector(s) for s in df["sequence"]]
        max_len  = max(len(v) for v in gap_vecs)
        X_gap    = np.array([
            np.pad(v, (0, max_len - len(v))) for v in gap_vecs
        ])

        # Combined matrix
        X_combo = np.hstack([X_kmer, X_gap])

        features[label] = {
            "X_kmer":  X_kmer,
            "X_gap":   X_gap,
            "X_combo": X_combo,
            "df":      df.copy(),
        }
        print(f"    k-mer matrix: {X_kmer.shape} | "
              f"gap matrix: {X_gap.shape} | "
              f"combined: {X_combo.shape}")

    return features

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Similarity to reference species
# ─────────────────────────────────────────────────────────────────────────────

def step4_reference_similarity(features, refs):
    section("STEP 4 — Similarity to reference species")

    if not refs:
        print("  Skipping — no reference sequences available")
        return features

    # Pre-compute reference k-mer vectors
    ref_vecs = {sp: kmer_vector(d["sequence"]) for sp, d in refs.items()}

    for label, feat in features.items():
        df     = feat["df"].copy()
        X_kmer = feat["X_kmer"]

        for species, rv in ref_vecs.items():
            col = "sim_" + species
            sims = [cosine_similarity(row, rv) for row in X_kmer]
            df[col] = sims
            print(f"  {label} vs {species[:25]:25s}: "
                  f"mean={np.mean(sims):.4f}  "
                  f"max={np.max(sims):.4f}  "
                  f"min={np.min(sims):.4f}")

        features[label]["df"] = df

    return features

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — K-Means clustering + cross-MSA consistency
# ─────────────────────────────────────────────────────────────────────────────

def step5_clustering(features):
    section("STEP 5 — K-Means clustering (k=9 clusters)")

    # Use 9 clusters to match your known haplotype count
    N_CLUSTERS    = 9
    cluster_store = {}

    for label, feat in features.items():
        print(f"  Clustering {label}...")
        X = StandardScaler().fit_transform(feat["X_combo"])
        km = KMeans(n_clusters=N_CLUSTERS, random_state=42,
                    n_init=20, max_iter=500)
        cl = km.fit_predict(X)
        features[label]["df"]["cluster"] = cl
        cluster_store[label] = cl

        sizes = dict(sorted(Counter(cl).items()))
        print(f"    Cluster sizes: {sizes}")

    # Cross-MSA agreement using Adjusted Rand Index
    # ARI = 1.0 means the two MSAs produce identical groupings
    # ARI = 0.0 means the groupings are random with respect to each other
    labels_list = list(cluster_store.items())
    print("\n  Cross-MSA agreement (Adjusted Rand Index):")
    print("  (1.0 = perfect agreement, 0.0 = no agreement)")
    for i in range(len(labels_list)):
        for j in range(i+1, len(labels_list)):
            la, a = labels_list[i]
            lb, b = labels_list[j]
            n   = min(len(a), len(b))
            ari = adjusted_rand_score(a[:n], b[:n])
            print(f"    {la} vs {lb}: ARI = {ari:.4f}")

    return features, cluster_store

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — Gap pattern classifier (Random Forest)
# ─────────────────────────────────────────────────────────────────────────────

def step6_gap_classifier(features):
    section("STEP 6 — Gap pattern classifier (Random Forest)")
    print("  This identifies which chromosomal positions best distinguish")
    print("  insertion types from each other.\n")

    gap_results = {}
    for label, feat in features.items():
        X_gap = feat["X_gap"]
        df    = feat["df"]

        # Label each sequence by which 200bp window has the most gaps
        # Sequences with no gaps get label -1 (clean / no insertion)
        row_gap_sum      = X_gap.sum(axis=1)
        dominant_window  = np.argmax(X_gap, axis=1)
        dominant_window[row_gap_sum == 0] = -1

        features[label]["df"]["dominant_gap_window"] = dominant_window

        # Only train on sequences that actually have gaps
        mask = dominant_window >= 0
        X_tr = X_gap[mask]
        y_tr = dominant_window[mask]

        n_seqs_with_gaps = mask.sum()
        print(f"  {label}:")
        print(f"    Sequences with gaps    : {n_seqs_with_gaps} / {len(mask)}")
        print(f"    Sequences without gaps : {(~mask).sum()} / {len(mask)}")

        if len(np.unique(y_tr)) < 2:
            print(f"    Not enough gap variation to train classifier — skipping\n")
            continue

        le   = LabelEncoder()
        y_e  = le.fit_transform(y_tr)

        rf = RandomForestClassifier(
            n_estimators=300,
            random_state=42,
            n_jobs=-1
        )

        # 5-fold cross validation
        cv = cross_val_score(
            rf, X_tr, y_e,
            cv=StratifiedKFold(n_splits=min(5, len(np.unique(y_e)))),
            scoring='f1_macro'
        )
        rf.fit(X_tr, y_e)

        importances  = rf.feature_importances_
        top_idx      = np.argsort(importances)[::-1][:10]
        top_bp_pos   = [i * 200 for i in top_idx]

        print(f"    CV F1 score            : {cv.mean():.3f} ± {cv.std():.3f}")
        print(f"    Top 10 diagnostic bp positions: {top_bp_pos}\n")

        gap_results[label] = {
            "rf":           rf,
            "importances":  importances,
            "top_positions":top_bp_pos,
            "cv_f1":        cv.mean(),
        }

    return features, gap_results

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7 — UMAP visualisation
# ─────────────────────────────────────────────────────────────────────────────

def step7_umap(features, refs):
    section("STEP 7 — UMAP visualisation")

    ref_vecs   = np.array([kmer_vector(d["sequence"]) for d in refs.values()])
    ref_names  = list(refs.keys())
    ref_colors = list(SPECIES_COLORS.values())

    n_msa = len(features)
    fig, axes = plt.subplots(1, n_msa, figsize=(8 * n_msa, 7))
    if n_msa == 1:
        axes = [axes]

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)

    for ax, (label, feat) in zip(axes, features.items()):
        X_kmer   = feat["X_kmer"]
        df       = feat["df"]
        clusters = df.get("cluster", pd.Series(np.zeros(len(df)))).values

        # Embed Micro-Tom + reference sequences in the same 2D space
        X_all = np.vstack([X_kmer, ref_vecs]) if len(ref_vecs) > 0 else X_kmer
        emb   = reducer.fit_transform(X_all)
        mt_e  = emb[:len(X_kmer)]
        rf_e  = emb[len(X_kmer):]

        sc = ax.scatter(
            mt_e[:, 0], mt_e[:, 1],
            c=clusters, cmap='tab10',
            alpha=0.5, s=15, zorder=2
        )

        # Reference species as large stars
        for i, (sp, c) in enumerate(zip(ref_names, ref_colors)):
            if i < len(rf_e):
                ax.scatter(
                    rf_e[i, 0], rf_e[i, 1],
                    c=c, marker='*', s=400,
                    edgecolors='black', linewidths=0.8,
                    zorder=10, label=sp.replace('_', ' ')
                )

        ax.set_title(f"{label}\n(colored by cluster)", fontsize=11,
                     fontweight='bold')
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
        ax.legend(fontsize=7, loc='lower right')
        plt.colorbar(sc, ax=ax, label="Cluster")

    fig.suptitle(
        "UMAP: Micro-Tom 45S rRNA Sequences — All Three Alignments\n"
        "Stars = Reference species  |  Color = K-Means cluster",
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout()
    save_figure("step7_umap.png")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8 — Conservation heatmap
# ─────────────────────────────────────────────────────────────────────────────

def step8_conservation_heatmap(features, refs):
    section("STEP 8 — Conservation heatmap vs reference species")

    if not refs:
        print("  Skipping — no reference sequences available")
        return

    WINDOW   = 500    # compare sequences in 500 bp windows
    N_SAMPLE = 200    # use 200 sequences per MSA for speed
    N_BINS   = 40     # divide gene into 40 position bins

    ref_vecs = {sp: kmer_vector(d["sequence"]) for sp, d in refs.items()}
    rows     = []

    for label, feat in features.items():
        df = feat["df"]
        print(f"  Processing {label} ({min(N_SAMPLE, len(df))} sequences)...")

        for sp, rv in ref_vecs.items():
            all_sims = []
            for seq in df["sequence"].iloc[:N_SAMPLE]:
                seq_clean = seq.replace('-', '').upper()
                for i in range(0, len(seq_clean) - WINDOW, WINDOW):
                    wv  = kmer_vector(seq_clean[i:i+WINDOW])
                    all_sims.append((i + WINDOW // 2, cosine_similarity(wv, rv)))

            if not all_sims:
                continue

            pos_df = pd.DataFrame(all_sims, columns=["pos", "sim"])
            binned = pos_df.groupby(
                pd.cut(pos_df["pos"], bins=N_BINS)
            )["sim"].mean()

            row = {"MSA": label, "Species": sp.replace('_', ' ')}
            for k, v in enumerate(binned.values):
                row[f"w{k}"] = v
            rows.append(row)

    if not rows:
        print("  No data generated for heatmap")
        return

    hm_df = pd.DataFrame(rows).set_index(["MSA", "Species"])
    hm_df = hm_df.apply(pd.to_numeric, errors='coerce').fillna(0)

    fig, ax = plt.subplots(figsize=(22, max(4, len(hm_df) * 0.6)))
    sns.heatmap(
        hm_df, cmap="YlOrRd", ax=ax,
        vmin=0, vmax=1,
        xticklabels=5,
        cbar_kws={"label": "k-mer cosine similarity to reference"}
    )
    ax.set_title(
        "Conservation of Micro-Tom 45S rRNA Along Gene Body\n"
        "vs Arabidopsis, Solanum SL4.0, and Oryza sativa\n"
        "(Low values = regions unique to Micro-Tom)",
        fontsize=12, fontweight='bold'
    )
    ax.set_xlabel("Position along gene (5' → 3', binned)")
    ax.set_ylabel("")
    plt.tight_layout()
    save_figure("step8_conservation_heatmap.png")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 9 — Positional analysis in NOR
# ─────────────────────────────────────────────────────────────────────────────

def step9_positional_analysis(features):
    section("STEP 9 — Positional analysis of clusters in NOR2")
    print("  Testing whether haplotypes are spatially organised in the NOR.")
    print("  If yes → connects to Arabidopsis silencing-by-position literature.\n")

    n_msa = len(features)
    fig, axes = plt.subplots(n_msa, 1,
                              figsize=(16, 5 * n_msa),
                              sharex=False)
    if n_msa == 1:
        axes = [axes]

    cmap = plt.get_cmap('tab10')

    for ax, (label, feat) in zip(axes, features.items()):
        df = feat["df"].copy()

        if "cluster" not in df.columns:
            ax.set_title(f"{label} — no cluster data")
            continue
        if df["position"].sum() == 0:
            ax.set_title(f"{label} — no position data in sequence IDs")
            continue

        for cl in sorted(df["cluster"].unique()):
            sub = df[df["cluster"] == cl]
            ax.scatter(
                sub["position"], [cl] * len(sub),
                color=cmap(cl % 10),
                alpha=0.4, s=12,
                label=f"Cluster {cl} (n={len(sub)})"
            )

        # Annotate known region sizes if positions are available
        pos_range = df["position"].max() - df["position"].min()
        ax.set_title(
            f"{label} — Cluster distribution across NOR2\n"
            f"(Position range: {df['position'].min():,} – {df['position'].max():,} bp)",
            fontweight='bold'
        )
        ax.set_xlabel("Chromosomal position in NOR2 (bp)")
        ax.set_ylabel("Cluster ID")
        ax.set_yticks(sorted(df["cluster"].unique()))
        ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc='upper left')

    fig.suptitle(
        "Spatial Organisation of Haplotypes Across the NOR2 Region\n"
        "Each dot = one rRNA gene copy | Y-axis = assigned cluster",
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout()
    save_figure("step9_positional_map.png")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 10 — Gap importance plot
# ─────────────────────────────────────────────────────────────────────────────

def step10_gap_importance(gap_results):
    section("STEP 10 — Gap feature importance plots")

    if not gap_results:
        print("  No gap classifier results to plot")
        return

    n = len(gap_results)
    fig, axes = plt.subplots(1, n, figsize=(9 * n, 5))
    if n == 1:
        axes = [axes]

    # Approximate positions of known regions (in bp along alignment)
    # Edit these if you know the exact positions in your alignment
    KNOWN_REGIONS = {
        "25S":  (4500,  7000,  "#2196F3"),
        "18S":  (11000, 12800, "#E63946"),
        "ITS1": (7000,  7600,  "#4CAF50"),
        "ITS2": (10200, 11000, "#FF9800"),
    }

    for ax, (label, res) in zip(axes, gap_results.items()):
        imp = res["importances"]
        pos = np.arange(len(imp)) * 200   # window index → bp position

        ax.fill_between(pos, imp, alpha=0.75, color='#E63946',
                        label='Gap importance')
        ax.plot(pos, imp, color='#B71C1C', linewidth=0.8)

        # Shade known functional regions
        for rname, (start, end, color) in KNOWN_REGIONS.items():
            ax.axvspan(start, end, alpha=0.12, color=color, label=rname)

        ax.set_title(
            f"{label}\nWhich positions identify insertion type?\n"
            f"(CV F1 = {res['cv_f1']:.3f})",
            fontweight='bold'
        )
        ax.set_xlabel("Approximate bp position in alignment")
        ax.set_ylabel("Feature importance")
        ax.legend(fontsize=8)

    fig.suptitle(
        "Random Forest Gap Feature Importance\n"
        "Peaks outside 18S/25S = novel Micro-Tom-specific variable regions",
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout()
    save_figure("step10_gap_importance.png")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 11 — Haplotype file analysis
# ─────────────────────────────────────────────────────────────────────────────

def step11_haplotype_analysis(all_hap, refs):
    section("STEP 11 — Haplotype file analysis")

    if not all_hap:
        print("  No haplotype files loaded — skipping")
        return

    ref_vecs = {sp: kmer_vector(d["sequence"]) for sp, d in refs.items()} if refs else {}

    all_rows = []
    for label, df in all_hap.items():
        print(f"\n  {label} ({len(df)} haplotypes):")
        for _, row in df.iterrows():
            seq  = row["sequence"]
            vec  = kmer_vector(seq)
            info = {
                "haplotype_file": label,
                "id":             row["id"],
                "length":         row["length"],
                "gap_fraction":   seq.count('-') / max(len(seq), 1),
            }
            for sp, rv in ref_vecs.items():
                info[f"sim_{sp}"] = cosine_similarity(vec, rv)
            all_rows.append(info)

            # Print per-haplotype similarity
            sims = [f"{sp[:12]}={cosine_similarity(vec, rv):.4f}"
                    for sp, rv in ref_vecs.items()]
            print(f"    {row['id'][:50]:50s} | " + "  ".join(sims))

    hap_df = pd.DataFrame(all_rows)
    path   = os.path.join(OUTPUT_DIR, "step11_haplotype_similarity.csv")
    hap_df.to_csv(path, index=False)
    print(f"\n  Full results saved → {path}")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 12 — Final summary report
# ─────────────────────────────────────────────────────────────────────────────

def step12_summary_report(features, refs):
    section("STEP 12 — Summary report")

    rows = []
    for label, feat in features.items():
        df  = feat["df"]
        row = {
            "MSA":           label,
            "N_sequences":   len(df),
            "Mean_length_bp":f"{df['length'].mean():.0f}",
            "N_clusters":    df["cluster"].nunique() if "cluster" in df else "—",
        }

        # Similarity columns
        for sp in refs.keys():
            col = f"sim_{sp}"
            if col in df.columns:
                row[f"mean_sim_{sp[:15]}"] = f"{df[col].mean():.4f}"

        # Gap statistics
        if "dominant_gap_window" in df.columns:
            n_gaps = (df["dominant_gap_window"] >= 0).sum()
            row["pct_with_gaps"] = f"{100*n_gaps/len(df):.1f}%"

        rows.append(row)

    report = pd.DataFrame(rows)
    print(report.to_string(index=False))

    path = os.path.join(OUTPUT_DIR, "step12_summary_report.csv")
    report.to_csv(path, index=False)
    print(f"\n  Summary saved → {path}")
    return report

# ─────────────────────────────────────────────────────────────────────────────
# MAIN — runs everything in order
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "█"*60)
    print("  45S rRNA NOR Pipeline")
    print("  Micro-Tom × Arabidopsis × SL4.0 × Oryza sativa")
    print(f"  Working folder: {ALIGNMENT_DIR}")
    print(f"  Output folder:  {OUTPUT_DIR}")
    print("█"*60)

    # Validate that the working folder exists
    if not os.path.isdir(ALIGNMENT_DIR):
        print(f"\nERROR: Folder not found: {ALIGNMENT_DIR}")
        print("Please edit the ALIGNMENT_DIR variable at the top of this file.")
        sys.exit(1)

    refs                     = step1_download_references()
    all_data, all_hap        = step2_load_alignments()
    features                 = step3_features(all_data)
    features                 = step4_reference_similarity(features, refs)
    features, cluster_store  = step5_clustering(features)
    features, gap_results    = step6_gap_classifier(features)

    print("\n" + "="*60)
    print("  Generating figures (this may take a few minutes)...")
    print("="*60)

    step7_umap(features, refs)
    step8_conservation_heatmap(features, refs)
    step9_positional_analysis(features)
    step10_gap_importance(gap_results)
    step11_haplotype_analysis(all_hap, refs)
    step12_summary_report(features, refs)

    # Final summary of all files produced
    print("\n" + "█"*60)
    print("  PIPELINE COMPLETE — files produced:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        size = os.path.getsize(os.path.join(OUTPUT_DIR, f))
        print(f"    {f:<55} {size/1024:>6.1f} KB")
    print(f"\n  Open the pipeline_output folder to see all results.")
    print("█"*60)

if __name__ == "__main__":
    main()
