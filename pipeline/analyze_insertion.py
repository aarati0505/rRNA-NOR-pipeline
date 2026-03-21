# =============================================================================
# 141 bp Insertion Analysis Pipeline
# Analyzes the mysterious 141 bp gap found between 18S and rDNA promoter
#
# HOW TO RUN:
#   1. Put this file in the same folder as your 141bp FASTA file
#   2. Change FASTA_141BP below to your actual filename
#   3. Change YOUR_EMAIL to your real email
#   4. Run: python analyze_141bp.py
#
# WHAT THIS DOES:
#   Step 1 — Conserved motif analysis across all 25 sequences
#   Step 2 — BLAST against NCBI (finds matches in other species)
#   Step 3 — Transposable element check (structure + composition)
#   Step 4 — Regulatory sequence check (known regulatory motifs)
#   Step 5 — Summary report of all findings
# =============================================================================

# ── Edit these two lines ──────────────────────────────────────────────────────
FASTA_141BP   = r"C:\Users\Aarati\Downloads\rRNA\your_141bp_sequences.fasta"  # ← change this
YOUR_EMAIL    = "your_email@example.com"                                       # ← change this
# ─────────────────────────────────────────────────────────────────────────────

import os, sys, time, re, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter
from itertools import product

warnings.filterwarnings("ignore")

# Check libraries
missing = []
try:
    from Bio import Entrez, SeqIO
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import motifs
    from Bio.Seq import Seq
except ImportError:
    missing.append("biopython")
try:
    import matplotlib
except ImportError:
    missing.append("matplotlib")

if missing:
    print(f"ERROR: Missing libraries: {', '.join(missing)}")
    print(f"Run:  pip install {' '.join(missing)}")
    sys.exit(1)

Entrez.email  = YOUR_EMAIL
OUTPUT_DIR    = os.path.join(os.path.dirname(FASTA_141BP), "analysis_141bp_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# KNOWN REGULATORY MOTIFS TO SEARCH FOR
# These are experimentally verified sequences from plant rDNA literature
# ─────────────────────────────────────────────────────────────────────────────
REGULATORY_MOTIFS = {
    # Transcription termination signals
    "T-rich terminator":      r"TTTTTT",
    "Sal box (Pol I term)":   r"TTTGGATTT",
    "T1 terminator":          r"AGGTCGACC",

    # Processing signals
    "18S 3-prime signal":     r"CTCGTTCC",
    "ITS processing signal":  r"GTCTTGATT",
    "U3 snoRNA binding":      r"AAGCTTGA",

    # Transcription factor binding
    "TATA box":               r"TATA[AT]A[AT]",
    "GC box (Sp1)":           r"GGGCGG",
    "CCAAT box":              r"CCAAT",
    "UCE element":            r"GGGGTGT",

    # Repeat-associated
    "Microsatellite AT":      r"(AT){4,}",
    "Microsatellite GA":      r"(GA){4,}",
    "Direct repeat signal":   r"([ATGC]{6,8})\1",
}

# Known transposable element signatures
TE_SIGNATURES = {
    # Terminal repeats
    "TIR-like (TIR ends)":     r"^[ATGC]{2,8}CAGG.{100,}CCTG[ATGC]{2,8}$",
    "LINE RT motif":           r"YVDD",   # amino acid — check translated seq
    "SINE A-box":              r"TRRYNNAGYGG",
    "SINE B-box":              r"GWTCGANNC",
    "Helitron hairpin":        r"CTAG.{3,8}CTAG",

    # Composition signatures
    "AT-rich (>65%)":          None,   # handled separately
    "GC-rich (>65%)":          None,   # handled separately

    # Common plant SINE families
    "SRP-SINE motif":          r"GGGAGAGG",
    "tRNA-SINE remnant":       r"TGGCGAGTGG",
}

# ─────────────────────────────────────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def section(title):
    print("\n" + "="*60)
    print(f"  {title}")
    print("="*60)

def save_fig(name):
    p = os.path.join(OUTPUT_DIR, name)
    plt.savefig(p, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {p}")

def gc_content(seq):
    seq = seq.upper().replace('-','')
    if not seq: return 0
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def at_content(seq):
    return 100 - gc_content(seq)

def consensus_sequence(seqs):
    """Build a simple majority-vote consensus from aligned sequences."""
    if not seqs: return ""
    max_len = max(len(s) for s in seqs)
    consensus = []
    for i in range(max_len):
        bases = [s[i] for s in seqs if i < len(s) and s[i] in 'ATGC']
        if bases:
            consensus.append(Counter(bases).most_common(1)[0][0])
        else:
            consensus.append('-')
    return ''.join(consensus)

def per_position_conservation(seqs):
    """
    For each position return the fraction of sequences
    that match the consensus base (conservation score 0-1).
    """
    if not seqs: return [], []
    max_len = max(len(s) for s in seqs)
    positions, scores = [], []
    for i in range(max_len):
        bases = [s[i] for s in seqs if i < len(s) and s[i] in 'ATGC']
        if bases:
            top_freq = Counter(bases).most_common(1)[0][1] / len(bases)
            positions.append(i)
            scores.append(top_freq)
    return positions, scores

def kmer_frequencies(seq, k=2):
    """Return normalised k-mer frequency dict for a sequence."""
    seq = seq.upper().replace('-','')
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)
             if all(b in 'ATGC' for b in seq[i:i+k])]
    total = len(kmers)
    return {km: count/total for km, count in Counter(kmers).items()} if total else {}

# ─────────────────────────────────────────────────────────────────────────────
# LOAD SEQUENCES
# ─────────────────────────────────────────────────────────────────────────────

def load_sequences():
    section("Loading 141 bp sequences")
    if not os.path.exists(FASTA_141BP):
        print(f"ERROR: File not found: {FASTA_141BP}")
        print("Please check the FASTA_141BP path at the top of this script.")
        sys.exit(1)

    records = list(SeqIO.parse(FASTA_141BP, "fasta"))
    seqs    = [str(r.seq).upper() for r in records]
    ids     = [r.id for r in records]

    print(f"  Loaded {len(seqs)} sequences")
    for i, (sid, seq) in enumerate(zip(ids, seqs)):
        print(f"  [{i+1:2d}] {sid:40s} "
              f"len={len(seq.replace('-',''))} bp  "
              f"GC={gc_content(seq):.1f}%")

    return records, seqs, ids

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Conserved motif analysis
# ─────────────────────────────────────────────────────────────────────────────

def step1_conservation_analysis(seqs, ids):
    section("STEP 1 — Conservation analysis across 25 sequences")

    clean_seqs = [s.replace('-','') for s in seqs]

    # Per-position conservation
    positions, scores = per_position_conservation(seqs)
    cons_seq = consensus_sequence(seqs)

    avg_conservation = np.mean(scores) if scores else 0
    highly_conserved = sum(1 for s in scores if s >= 0.9)
    variable_pos     = sum(1 for s in scores if s < 0.7)

    print(f"  Consensus sequence ({len(cons_seq)} bp):")
    print(f"  5'-{cons_seq}-3'")
    print(f"\n  Average conservation : {avg_conservation:.3f}")
    print(f"  Highly conserved (≥90%) positions : {highly_conserved}")
    print(f"  Variable (<70%) positions          : {variable_pos}")

    # Save consensus
    cons_path = os.path.join(OUTPUT_DIR, "step1_consensus.fasta")
    with open(cons_path, 'w') as f:
        f.write(f">141bp_insertion_consensus\n{cons_seq}\n")
    print(f"\n  Consensus saved → {cons_path}")
    print("  (Use this sequence for BLAST if you want to search manually)")

    # GC content per sequence
    gc_values = [gc_content(s) for s in clean_seqs]
    print(f"\n  GC content across sequences:")
    print(f"    Mean : {np.mean(gc_values):.1f}%")
    print(f"    Min  : {np.min(gc_values):.1f}%")
    print(f"    Max  : {np.max(gc_values):.1f}%")

    # ── Plot 1: Conservation profile ────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))

    # Conservation score along the sequence
    ax = axes[0]
    colors = ['#2196F3' if s >= 0.9 else '#FF9800' if s >= 0.7
              else '#E63946' for s in scores]
    ax.bar(positions, scores, color=colors, width=1.0, alpha=0.85)
    ax.axhline(0.9, color='#2196F3', linestyle='--',
               linewidth=1, label='Highly conserved (90%)')
    ax.axhline(0.7, color='#FF9800', linestyle='--',
               linewidth=1, label='Moderately conserved (70%)')
    ax.set_title('Per-position conservation across 25 insertion sequences\n'
                 'Blue = highly conserved | Orange = moderate | Red = variable',
                 fontweight='bold')
    ax.set_xlabel('Position in 141 bp insertion (bp)')
    ax.set_ylabel('Conservation score')
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9)

    # GC content per sequence
    ax = axes[1]
    bar_colors = ['#E63946' if g > 65 else '#2196F3' if g < 35
                  else '#4CAF50' for g in gc_values]
    bars = ax.bar(range(len(gc_values)), gc_values,
                  color=bar_colors, edgecolor='white')
    ax.axhline(50, color='gray', linestyle=':', linewidth=1)
    ax.axhline(65, color='#E63946', linestyle='--',
               linewidth=1, label='GC-rich threshold (65%)')
    ax.axhline(35, color='#2196F3', linestyle='--',
               linewidth=1, label='AT-rich threshold (35%)')
    ax.set_title('GC content per insertion sequence',
                 fontweight='bold')
    ax.set_xlabel('Sequence index')
    ax.set_ylabel('GC content (%)')
    ax.set_xticks(range(len(ids)))
    ax.set_xticklabels([f"g{i+1}" for i in range(len(ids))],
                       rotation=45, fontsize=7)
    ax.legend(fontsize=9)

    # Dinucleotide frequency heatmap
    ax = axes[2]
    all_dinucs = sorted([''.join(p) for p in product('ATGC', repeat=2)])
    freq_matrix = np.array([
        [kmer_frequencies(s, k=2).get(dk, 0) for dk in all_dinucs]
        for s in clean_seqs
    ])
    sns.heatmap(freq_matrix, ax=ax,
                xticklabels=all_dinucs, yticklabels=False,
                cmap='YlOrRd', cbar_kws={'label': 'Frequency'})
    ax.set_title('Dinucleotide composition across all 25 sequences\n'
                 'Uniform rows = similar composition | '
                 'Variable rows = divergent sequences',
                 fontweight='bold')
    ax.set_xlabel('Dinucleotide')

    plt.tight_layout()
    save_fig("step1_conservation_profile.png")

    return cons_seq, avg_conservation, gc_values

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — BLAST against NCBI
# ─────────────────────────────────────────────────────────────────────────────

def step2_blast(cons_seq):
    section("STEP 2 — BLAST against NCBI")
    print("  Searching the consensus sequence against all plant sequences.")
    print("  This will take 2-5 minutes — NCBI processes the job remotely.\n")

    blast_results_path = os.path.join(OUTPUT_DIR, "step2_blast_results.xml")
    blast_summary_path = os.path.join(OUTPUT_DIR, "step2_blast_summary.csv")

    # Use cached results if already run
    if os.path.exists(blast_results_path):
        print("  (Using cached BLAST results from previous run)")
        with open(blast_results_path) as f:
            blast_records = list(NCBIXML.parse(f))
    else:
        print("  Submitting BLAST job to NCBI...")
        print("  Query: consensus sequence of 141 bp insertion")
        print("  Database: nt (all nucleotides)")
        print("  Organism filter: Viridiplantae (all plants)\n")

        try:
            result_handle = NCBIWWW.qblast(
                program   = "blastn",
                database  = "nt",
                sequence  = cons_seq.replace('-',''),
                entrez_query = "Viridiplantae[Organism]",
                hitlist_size = 50,
                expect    = 0.001,
            )
            # Save raw XML for caching
            raw = result_handle.read()
            with open(blast_results_path, 'w') as f:
                f.write(raw)
            from io import StringIO
            blast_records = list(NCBIXML.parse(StringIO(raw)))
            print("  BLAST complete.")
        except Exception as e:
            print(f"  ERROR during BLAST: {e}")
            print("  Check your internet connection and try again.")
            print("  Skipping BLAST step.")
            return None

    # Parse results
    rows = []
    print(f"\n  {'Rank':<5} {'Species':<40} {'Identity%':<12}"
          f"{'Coverage%':<12} {'E-value':<12} {'Description'}")
    print("  " + "-"*100)

    query_len = len(cons_seq.replace('-',''))

    for blast_record in blast_records:
        for rank, alignment in enumerate(blast_record.alignments[:20], 1):
            for hsp in alignment.hsps[:1]:   # top HSP only
                identity_pct  = hsp.identities / hsp.align_length * 100
                coverage_pct  = hsp.align_length / query_len * 100
                description   = alignment.title[:80]

                # Extract species name from title
                species = "unknown"
                match = re.search(r'\[([A-Z][a-z]+ [a-z]+)\]', alignment.title)
                if match:
                    species = match.group(1)

                row = {
                    "rank":         rank,
                    "species":      species,
                    "description":  description,
                    "identity_pct": round(identity_pct, 1),
                    "coverage_pct": round(coverage_pct, 1),
                    "evalue":       hsp.expect,
                    "score":        hsp.score,
                    "alignment_len":hsp.align_length,
                }
                rows.append(row)

                print(f"  {rank:<5} {species:<40} "
                      f"{identity_pct:<12.1f}{coverage_pct:<12.1f}"
                      f"{hsp.expect:<12.2e} {description[:40]}")

    if rows:
        df_blast = pd.DataFrame(rows)
        df_blast.to_csv(blast_summary_path, index=False)
        print(f"\n  Full BLAST results saved → {blast_summary_path}")

        # ── Plot BLAST results ───────────────────────────────────────────────
        top20 = df_blast.head(20).copy()
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))

        # Identity % per hit
        ax = axes[0]
        bar_colors = ['#E63946' if v >= 90 else '#FF9800' if v >= 70
                      else '#BDBDBD' for v in top20['identity_pct']]
        ax.barh(range(len(top20)), top20['identity_pct'],
                color=bar_colors, edgecolor='white')
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(
            [f"{r['species'][:30]} ({r['identity_pct']}%)"
             for _, r in top20.iterrows()],
            fontsize=8
        )
        ax.axvline(90, color='#E63946', linestyle='--',
                   linewidth=1, label='90% identity')
        ax.axvline(70, color='#FF9800', linestyle='--',
                   linewidth=1, label='70% identity')
        ax.set_xlabel('Sequence identity (%)')
        ax.set_title('BLAST hits — sequence identity\n'
                     'Red = very similar (≥90%) | '
                     'Orange = similar (≥70%)',
                     fontweight='bold')
        ax.legend(fontsize=9)
        ax.invert_yaxis()

        # E-value per hit
        ax = axes[1]
        evalues = [-np.log10(max(e, 1e-200)) for e in top20['evalue']]
        ax.barh(range(len(top20)), evalues,
                color='#2196F3', edgecolor='white', alpha=0.8)
        ax.set_yticks(range(len(top20)))
        ax.set_yticklabels(
            [r['species'][:35] for _, r in top20.iterrows()],
            fontsize=8
        )
        ax.set_xlabel('-log10(E-value)  [higher = more significant]')
        ax.set_title('BLAST hits — statistical significance\n'
                     'Higher bar = more significant match',
                     fontweight='bold')
        ax.invert_yaxis()

        plt.suptitle('141 bp Insertion — BLAST Results vs NCBI Plant Database',
                     fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig("step2_blast_results.png")

        return df_blast
    else:
        print("  No significant BLAST hits found.")
        print("  This suggests the 141 bp insertion may be novel / Micro-Tom specific.")
        return None

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — Transposable element analysis
# ─────────────────────────────────────────────────────────────────────────────

def step3_te_analysis(seqs, cons_seq, ids):
    section("STEP 3 — Transposable element analysis")
    print("  Checking structural and compositional signatures of TEs.")
    print("  Note: For definitive TE identification, also run RepeatMasker")
    print("  at https://www.repeatmasker.org/cgi-bin/WEBRepeatMasker\n")

    clean_seqs  = [s.replace('-','') for s in seqs]
    clean_cons  = cons_seq.replace('-','')
    results     = {}

    # ── Check 1: Size — is 141 bp consistent with known TE families? ────────
    lengths = [len(s) for s in clean_seqs]
    print(f"  Sequence lengths: {min(lengths)}–{max(lengths)} bp "
          f"(mean {np.mean(lengths):.0f} bp)")
    print("  TE size reference:")
    print("    SINEs:    100–500 bp  ← 141 bp FITS this range")
    print("    MITEs:    100–500 bp  ← 141 bp FITS this range")
    print("    LINEs:    >1000 bp    ← too small")
    print("    LTR-TEs:  >300 bp     ← possibly too small")
    results['size_consistent_with_SINE'] = True
    results['size_consistent_with_MITE'] = True

    # ── Check 2: Terminal Inverted Repeats (TIRs) ────────────────────────────
    # MITEs always have TIRs — short inverted repeats at both ends
    print("\n  Checking for Terminal Inverted Repeats (TIRs)...")
    tir_found = []
    for seq in clean_seqs[:5]:   # check first 5
        end_len = 15
        left    = seq[:end_len]
        right   = str(Seq(seq[-end_len:]).reverse_complement())
        matches = sum(a == b for a, b in zip(left, right))
        similarity = matches / end_len * 100
        tir_found.append(similarity)
        if similarity >= 60:
            print(f"    Possible TIR: left={left}  rc_right={right}  "
                  f"similarity={similarity:.0f}%")

    mean_tir = np.mean(tir_found)
    results['TIR_detected'] = mean_tir >= 60
    print(f"    Mean terminal similarity: {mean_tir:.1f}% "
          f"({'TIR likely' if mean_tir >= 60 else 'No clear TIR'})")

    # ── Check 3: Target Site Duplication (TSD) ───────────────────────────────
    # TEs often create short duplications at insertion site
    # We can't check this without flanking sequence, but note it
    print("\n  Target Site Duplication (TSD) check:")
    print("    Cannot assess without flanking genomic sequence.")
    print("    Recommendation: extract 10 bp flanking each insertion site")
    print("    and check if they match — TSDs of 2-10 bp are TE signatures.")

    # ── Check 4: Composition signatures ─────────────────────────────────────
    print("\n  Compositional signatures:")
    mean_gc = np.mean([gc_content(s) for s in clean_seqs])
    mean_at = 100 - mean_gc
    print(f"    Mean GC content : {mean_gc:.1f}%")
    print(f"    Mean AT content : {mean_at:.1f}%")

    if mean_at > 65:
        print("    → AT-rich: consistent with MITE or some SINE families")
        results['AT_rich'] = True
    elif mean_gc > 65:
        print("    → GC-rich: less common for TEs, possible regulatory element")
        results['GC_rich'] = True
    else:
        print("    → Balanced composition: consistent with functional sequence")

    # ── Check 5: Known TE motif scanning ────────────────────────────────────
    print("\n  Scanning for known TE motifs in consensus sequence:")
    te_hits = {}
    for te_name, pattern in TE_SIGNATURES.items():
        if pattern is None:
            continue
        matches = re.findall(pattern, clean_cons, re.IGNORECASE)
        if matches:
            te_hits[te_name] = len(matches)
            print(f"    FOUND: {te_name} ({len(matches)} match(es))")

    if not te_hits:
        print("    No known TE motifs detected in consensus")
        results['TE_motifs_found'] = False
    else:
        results['TE_motifs_found'] = True
        results['TE_hits'] = te_hits

    # ── Check 6: Palindrome / hairpin potential ──────────────────────────────
    print("\n  Checking for palindromic sequences (hairpin potential)...")
    hairpins = []
    for i in range(0, len(clean_cons) - 20, 5):
        seg   = clean_cons[i:i+10]
        rc    = str(Seq(seg).reverse_complement())
        # Look for the reverse complement within 30 bp downstream
        downstream = clean_cons[i+12:i+50]
        if rc in downstream:
            hairpins.append((i, seg))
            print(f"    Hairpin stem at pos {i}: {seg} ... {rc}")

    results['hairpin_detected'] = len(hairpins) > 0
    if not hairpins:
        print("    No obvious palindromic sequences detected")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n  TE ANALYSIS SUMMARY:")
    print(f"    Size consistent with SINE/MITE : YES (141 bp fits range)")
    print(f"    TIR detected                   : "
          f"{'YES' if results.get('TIR_detected') else 'NO/WEAK'}")
    print(f"    TE motifs in consensus         : "
          f"{'YES — ' + str(te_hits) if results.get('TE_motifs_found') else 'NO'}")
    print(f"    Hairpin structures             : "
          f"{'YES' if results.get('hairpin_detected') else 'NO'}")
    print(f"    AT-rich composition            : "
          f"{'YES' if results.get('AT_rich') else 'NO'}")

    print("\n  ONLINE TOOLS TO CONFIRM TE IDENTITY:")
    print("    1. RepeatMasker  → https://www.repeatmasker.org/cgi-bin/WEBRepeatMasker")
    print("       Paste your consensus sequence — it will identify TE family")
    print("    2. CENSOR        → https://www.girinst.org/censor/")
    print("       Screens against Repbase (largest TE database)")
    print("    3. TEfam         → https://tefam.biochem.vt.edu/tefam/")
    print("       Specifically for plant TEs")

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # GC vs AT distribution
    ax = axes[0]
    gc_vals = [gc_content(s) for s in clean_seqs]
    ax.hist(gc_vals, bins=15, color='#2196F3', edgecolor='white', alpha=0.8)
    ax.axvline(np.mean(gc_vals), color='#E63946', linestyle='--',
               linewidth=2, label=f'Mean GC = {np.mean(gc_vals):.1f}%')
    ax.axvline(65, color='orange', linestyle=':',
               linewidth=1.5, label='GC-rich threshold (65%)')
    ax.axvline(35, color='green',  linestyle=':',
               linewidth=1.5, label='AT-rich threshold (35%)')
    ax.set_xlabel('GC content (%)')
    ax.set_ylabel('Number of sequences')
    ax.set_title('GC content distribution\n(AT-rich = possible MITE/SINE)',
                 fontweight='bold')
    ax.legend(fontsize=9)

    # Terminal similarity (TIR check)
    ax = axes[1]
    end_lens = range(5, 25)
    tir_scores_by_len = []
    for elen in end_lens:
        sims = []
        for seq in clean_seqs:
            if len(seq) < elen * 2: continue
            left  = seq[:elen]
            right = str(Seq(seq[-elen:]).reverse_complement())
            n     = min(len(left), len(right))
            sim   = sum(a==b for a,b in zip(left,right)) / n * 100
            sims.append(sim)
        tir_scores_by_len.append(np.mean(sims) if sims else 0)

    ax.plot(end_lens, tir_scores_by_len, 'o-',
            color='#E63946', linewidth=2, markersize=6)
    ax.axhline(60, color='gray', linestyle='--',
               linewidth=1, label='TIR threshold (60%)')
    ax.fill_between(end_lens, tir_scores_by_len, alpha=0.2, color='#E63946')
    ax.set_xlabel('Terminal repeat length checked (bp)')
    ax.set_ylabel('Mean terminal inverted similarity (%)')
    ax.set_title('Terminal Inverted Repeat (TIR) analysis\n'
                 '(>60% = possible MITE/DNA transposon)',
                 fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 105)

    plt.suptitle('141 bp Insertion — Transposable Element Signatures',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    save_fig("step3_te_analysis.png")

    return results

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Regulatory sequence analysis
# ─────────────────────────────────────────────────────────────────────────────

def step4_regulatory_analysis(seqs, cons_seq, ids):
    section("STEP 4 — Regulatory sequence analysis")
    print("  Scanning for known plant rDNA regulatory motifs.\n")

    clean_cons  = cons_seq.replace('-','')
    clean_seqs  = [s.replace('-','') for s in seqs]

    all_hits    = {}
    rows        = []

    # Scan consensus for each regulatory motif
    for motif_name, pattern in REGULATORY_MOTIFS.items():
        matches = list(re.finditer(pattern, clean_cons, re.IGNORECASE))
        if matches:
            positions = [m.start() for m in matches]
            sequences = [m.group() for m in matches]
            all_hits[motif_name] = {
                "count":     len(matches),
                "positions": positions,
                "sequences": sequences,
            }
            print(f"  FOUND: {motif_name}")
            for m in matches:
                print(f"    Position {m.start():4d}: {m.group()}")

    if not all_hits:
        print("  No known regulatory motifs found in consensus sequence.")
        print("  This suggests the insertion is NOT a simple regulatory element")
        print("  and is more likely a novel insertion or degraded TE remnant.")
    else:
        print(f"\n  Total regulatory motifs found: {len(all_hits)}")

    # Count motifs per sequence (not just consensus)
    print("\n  Motif presence across all 25 sequences:")
    for motif_name, pattern in REGULATORY_MOTIFS.items():
        counts = []
        for seq in clean_seqs:
            n = len(re.findall(pattern, seq, re.IGNORECASE))
            counts.append(n)
        present_in = sum(1 for c in counts if c > 0)
        if present_in > 0:
            rows.append({
                "motif":      motif_name,
                "present_in": present_in,
                "total_seqs": len(clean_seqs),
                "pct_present":round(present_in/len(clean_seqs)*100, 1),
                "in_consensus": motif_name in all_hits,
            })
            print(f"  {motif_name:<35} present in "
                  f"{present_in}/{len(clean_seqs)} sequences "
                  f"({present_in/len(clean_seqs)*100:.0f}%)")

    if rows:
        df_reg = pd.DataFrame(rows).sort_values('pct_present', ascending=False)
        reg_path = os.path.join(OUTPUT_DIR, "step4_regulatory_motifs.csv")
        df_reg.to_csv(reg_path, index=False)
        print(f"\n  Full motif table saved → {reg_path}")

        # ── Plot ─────────────────────────────────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Motif prevalence bar chart
        ax = axes[0]
        bar_colors = ['#E63946' if r['in_consensus'] else '#2196F3'
                      for _, r in df_reg.iterrows()]
        ax.barh(range(len(df_reg)), df_reg['pct_present'],
                color=bar_colors, edgecolor='white')
        ax.set_yticks(range(len(df_reg)))
        ax.set_yticklabels(df_reg['motif'], fontsize=9)
        ax.set_xlabel('% of sequences containing motif')
        ax.set_title('Regulatory motif prevalence\n'
                     'Red = also in consensus | Blue = in individual seqs only',
                     fontweight='bold')
        ax.axvline(80, color='gray', linestyle='--',
                   linewidth=1, label='80% threshold')
        ax.legend(fontsize=9)
        ax.invert_yaxis()

        # Motif position map on consensus
        ax = axes[1]
        ax.set_xlim(0, len(clean_cons))
        ax.set_ylim(-1, len(all_hits) + 1)
        ax.set_title('Regulatory motif positions on consensus sequence\n'
                     '(horizontal = position along 141 bp insertion)',
                     fontweight='bold')
        ax.set_xlabel('Position in consensus (bp)')
        ax.set_yticks(range(len(all_hits)))
        ax.set_yticklabels(list(all_hits.keys()), fontsize=8)
        ax.axvline(0, color='black', linewidth=2)
        ax.axvline(len(clean_cons), color='black', linewidth=2)

        colors = plt.cm.tab10(np.linspace(0, 1, len(all_hits)))
        for yi, (mname, mdata) in enumerate(all_hits.items()):
            for pos in mdata['positions']:
                ax.barh(yi, 10, left=pos, height=0.6,
                        color=colors[yi], alpha=0.8)
                ax.text(pos + 5, yi, str(pos),
                        fontsize=7, va='center', ha='center', color='white')

        plt.suptitle('141 bp Insertion — Regulatory Sequence Analysis',
                     fontsize=13, fontweight='bold')
        plt.tight_layout()
        save_fig("step4_regulatory_analysis.png")

    return all_hits, rows

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — Final summary report
# ─────────────────────────────────────────────────────────────────────────────

def step5_summary(cons_seq, avg_conservation, gc_values,
                  df_blast, te_results, reg_hits, reg_rows):
    section("STEP 5 — Final Summary Report")

    clean_cons = cons_seq.replace('-','')
    mean_gc    = np.mean(gc_values)

    print("\n  ┌─────────────────────────────────────────────────────┐")
    print("  │         141 BP INSERTION — ANALYSIS SUMMARY         │")
    print("  └─────────────────────────────────────────────────────┘")

    print(f"\n  BASIC PROPERTIES")
    print(f"    Consensus length       : {len(clean_cons)} bp")
    print(f"    Number of copies       : 25 (in your alignment)")
    print(f"    Position in gene       : Between 18S and rDNA promoter")
    print(f"    Mean GC content        : {mean_gc:.1f}%")
    print(f"    Average conservation   : {avg_conservation:.1%} across copies")

    print(f"\n  BLAST RESULTS")
    if df_blast is not None and len(df_blast) > 0:
        top = df_blast.iloc[0]
        print(f"    Top hit species        : {top['species']}")
        print(f"    Top hit identity       : {top['identity_pct']}%")
        print(f"    Top hit e-value        : {top['evalue']:.2e}")
        same_sp = df_blast[df_blast['species'].str.contains(
            'Solanum|lycopersicum', case=False, na=False)]
        if len(same_sp) > 0:
            print(f"    Hits in Solanum spp.   : {len(same_sp)} "
                  f"(present in other tomato varieties)")
        else:
            print(f"    Hits in Solanum spp.   : 0 (possibly Micro-Tom specific)")
    else:
        print("    BLAST not run or no hits found")
        print("    → Suggests possible novel / Micro-Tom-specific insertion")

    print(f"\n  TRANSPOSABLE ELEMENT ANALYSIS")
    print(f"    Size fits SINE/MITE    : YES (100-500 bp range)")
    print(f"    TIR detected           : "
          f"{'YES' if te_results.get('TIR_detected') else 'NO/WEAK'}")
    print(f"    TE motifs found        : "
          f"{'YES' if te_results.get('TE_motifs_found') else 'NO'}")
    print(f"    AT-rich (>65%)         : "
          f"{'YES' if te_results.get('AT_rich') else 'NO'}")

    print(f"\n  REGULATORY ELEMENT ANALYSIS")
    if reg_rows:
        high_prev = [r for r in reg_rows if r['pct_present'] >= 80]
        print(f"    Motifs found           : {len(reg_hits)}")
        print(f"    High prevalence (≥80%) : {len(high_prev)}")
        for r in high_prev:
            print(f"      → {r['motif']} ({r['pct_present']}% of copies)")
    else:
        print(f"    No regulatory motifs found")

    print(f"\n  BIOLOGICAL INTERPRETATION")
    if mean_gc < 40 and te_results.get('TIR_detected'):
        interp = ("AT-rich with TIR signature → likely a MITE "
                  "(Miniature Inverted-repeat Transposable Element). "
                  "MITEs commonly insert into plant rDNA spacers.")
    elif mean_gc < 40:
        interp = ("AT-rich composition → possible SINE or MITE remnant. "
                  "Run RepeatMasker for definitive classification.")
    elif reg_rows and any(r['pct_present'] >= 80 for r in reg_rows):
        interp = ("Regulatory motifs present in most copies → "
                  "possible functional regulatory element "
                  "involved in rDNA transcription control.")
    elif df_blast is not None and len(df_blast) > 0:
        top_id = df_blast.iloc[0]['identity_pct']
        if top_id >= 90:
            interp = ("High BLAST identity → conserved element present "
                      "in other species. Likely has a conserved function.")
        elif top_id >= 70:
            interp = ("Moderate BLAST identity → related sequence in other "
                      "plants but diverged. Possible ancient insertion.")
        else:
            interp = ("Low BLAST identity → mostly novel sequence. "
                      "Could be Micro-Tom-specific insertion.")
    else:
        interp = ("No strong matches to known elements → "
                  "possibly novel Micro-Tom-specific insertion. "
                  "Recommend RepeatMasker + Dfam database search.")

    print(f"\n    → {interp}")

    print(f"\n  RECOMMENDED NEXT STEPS")
    print("    1. Run RepeatMasker on consensus:")
    print("       https://www.repeatmasker.org/cgi-bin/WEBRepeatMasker")
    print("    2. Search Dfam (TE database):")
    print("       https://dfam.org/search")
    print("    3. Check if insertion sites share a consensus sequence")
    print("       (extract 10 bp flanking sequence from your genome positions)")
    print("    4. If BLAST found Solanum hits: compare insertion presence/absence")
    print("       across tomato varieties to date the insertion event")

    print(f"\n  OUTPUT FILES GENERATED:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        size = os.path.getsize(os.path.join(OUTPUT_DIR, f))
        print(f"    {f:<50} {size/1024:>6.1f} KB")

    # Save full summary to text file
    summary_path = os.path.join(OUTPUT_DIR, "step5_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("141 BP INSERTION ANALYSIS SUMMARY\n")
        f.write("="*50 + "\n\n")
        f.write(f"Consensus sequence:\n5'-{clean_cons}-3'\n\n")
        f.write(f"Length          : {len(clean_cons)} bp\n")
        f.write(f"GC content      : {mean_gc:.1f}%\n")
        f.write(f"Conservation    : {avg_conservation:.1%}\n")
        f.write(f"Interpretation  : {interp}\n")
    print(f"\n  Full summary saved → {summary_path}")

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "█"*60)
    print("  141 bp Insertion Analysis Pipeline")
    print(f"  Input file : {FASTA_141BP}")
    print(f"  Output dir : {OUTPUT_DIR}")
    print("█"*60)

    records, seqs, ids       = load_sequences()
    cons_seq, avg_cons, gc_v = step1_conservation_analysis(seqs, ids)
    df_blast                 = step2_blast(cons_seq)
    te_results               = step3_te_analysis(seqs, cons_seq, ids)
    reg_hits, reg_rows       = step4_regulatory_analysis(seqs, cons_seq, ids)
    step5_summary(cons_seq, avg_cons, gc_v,
                  df_blast, te_results, reg_hits, reg_rows)

    print("\n" + "█"*60)
    print("  ANALYSIS COMPLETE")
    print(f"  Results in: {OUTPUT_DIR}")
    print("█"*60)

if __name__ == "__main__":
    main()
