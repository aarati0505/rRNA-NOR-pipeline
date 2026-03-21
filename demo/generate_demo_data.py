"""
generate_demo_data.py
─────────────────────
Generates synthetic 45S rRNA-like FASTA files for testing the pipeline.

These sequences are NOT real biological data — they are randomly generated
with composition profiles loosely resembling plant rDNA regions.
They exist solely so users can verify the pipeline runs correctly
before using their own data.

Run:
    python demo/generate_demo_data.py
"""

import numpy as np
import os

np.random.seed(42)

OUTPUT_DIR  = os.path.dirname(os.path.abspath(__file__))
N_SEQUENCES = 50       # sequences per alignment file
GENE_LENGTH = 10000    # approximate length of one rDNA unit

# Approximate nucleotide compositions for each functional region
# [A, T, G, C] probabilities
REGION_PROFILES = {
    "5ETS":     ([0.28, 0.22, 0.27, 0.23], 500),   # (composition, length bp)
    "18S":      ([0.25, 0.25, 0.25, 0.25], 1800),
    "ITS1":     ([0.30, 0.20, 0.25, 0.25], 300),
    "5_8S":     ([0.26, 0.24, 0.26, 0.24], 160),
    "ITS2":     ([0.28, 0.22, 0.25, 0.25], 250),
    "25S":      ([0.25, 0.25, 0.25, 0.25], 3400),
    "3ETS":     ([0.27, 0.23, 0.26, 0.24], 790),
    "IGS":      ([0.30, 0.30, 0.20, 0.20], 2800),
}

BASES = ['A', 'T', 'G', 'C']

def make_region(composition, length, mutation_rate=0.02):
    """Generate a random DNA sequence with given nucleotide composition."""
    seq = np.random.choice(BASES, size=length, p=composition)
    # Add random mutations
    mask = np.random.random(length) < mutation_rate
    seq[mask] = np.random.choice(BASES, size=mask.sum())
    return ''.join(seq)

def make_gene(with_insertion=False, with_large_gap=False):
    """
    Build one synthetic rDNA unit from regional profiles.
    Optionally add a small insertion or large gap to simulate
    the variation seen in real rDNA arrays.
    """
    parts = []
    for region, (comp, length) in REGION_PROFILES.items():
        # Add some length variation between copies
        var_len = max(50, int(length * np.random.uniform(0.95, 1.05)))
        parts.append(make_region(comp, var_len))

    gene = ''.join(parts)

    # Optionally insert a small sequence (simulates the 141bp-type insertion)
    if with_insertion:
        # Insert a GC-rich 141 bp sequence between ITS2 and 25S
        insert_pos = sum(v[1] for k, v in REGION_PROFILES.items()
                         if k in ['5ETS','18S','ITS1','5_8S','ITS2'])
        insert_seq = make_region([0.18, 0.18, 0.32, 0.32], 141)
        gene = gene[:insert_pos] + insert_seq + gene[insert_pos:]

    # Optionally add gaps to simulate large insertions in alignment
    if with_large_gap:
        gap_pos = int(len(gene) * 0.3)
        gap_len = np.random.randint(800, 3000)
        gap_seq = make_region([0.25, 0.25, 0.25, 0.25], gap_len)
        gene = gene[:gap_pos] + gap_seq + gene[gap_pos:]

    return gene

def make_sequence_id(i, start_pos):
    """Generate a sequence ID mimicking real NOR region IDs."""
    end_pos = start_pos + np.random.randint(8000, 12000)
    return f"NOR2_region_(15.8_Mb)/{start_pos}_{end_pos}"

def write_fasta(sequences, ids, filepath):
    """Write sequences to a FASTA file."""
    with open(filepath, 'w') as f:
        for seq_id, seq in zip(ids, sequences):
            f.write(f">{seq_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    print(f"  Written: {filepath}  ({len(sequences)} sequences)")

def main():
    print("Generating synthetic demo sequences for rRNA NOR pipeline testing")
    print("NOTE: These are NOT real biological sequences\n")

    # ── Alignment 1: all gaps present ────────────────────────────────────────
    seqs1, ids1 = [], []
    pos = 244
    for i in range(N_SEQUENCES):
        has_insertion  = (i % 20 == 0)    # ~5% have the small insertion
        has_large_gap  = (i % 10 == 0)    # ~10% have a large gap
        gene = make_gene(with_insertion=has_insertion,
                         with_large_gap=has_large_gap)
        seqs1.append(gene)
        ids1.append(make_sequence_id(i, pos))
        pos += np.random.randint(9000, 11000)

    write_fasta(seqs1, ids1,
                os.path.join(OUTPUT_DIR, "Alignment_1_all_gaps_DEMO.fasta"))

    # ── Alignment 2: with 18S insertion variant ───────────────────────────────
    seqs2, ids2 = [], []
    pos = 244
    for i in range(N_SEQUENCES):
        has_insertion = (i % 20 == 0)
        gene = make_gene(with_insertion=has_insertion, with_large_gap=False)
        seqs2.append(gene)
        ids2.append(make_sequence_id(i, pos))
        pos += np.random.randint(9000, 11000)

    write_fasta(seqs2, ids2,
                os.path.join(OUTPUT_DIR, "Alignment_2_18S_gap_DEMO.fasta"))

    # ── Alignment 3: small gaps only ─────────────────────────────────────────
    seqs3, ids3 = [], []
    pos = 244
    for i in range(N_SEQUENCES):
        gene = make_gene(with_insertion=False, with_large_gap=False)
        seqs3.append(gene)
        ids3.append(make_sequence_id(i, pos))
        pos += np.random.randint(9000, 11000)

    write_fasta(seqs3, ids3,
                os.path.join(OUTPUT_DIR, "Alignment_3_small_gaps_DEMO.fasta"))

    # ── Combined demo file ────────────────────────────────────────────────────
    all_seqs = seqs1 + seqs2 + seqs3
    all_ids  = ids1  + ids2  + ids3
    write_fasta(all_seqs, all_ids,
                os.path.join(OUTPUT_DIR, "demo_sequences.fasta"))

    # ── Haplotype demo files ──────────────────────────────────────────────────
    for hap_num in range(1, 4):
        hap_seqs, hap_ids = [], []
        for h in range(9):
            gene = make_gene(with_insertion=(h == 3),
                             with_large_gap=False)
            hap_seqs.append(gene)
            hap_ids.append(f"H{h+1:04d}/without_removing_{hap_num}_DEMO")
        fname = (f"mafft_aligned_Summary_Filtered_Haplotype_"
                 f"file_without_removing_{hap_num}_DEMO.fasta")
        write_fasta(hap_seqs, hap_ids, os.path.join(OUTPUT_DIR, fname))

    print("\nDemo files generated successfully.")
    print("To run the pipeline on these files:")
    print("  1. Set ALIGNMENT_DIR to the demo/ folder in rRNA_pipeline.py")
    print("  2. Update ALIGNMENT_FILES dict to use the _DEMO filenames")
    print("  3. Run: python pipeline/rRNA_pipeline.py")

if __name__ == "__main__":
    main()
