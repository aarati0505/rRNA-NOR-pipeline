# =============================================================================
# Haplotype Analysis Scripts for 45S rRNA Gene Sequences
# Two methods: standard R libraries and geneHapR package
#
# HOW TO RUN:
#   1. Open RStudio
#   2. Install required packages (run the install block once)
#   3. Run Method 1 or Method 2
# =============================================================================

# ── Install packages (run once) ───────────────────────────────────────────────
# install.packages("ape")
# install.packages("pegas")
# install.packages("geneHapR")
# install.packages("Biostrings")   # from Bioconductor:
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("Biostrings")

# =============================================================================
# METHOD 1 — Standard R libraries (ape + pegas)
# Best for: quick haplotyping, haplotype network plots
# Limitation: may produce many haplotypes if gaps are present
# =============================================================================

method1_haplotyping <- function(fasta_path = NULL) {

  library(ape)
  library(pegas)

  # Load FASTA file — opens file picker if no path given
  if (is.null(fasta_path)) {
    fasta_path <- file.choose()
  }

  cat("Loading sequences from:", fasta_path, "\n")
  dna_sequences <- read.dna(fasta_path, format = "fasta")
  cat("Loaded", nrow(dna_sequences), "sequences\n")

  # Identify haplotypes
  cat("Identifying haplotypes...\n")
  haplotypes  <- haplotype(dna_sequences)
  hap_summary <- haploNet(haplotypes)

  # Print results
  print(haplotypes)
  cat("\nSummary:\n")
  summary(haplotypes)

  # Plot haplotype network
  # Circle size proportional to haplotype frequency
  plot(hap_summary,
       size        = attr(hap_summary, "freq"),
       scale.ratio = 0.5,
       pie         = attr(hap_summary, "freq"))

  # Save haplotype sequences
  write.dna(haplotypes,
            file   = "haplotypes_output_method1.fasta",
            format = "fasta")
  cat("Haplotype sequences saved to: haplotypes_output_method1.fasta\n")

  return(haplotypes)
}

# Run Method 1:
# hap1 <- method1_haplotyping()
# Or with a specific path:
# hap1 <- method1_haplotyping("C:/path/to/your/alignment.fasta")


# =============================================================================
# METHOD 2 — geneHapR package
# Best for: controlling gap tolerance, more flexible filtering
# Produces: haplotype CSV + assignment CSV
# =============================================================================

method2_haplotyping <- function(fasta_path, max_gap_freq = 0.8) {

  library(geneHapR)
  library(Biostrings)

  cat("Loading sequences from:", fasta_path, "\n")
  dna_seqs <- readDNAStringSet(fasta_path)
  cat("Loaded", length(dna_seqs), "sequences\n")

  # Assign names if missing or empty
  if (is.null(names(dna_seqs)) ||
      any(is.na(names(dna_seqs))) ||
      any(names(dna_seqs) == "")) {
    cat("Warning: missing sequence names — assigning seq1, seq2...\n")
    names(dna_seqs) <- paste0("seq", seq_along(dna_seqs))
  }

  # Convert to matrix for gap filtering
  cat("Converting to alignment matrix...\n")
  seq_matrix <- as.matrix(dna_seqs)

  # Remove columns where too many sequences have gaps
  cat(sprintf("Filtering columns with >%.0f%% gaps...\n", max_gap_freq * 100))
  gap_freq            <- colSums(seq_matrix == "-") / nrow(seq_matrix)
  keep_cols           <- which(gap_freq <= max_gap_freq)
  seq_matrix_filtered <- seq_matrix[, keep_cols]
  cat(sprintf("Kept %d / %d alignment columns\n",
              length(keep_cols), ncol(seq_matrix)))

  # Convert filtered matrix back to DNAStringSet
  filtered_seqs        <- DNAStringSet(apply(
    seq_matrix_filtered, 1, paste0, collapse = ""
  ))
  names(filtered_seqs) <- names(dna_seqs)

  # Identify haplotypes
  cat("Identifying haplotypes...\n")
  hap_result <- seqs2hap(filtered_seqs, max.gap.freq = 2)

  # Print results
  cat("\nHaplotype results:\n")
  print(hap_result$hap)
  cat("\nHaplotype information:\n")
  print(hap_result$hapInfo)

  # Save results
  write.csv(hap_result$hap,
            "haplotypes_method2.csv",
            row.names = FALSE)
  write.csv(hap_result$hapInfo,
            "hap_assignments_method2.csv",
            row.names = FALSE)
  cat("Results saved to: haplotypes_method2.csv and hap_assignments_method2.csv\n")

  return(hap_result)
}

# Run Method 2:
# hap2 <- method2_haplotyping("C:/path/to/Alignment_3.fasta")
# Adjust max_gap_freq: lower = stricter (fewer columns kept)
# hap2 <- method2_haplotyping("C:/path/to/alignment.fasta", max_gap_freq = 0.5)


# =============================================================================
# SEQUENCE CLUSTERING — K-Means alternative to haplotyping
# Use this if the above methods produce too many haplotypes
# due to gaps. Clustering groups sequences by similarity
# rather than exact haplotype identity.
# =============================================================================

cluster_sequences <- function(fasta_path, n_clusters = 9) {

  library(ape)
  library(stats)

  cat("Loading sequences...\n")
  dna_seqs <- read.dna(fasta_path, format = "fasta")

  # Compute pairwise distance matrix
  cat("Computing pairwise distances (this may take a few minutes)...\n")
  dist_matrix <- dist.dna(dna_seqs,
                           model    = "raw",
                           pairwise.deletion = TRUE)

  # K-Means clustering on distance matrix
  cat(sprintf("Clustering into %d groups...\n", n_clusters))
  km <- kmeans(as.matrix(dist_matrix),
               centers  = n_clusters,
               nstart   = 20,
               iter.max = 500)

  # Report cluster sizes
  cat("\nCluster sizes:\n")
  print(table(km$cluster))

  # Save cluster assignments
  cluster_df <- data.frame(
    sequence_id = rownames(as.matrix(dist_matrix)),
    cluster     = km$cluster
  )
  write.csv(cluster_df, "cluster_assignments.csv", row.names = FALSE)
  cat("Cluster assignments saved to: cluster_assignments.csv\n")

  return(km)
}

# Run clustering:
# km <- cluster_sequences("C:/path/to/alignment.fasta", n_clusters = 9)


# =============================================================================
# FIND OPTIMAL CLUSTER NUMBER (elbow method)
# Run this before cluster_sequences() to choose n_clusters
# =============================================================================

find_optimal_clusters <- function(fasta_path,
                                  max_k   = 20,
                                  sample_n = 200) {

  library(ape)

  cat("Loading sequences...\n")
  dna_seqs <- read.dna(fasta_path, format = "fasta")

  # Subsample for speed if many sequences
  if (nrow(dna_seqs) > sample_n) {
    cat(sprintf("Subsampling %d sequences for speed...\n", sample_n))
    idx      <- sample(nrow(dna_seqs), sample_n)
    dna_seqs <- dna_seqs[idx, ]
  }

  cat("Computing distance matrix...\n")
  dist_mat <- as.matrix(dist.dna(dna_seqs, model = "raw",
                                  pairwise.deletion = TRUE))

  # Compute within-cluster sum of squares for k=2..max_k
  cat(sprintf("Testing k = 2 to %d clusters...\n", max_k))
  wss <- sapply(2:max_k, function(k) {
    km <- kmeans(dist_mat, centers = k, nstart = 10, iter.max = 200)
    km$tot.withinss
  })

  # Plot elbow curve
  plot(2:max_k, wss,
       type  = "b",
       pch   = 19,
       xlab  = "Number of clusters (k)",
       ylab  = "Total within-cluster sum of squares",
       main  = "Elbow method for optimal k\n(look for the 'elbow' point)")

  # Simple elbow detection
  diffs       <- diff(wss)
  diffs2      <- diff(diffs)
  elbow_k     <- which.max(diffs2) + 2
  cat(sprintf("\nSuggested optimal k: %d\n", elbow_k))
  abline(v = elbow_k, col = "red", lty = 2)

  return(elbow_k)
}

# Run elbow analysis:
# optimal_k <- find_optimal_clusters("C:/path/to/alignment.fasta")
