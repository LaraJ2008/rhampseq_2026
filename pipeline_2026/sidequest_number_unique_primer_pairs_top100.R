# of our top 100 abundance sequences, are there 100 primer pairs being used? Or less?

if (!requireNamespace("pwalign", quietly = TRUE)) {
  BiocManager::install("pwalign")
}
library(pwalign)
library(ggplot2)


results <- fread("dada2_bypass_check/top100_combined_results.csv")

# split the semicolon-separated target lists into one row per read-target pair
target_long <- results[, .(target_id = str_trim(str_split(matched_primer_targets, ";")[[1]])),
                       by = qseqid]

cat("Total (read, target) pairings:", nrow(target_long), "\n")
cat("Distinct target_ids used across the top 100:", uniqueN(target_long$target_id), "\n\n")

# which targets appear more than once -- i.e., captured by MULTIPLE different
# abundant reads, not just one
target_counts <- target_long[, .N, by = target_id][order(-N)]
print(target_counts[N > 1])

#~~~~~~~~~~~~Ok how similar are the sequences in these groupings?~~~~~~~~~~
seqs <- readDNAStringSet("dada2_bypass_check/top100_abundant.fasta")
seq_ids_clean <- str_remove(names(seqs), ";size=\\d+$")

# pull the reads assigned to your most heavily-repeated target
target_reads <- target_long[target_id == "RH.0BD19D1A26A74EBZ0Z", qseqid]
these_seqs <- seqs[seq_ids_clean %in% target_reads]

length(these_seqs)  # sanity check, should be 16
width(these_seqs)   # check they're comparable lengths

# quick pairwise identity check
library(Biostrings)
#aln <- pairwiseAlignment(these_seqs[1], these_seqs[-1])
#pid(aln)  # percent identity of read #1 against each of the other 15


ref_repeated <- rep(these_seqs[1], length(these_seqs) - 1)
aln <- pwalign::pairwiseAlignment(ref_repeated, these_seqs[-1])
pid(aln)

library(data.table)
library(stringr)
library(Biostrings)
library(pwalign)

## Assumes `target_long` and `seqs`/`seq_ids_clean` already exist from your
## session (target_long from splitting matched_primer_targets; seqs from
## top100_abundant.fasta). Recreating them here in case you're starting fresh:

results <- fread("dada2_bypass_check/top100_combined_results.csv")
target_long <- results[, .(target_id = str_trim(str_split(matched_primer_targets, ";")[[1]])),
                       by = qseqid]

seqs <- readDNAStringSet("dada2_bypass_check/top100_abundant.fasta")
seq_ids_clean <- str_remove(names(seqs), ";size=\\d+$")

target_counts <- target_long[, .N, by = target_id][order(-N)]
repeated_targets <- target_counts[N > 1, target_id]

cat("Checking", length(repeated_targets), "targets with multiple abundant reads...\n\n")

## ---- For each repeated target, compute all pairwise identities among its ----
## reads, and summarize the range/mean

identity_summary <- rbindlist(lapply(repeated_targets, function(tid) {
  
  target_reads <- target_long[target_id == tid, qseqid]
  these_seqs <- seqs[seq_ids_clean %in% target_reads]
  
  n <- length(these_seqs)
  if (n < 2) return(NULL)
  
  # all pairwise comparisons among this group's reads
  pident_values <- c()
  for (i in 1:(n-1)) {
    ref_repeated <- rep(these_seqs[i], n - i)
    aln <- pairwiseAlignment(ref_repeated, these_seqs[(i+1):n])
    pident_values <- c(pident_values, pid(aln))
  }
  
  data.table(
    target_id = tid,
    n_reads = n,
    n_comparisons = length(pident_values),
    min_pident = round(min(pident_values), 1),
    mean_pident = round(mean(pident_values), 1),
    max_pident = round(max(pident_values), 1),
    interpretation = if (min(pident_values) >= 90) {
      "High identity throughout -- likely real strain-level SNP diversity, same gene"
    } else if (min(pident_values) >= 75) {
      "Moderate divergence -- worth a closer look, could be allelic OR closely related paralog"
    } else {
      "Substantial divergence -- consistent with primer catching a DIFFERENT related P450"
    }
  )
}))

setorder(identity_summary, min_pident)

print(identity_summary)

fwrite(identity_summary, "dada2_bypass_check/repeated_target_identity_check.csv")

cat("\n=== Summary ===\n")
print(table(identity_summary$interpretation))



library(data.table)
library(stringr)
library(Biostrings)
library(pwalign)
library(ggplot2)

## Assumes target_long, seqs, seq_ids_clean, and repeated_targets already
## exist in your session from the previous script. Recreating here in case
## you're starting fresh in a new session:

results <- fread("dada2_bypass_check/top100_combined_results.csv")
target_long <- results[, .(target_id = str_trim(str_split(matched_primer_targets, ";")[[1]])),
                       by = qseqid]

seqs <- readDNAStringSet("dada2_bypass_check/top100_abundant.fasta")
seq_ids_clean <- str_remove(names(seqs), ";size=\\d+$")

target_counts <- target_long[, .N, by = target_id][order(-N)]
repeated_targets <- target_counts[N > 1, target_id]

## ---- Build a FULL (symmetric, n x n) pairwise identity matrix per target ----

all_pairs_long <- rbindlist(lapply(repeated_targets, function(tid) {
  
  target_reads <- target_long[target_id == tid, qseqid]
  these_seqs <- seqs[seq_ids_clean %in% target_reads]
  n <- length(these_seqs)
  if (n < 2) return(NULL)
  
  # short labels for plotting instead of full read IDs
  labels <- paste0("r", seq_len(n))
  
  # build the full n x n matrix, diagonal = 100 (self-identity)
  mat <- matrix(100, nrow = n, ncol = n, dimnames = list(labels, labels))
  
  for (i in 1:(n-1)) {
    ref_repeated <- rep(these_seqs[i], n - i)
    aln <- pairwiseAlignment(ref_repeated, these_seqs[(i+1):n])
    vals <- pid(aln)
    mat[i, (i+1):n] <- vals
    mat[(i+1):n, i] <- vals   # mirror for symmetry
  }
  
  # melt to long format for ggplot
  long_dt <- as.data.table(as.table(mat))
  setnames(long_dt, c("read_i", "read_j", "pident"))
  long_dt[, target_id := tid]
  long_dt[, n_reads := n]
  long_dt
}))

## order facets by how divergent they are (lowest min identity first), so the
## most "interesting" (potentially promiscuous) targets appear first in the grid
facet_order <- all_pairs_long[, .(min_pident = min(pident[pident < 100])), by = target_id]
setorder(facet_order, min_pident)
all_pairs_long[, target_id := factor(target_id, levels = facet_order$target_id)]

## ---- Plot: faceted heatmap grid ------------------------------------------------

p <- ggplot(all_pairs_long, aes(x = read_i, y = read_j, fill = pident)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "firebrick", mid = "gold", high = "steelblue",
    midpoint = 85, limits = c(min(all_pairs_long$pident), 100),
    name = "% identity"
  ) +
  facet_wrap(~ target_id, scales = "free", ncol = 5) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(title = "Pairwise sequence identity among reads sharing a primer pair",
       subtitle = "Ordered by lowest observed identity (most divergent groups first) -- low identity (red) suggests a primer may be capturing more than one distinct sequence")

ggsave("dada2_bypass_check/repeated_target_identity_heatmap.png",
       plot = p, width = 14, height = 10, dpi = 300)

cat("Saved heatmap to dada2_bypass_check/repeated_target_identity_heatmap.png\n")
print(p)

###WHICH ARE DIVERGENT~~~~~~~~~~~~~~~~~~~~~~~~~~

print(identity_summary[interpretation == "Substantial divergence -- consistent with primer catching a DIFFERENT related P450"])

#~~~~~~~~~~~~~~~~~~~~~look up actual sequences~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NOT YET RUN

divergent_targets <- identity_summary[
  interpretation == "Substantial divergence -- consistent with primer catching a DIFFERENT related P450",
  target_id
]

for (tid in divergent_targets) {
  reads_for_target <- target_long[target_id == tid, qseqid]
  these_seqs <- seqs[seq_ids_clean %in% reads_for_target]
  
  # cluster reads by sequence similarity to identify the distinct sub-groups
  n <- length(these_seqs)
  dist_mat <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    ref_repeated <- rep(these_seqs[i], n - i)
    aln <- pairwiseAlignment(ref_repeated, these_seqs[(i+1):n])
    dist_mat[i, (i+1):n] <- 100 - pid(aln)
    dist_mat[(i+1):n, i] <- dist_mat[i, (i+1):n]
  }
  
  hc <- hclust(as.dist(dist_mat), method = "average")
  clusters <- cutree(hc, h = 15)  # cut at 85% identity threshold
  
  cat("\n", tid, "- clusters found:", length(unique(clusters)), "\n")
  print(table(clusters))
  
  # save one representative sequence per cluster for BLASTing
  for (cl in unique(clusters)) {
    rep_seq <- these_seqs[which(clusters == cl)[1]]
    writeXStringSet(rep_seq,
                    paste0("dada2_bypass_check/", tid, "_cluster", cl, "_representative.fasta"))
  }
}
