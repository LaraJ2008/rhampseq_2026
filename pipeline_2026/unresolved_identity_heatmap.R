library(Biostrings)
library(pwalign)
library(data.table)
library(ggplot2)
library(stringr)

seqs <- readDNAStringSet("/mnt/d/NSF_FAST_poc_extracted/unexplained_chimera_for_uniprot.fasta")
names(seqs) <- str_remove(names(seqs), ";size=\\d+$")
n <- length(seqs)
cat("Sequences:", n, "\n")

labels <- paste0("r", seq_len(n))
mat <- matrix(100, nrow = n, ncol = n, dimnames = list(labels, labels))

for (i in 1:(n-1)) {
  ref_repeated <- rep(seqs[i], n - i)
  aln <- pairwiseAlignment(ref_repeated, seqs[(i+1):n])
  vals <- pid(aln)
  mat[i, (i+1):n] <- vals
  mat[(i+1):n, i] <- vals
}

long_dt <- as.data.table(as.table(mat))
setnames(long_dt, c("read_i", "read_j", "pident"))

## hierarchical clustering to reveal any block structure -- order the heatmap
## axes by cluster rather than arbitrary read order, so real groups (if any)
## become visually obvious as blocks along the diagonal
dist_mat <- as.dist(100 - mat)
hc <- hclust(dist_mat, method = "average")
ordered_labels <- labels[hc$order]
long_dt[, read_i := factor(read_i, levels = ordered_labels)]
long_dt[, read_j := factor(read_j, levels = ordered_labels)]

p <- ggplot(long_dt, aes(x = read_i, y = read_j, fill = pident)) +
  geom_tile() +
  scale_fill_gradient2(low = "firebrick", mid = "gold", high = "steelblue",
                        midpoint = 65, name = "% identity") +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Pairwise identity among the 91 unresolved ('non-P450') reads",
       subtitle = "Ordered by hierarchical clustering -- visible blocks indicate related sequences",
       x = NULL, y = NULL)

ggsave("unresolved_reads_identity_heatmap.png", p, width = 8, height = 7, dpi = 300)
print(p)

## report cluster membership at a reasonable threshold (e.g. 70% identity)
clusters <- cutree(hc, h = 30)  # cut at 70% identity
cat("\n=== Cluster sizes at 70% identity threshold ===\n")
print(table(clusters))

cluster_dt <- data.table(label = labels, cluster = clusters)
seq_id_map <- data.table(label = labels, qseqid = names(seqs))
cluster_dt <- merge(cluster_dt, seq_id_map, by = "label")
fwrite(cluster_dt, "unresolved_reads_clusters.csv")

n_singleton_clusters <- sum(table(clusters) == 1)
cat("\nSingleton clusters (no close relative among the 91):", n_singleton_clusters,
    "out of", n, "\n")
cat("If most reads are singletons -> heterogeneous, supports general filtering advice.\n")
cat("If a few large clusters dominate -> consistent gene(s), worth targeted BLAST.\n")
