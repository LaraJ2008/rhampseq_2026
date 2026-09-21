library(data.table)
library(stringr)
library(Biostrings)

## =============================================================================
## PART 1: run this first, BEFORE the bash script
## =============================================================================
## Combines three sources into one fasta + metadata table:
##   - 37 identified targets (full labels, colored by taxonomic order)
##   - 91 confirmed-unresolved chimera reads (no label, grey points)
##   - 40 single-primed reads, 20 F_only + 20 R_only (no label, distinct color)

best_hit_lineage <- fread("best_hit_per_target_with_lineage.csv")

## ---- 1. Target sequences (unchanged from before) -----------------------------

targets <- readDNAStringSet("P450_final_targets_IDT.fasta")
target_seqs <- targets[names(targets) %in% best_hit_lineage$qseqid]

target_meta <- data.table(
  label = names(target_seqs),
  group = "Identified target",
  order = best_hit_lineage$order[match(names(target_seqs), best_hit_lineage$qseqid)]
)

## ---- 2. The 91 confirmed-unresolved chimera reads ----------------------------

chimera_seqs <- readDNAStringSet("D:/NSF_FAST_poc_extracted/unexplained_chimera_for_uniprot.fasta")
names(chimera_seqs) <- str_remove(names(chimera_seqs), ";size=\\d+$")

chimera_meta <- data.table(
  label = names(chimera_seqs),
  group = "Unresolved chimera",
  order = NA_character_
)

## ---- 3. Representative F_only / R_only reads (20 each) -----------------------

double_primed_meta <- fread("double_primed_sample_metadata.csv")
double_primed_seqs <- readDNAStringSet("double_primed_sample.fasta")
names(double_primed_seqs) <- str_remove(names(double_primed_seqs), ";size=\\d+$")

set.seed(456)
f_only_ids <- double_primed_meta[pair_status == "F_only"][sample(.N, min(20, .N)), qseqid]
r_only_ids <- double_primed_meta[pair_status == "R_only"][sample(.N, min(20, .N)), qseqid]

single_primed_seqs <- double_primed_seqs[names(double_primed_seqs) %in% c(f_only_ids, r_only_ids)]
single_primed_meta <- data.table(
  label = names(single_primed_seqs),
  group = "Single-primed (F/R only)",
  order = NA_character_
)

## ---- 4. Combine everything, write out ----------------------------------------

all_seqs <- c(target_seqs, chimera_seqs, single_primed_seqs)
all_meta <- rbindlist(list(target_meta, chimera_meta, single_primed_meta))

## Sanitize names for Newick compatibility: colons (":") are the branch-length
## separator in Newick format, and FastTree truncates any name at the first
## colon -- causing ALL MiSeq-derived read IDs (which start "M03520:...") to
## collapse to the identical truncated name. Also strip other Newick-special
## characters (parens, commas, semicolons) just in case, and keep the
## metadata label column in sync so tip-matching still works downstream.

sanitize_name <- function(x) str_replace_all(x, "[:,;()]", "_")

names(all_seqs) <- sanitize_name(names(all_seqs))
all_meta[, label := sanitize_name(label)]

stopifnot(!anyDuplicated(names(all_seqs)))  # confirm names are still unique
                                              # after sanitizing, fail loudly
                                              # here rather than silently later

cat("Total sequences for tree:\n")
print(table(all_meta$group))

writeXStringSet(all_seqs, "tree_input_all_groups.fasta")
fwrite(all_meta, "tree_input_metadata.csv")

writeLines(names(all_seqs), "target_ids_for_tree.txt")
cat("\nWrote tree_input_all_groups.fasta and tree_input_metadata.csv\n")
cat("Now run build_phylo_tree.sh in bash, then come back and run Part 2 below.\n")

## =============================================================================
## PART 2: run this AFTER bash finishes (tree_targets.nwk exists)
## =============================================================================

## Uncomment and run once the tree file exists:

if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}
library(ggtree)
library(treeio)
library(ggplot2)
library(stringr)

tree <- read.tree("tree_targets.nwk")

tip_data <- fread("tree_input_metadata.csv")
tip_data[is.na(order) & group == "Identified target", order := "Unclassified"]

## color scheme: targets colored by order (Caulobacterales highlighted),
# chimeras and single-primed reads get their own flat, muted colors
orders <- unique(tip_data[group == "Identified target", order])
order_colors <- setNames(scales::hue_pal()(length(orders)), orders)
order_colors["Caulobacterales"] <- "#2166AC"

tip_data[, plot_color := fifelse(group == "Identified target", order, group)]
color_map <- c(order_colors,
               "Unresolved chimera" = "grey60",
               "Single-primed (F/R only)" = "#B2182B")

p <- ggtree(tree, layout = "rectangular") %<+% tip_data +
  ## small unlabeled points for chimera/single-primed reads
  geom_tippoint(data = td_filter(group != "Identified target"),
                aes(color = plot_color), size = 1, alpha = 0.6) +
  ## full labels only for the 37 identified targets
  geom_tiplab(data = td_filter(group == "Identified target"),
              aes(color = plot_color), size = 2.5, align = TRUE, linesize = 0.2) +
  geom_tippoint(data = td_filter(group == "Identified target"),
                aes(color = plot_color), size = 2) +
  scale_color_manual(values = color_map, name = NULL) +
  theme(legend.position = "right") +
  labs(title = "Phylogenetic placement of identified targets vs. chimera and single-primed reads",
       subtitle = "Labeled tips = identified targets (colored by order); grey = unresolved chimeras; red = single-primed reads")

ggsave("target_phylogenetic_tree_with_chimeras.png", p, width = 12, height = 16, dpi = 300)
print(p)
