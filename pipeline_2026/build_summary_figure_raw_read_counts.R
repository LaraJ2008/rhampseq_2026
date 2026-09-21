library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)

## =============================================================================
## FIGURE PART A: Overall read composition
## =============================================================================
## Pulls totals directly from files already generated earlier -- no need to
## recompute anything, just assemble what we already have.

## ---- 1. Abundant (>=5) vs Rare (<5) totals, and the abundant off-target rate

off_target <- fread("C_off_target_check/final_off_target_rates.csv")

total_reads      <- sum(off_target$total_reads, na.rm = TRUE)
total_abundant   <- sum(off_target$reads_examined, na.rm = TRUE)
abundant_matched <- sum(off_target$reads_matched, na.rm = TRUE)
abundant_unmatched <- sum(off_target$reads_unmatched, na.rm = TRUE)
total_rare       <- total_reads - total_abundant

cat("Total reads:", format(total_reads, big.mark = ","), "\n")
cat("Abundant (>=5 copies):", format(total_abundant, big.mark = ","),
    sprintf("(%.1f%%)", 100 * total_abundant / total_reads), "\n")
cat("Rare (<5 copies):", format(total_rare, big.mark = ","),
    sprintf("(%.1f%%)", 100 * total_rare / total_reads), "\n")

## ---- 2. Rare-read breakdown: matched_pair / chimera / F_only / R_only / -----
## no primer match at all

pair_status <- fread("unexamined_check/matched_but_rare_pair_status.csv")
rare_classified <- pair_status[, .N, by = pair_status][, setNames(N, pair_status)]

n_matched_pair <- rare_classified[["matched_pair"]]
n_chimera      <- rare_classified[["chimera"]]
n_f_only       <- rare_classified[["F_only"]]
n_r_only       <- rare_classified[["R_only"]]
n_rare_any_primer <- n_matched_pair + n_chimera + n_f_only + n_r_only
n_rare_no_match    <- total_rare - n_rare_any_primer

cat("\nRare-read breakdown:\n")
cat("  matched_pair:", format(n_matched_pair, big.mark = ","), "\n")
cat("  chimera:", format(n_chimera, big.mark = ","), "\n")
cat("  F_only:", format(n_f_only, big.mark = ","), "\n")
cat("  R_only:", format(n_r_only, big.mark = ","), "\n")
cat("  no primer match at all:", format(n_rare_no_match, big.mark = ","), "\n")

## ---- 3. Assemble one composition table for plotting -------------------------

composition <- data.table(
  group = c("Abundant\n(matches primer)", "Abundant\n(no primer match)",
            "Rare: matched pair", "Rare: chimera",
            "Rare: single-primed (F only)", "Rare: single-primed (R only)",
            "Rare: no primer match"),
  category = c("Abundant (\u22655 copies)", "Abundant (\u22655 copies)",
               "Rare (<5 copies)", "Rare (<5 copies)",
               "Rare (<5 copies)", "Rare (<5 copies)", "Rare (<5 copies)"),
  count = c(abundant_matched, abundant_unmatched,
            n_matched_pair, n_chimera, n_f_only, n_r_only, n_rare_no_match)
)
composition[, pct := 100 * count / total_reads]
composition[, group := factor(group, levels = rev(group))]

fwrite(composition, "read_composition_summary.csv")

## ---- 4. Plot: stacked composition bar -----------------------------------------

p_composition <- ggplot(composition, aes(x = category, y = count, fill = group)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = ifelse(pct > 2, scales::comma(count), "")),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = c(
    "Abundant\n(matches primer)"          = "#2166AC",
    "Abundant\n(no primer match)"         = "#B2182B",
    "Rare: matched pair"                  = "#66C2A5",
    "Rare: chimera"                       = "#FC8D62",
    "Rare: single-primed (F only)"        = "#8DA0CB",
    "Rare: single-primed (R only)"        = "#A6D854",
    "Rare: no primer match"               = "#E5C494"
  )) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(title = "Read composition across all samples",
       x = NULL, y = "Number of reads", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 10))

## =============================================================================
## FIGURE PART B: Taxonomic order of the highest-abundance reads
## =============================================================================
## Uses order-level classification (from the lineage lookup) rather than
## genus -- far fewer NAs/singletons, and reveals the dominant Sphingomonadales
## signal much more clearly than the scattered genus-level view did.

best_hit_lineage <- fread("best_hit_per_target_with_lineage.csv")

order_counts <- best_hit_lineage[, .N, by = order][order(-N)]
order_counts[is.na(order), order := "Unclassified"]

## flag Caulobacterales for a highlighted color, same as before
order_counts[, highlight := order == "Caulobacterales"]

p_identity <- ggplot(order_counts, aes(x = reorder(order, N), y = N, fill = highlight)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey70"),
                     guide = "none") +
  scale_y_continuous(limits = c(0, max(order_counts$N) * 1.15)) +
  coord_flip() +
  labs(title = "Taxonomic order of best-hit matches for highest-abundance targets",
       subtitle = "Blue = Caulobacterales  |  Sphingomonadales dominates overall",
       x = NULL, y = "Number of targets (n = 37)") +
  theme_minimal(base_size = 11)

## =============================================================================
## Combine and save
## =============================================================================

combined <- p_composition / p_identity +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(title = "P450 rhAmpSeq panel: read-level summary",
                   theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave("read_characteristics_summary_figure_raw_counts.png", combined,
       width = 10, height = 11, dpi = 300)

cat("\nSaved figure to read_characteristics_summary_figure.png\n")
print(combined)
