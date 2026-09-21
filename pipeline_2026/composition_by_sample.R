library(data.table)
library(ggplot2)
library(stringr)

## =============================================================================
## Per-sample version of the read composition figure
## =============================================================================

## ---- 1. Per-sample abundant matched/unmatched counts -------------------------

off_target <- fread("C_off_target_check/final_off_target_rates.csv")
off_target <- off_target[, .(sample, total_reads, reads_examined,
                              abundant_matched = reads_matched,
                              abundant_unmatched = reads_unmatched)]

## ---- 2. Per-sample rare-read breakdown ---------------------------------------

pair_status <- fread("matched_but_rare_pair_status.csv")

rare_by_sample <- dcast(pair_status[, .N, by = .(sample, pair_status)],
                         sample ~ pair_status, value.var = "N", fill = 0)

## ensure all four category columns exist even if a sample had zero of one
for (cat_col in c("matched_pair", "chimera", "F_only", "R_only")) {
  if (!cat_col %in% names(rare_by_sample)) rare_by_sample[[cat_col]] <- 0
}

## ---- 3. Combine and compute the remaining "no primer match at all" ---------

composition_by_sample <- merge(off_target, rare_by_sample, by = "sample", all.x = TRUE)
composition_by_sample[is.na(matched_pair), matched_pair := 0]
composition_by_sample[is.na(chimera), chimera := 0]
composition_by_sample[is.na(F_only), F_only := 0]
composition_by_sample[is.na(R_only), R_only := 0]

composition_by_sample[, total_rare := total_reads - reads_examined]
composition_by_sample[, rare_no_match := total_rare -
                         (matched_pair + chimera + F_only + R_only)]

## ---- 3b. Apply a minimum total-reads cutoff -----------------------------------
## Samples with very low total read depth (e.g. failed library prep, a
## near-empty well) produce misleading proportions -- a single read can swing
## the percentage dramatically -- and add little real signal to the raw-count
## view either. Adjust MIN_TOTAL_READS as needed.
off_target <- fread("C_off_target_check/final_off_target_rates.csv")
print(off_target[order(total_reads), .(sample, total_reads)])


MIN_TOTAL_READS <- 50000

excluded <- composition_by_sample[total_reads < MIN_TOTAL_READS]
composition_by_sample <- composition_by_sample[total_reads >= MIN_TOTAL_READS]

cat("=== Sample filtering (MIN_TOTAL_READS =", format(MIN_TOTAL_READS, big.mark = ","), ") ===\n")
cat("Samples kept:", nrow(composition_by_sample), "\n")
cat("Samples excluded (likely failed):", nrow(excluded), "\n")
if (nrow(excluded) > 0) {
  print(excluded[order(total_reads), .(sample, total_reads)])
}
cat("\n")

## ---- 4. Reshape to long format for plotting ----------------------------------

long_comp <- melt(composition_by_sample,
                   id.vars = "sample",
                   measure.vars = c("abundant_matched", "abundant_unmatched",
                                     "matched_pair", "chimera", "F_only", "R_only",
                                     "rare_no_match"),
                   variable.name = "group", value.name = "count")

## human-readable labels, same categories/colors as the aggregate figure
group_labels <- c(
  abundant_matched   = "Abundant (matches primer)",
  abundant_unmatched = "Abundant (no primer match)",
  matched_pair       = "Rare: matched pair",
  chimera            = "Rare: chimera",
  F_only             = "Rare: single-primed (F only)",
  R_only             = "Rare: single-primed (R only)",
  rare_no_match      = "Rare: no primer match"
)
long_comp[, group := factor(group_labels[as.character(group)],
                             levels = group_labels)]

## ---- 5. Parse sample names into treatment series + dose + replicate --------
## e.g. "as10-1" -> treatment "as", dose 10, replicate 1

long_comp[, treatment := str_extract(sample, "^[a-z]+")]
long_comp[, dose := as.integer(str_extract(sample, "(?<=[a-z])\\d+"))]
long_comp[, replicate := str_extract(sample, "(?<=-)\\d+$")]

long_comp[, sample := factor(sample, levels = unique(sample[order(treatment, dose, replicate)]))]

fwrite(long_comp, "read_composition_by_sample_50000.csv")

## ---- 6. Plot: faceted by treatment series, ordered by dose -------------------

color_map <- c(
  "Abundant (matches primer)"          = "#2166AC",
 # "Abundant (no primer match)"         = "#B2182B",
  "Rare: matched pair"                 = "#66C2A5",
  "Rare: chimera"                      = "#FC8D62",
  "Rare: single-primed (F only)"       = "#8DA0CB",
  "Rare: single-primed (R only)"       = "#A6D854",
  "Rare: no primer match"              = "#E5C494"
)

p <- ggplot(long_comp, aes(x = sample, y = count, fill = group)) +
  geom_col(width = 0.75) +
  facet_wrap(~ treatment, scales = "free_x", ncol = 1,
             labeller = as_labeller(c(as = "AS treatment series",
                                       sd = "SD treatment series"))) +
  scale_fill_manual(values = color_map, name = NULL) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(title = "Read composition by sample",
       subtitle = paste0("Ordered by dose within each treatment series  |  ",
                          nrow(excluded), " low-depth sample(s) excluded (< ",
                          format(MIN_TOTAL_READS, big.mark = ","), " total reads)"),
       x = NULL, y = "Number of reads") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))

ggsave("read_composition_by_sample_50000.png", p, width = 13, height = 10, dpi = 300)
cat("Saved to read_composition_by_sample.png\n")
print(p)

## ---- 7. Also produce a PROPORTION (100%-stacked) version -------------------
## Raw counts make small samples (e.g. as30-1, ~330 total reads) invisible next
## to large ones (~600k). This version normalizes each sample to 100% so
## COMPOSITION is comparable regardless of sequencing depth -- useful as a
## complementary view, not a replacement.

p_prop <- ggplot(long_comp, aes(x = sample, y = count, fill = group)) +
  geom_col(width = 0.75, position = "fill") +
  facet_wrap(~ treatment, scales = "free_x", ncol = 1,
             labeller = as_labeller(c(as = "AS treatment series",
                                       sd = "SD treatment series"))) +
  scale_fill_manual(values = color_map, name = NULL) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(title = "Read composition by sample (proportion)",
       subtitle = paste0("Normalized to 100%  |  ", nrow(excluded),
                          " low-depth sample(s) excluded (< ",
                          format(MIN_TOTAL_READS, big.mark = ","), " total reads)"),
       x = NULL, y = "Proportion of reads") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))

ggsave("read_composition_by_sample_proportion_5000.png", p_prop, width = 13, height = 10, dpi = 300)
cat("Saved to read_composition_by_sample_proportion.png\n")
print(p_prop)
