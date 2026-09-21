library(data.table)
library(stringr)

out_dir <- "dada2_bypass_check"

## ---- 1. Build the target_id -> true source candidate ID lookup ---------------

primer_source_hits <- fread("all_primers_vs_all_candidates.tsv", header = FALSE,
                             col.names = c("qseqid","sseqid","pident","length",
                                           "mismatch","qstart","qend","sstart",
                                           "send","evalue","bitscore"))

# same filtering threshold used throughout this whole investigation: require
# near-full-length, high-identity matches to trust a primer<->source link
primer_source_filtered <- primer_source_hits[length >= 18 & pident >= 90]

cat("Primer-vs-candidate hits before filtering:", nrow(primer_source_hits), "\n")
cat("After filtering (length>=18, pident>=90):", nrow(primer_source_filtered), "\n\n")

primer_source_filtered[, target_id := str_remove(qseqid, "_[FR]$")]
primer_source_filtered[, orientation := str_extract(qseqid, "[FR]$")]

# best (highest bitscore) source match per individual primer
best_source_per_primer <- primer_source_filtered[order(-bitscore)][, .SD[1], by = qseqid]

# For each target_id, F and R SHOULD point to the same source candidate if the
# primer pair was genuinely designed from one contiguous piece of source
# sequence. Check whether they agree, and flag any that don't as worth a
# second look (could indicate a target assembled from two overlapping reads,
# or a real discrepancy worth investigating).
target_source_map <- dcast(best_source_per_primer, target_id ~ orientation,
                            value.var = "sseqid")

target_source_map[, source_agrees := (F == R)]
n_disagree <- target_source_map[source_agrees == FALSE, .N]
cat("Targets where F and R primers trace to DIFFERENT source candidates:", n_disagree, "\n")
if (n_disagree > 0) {
  cat("(worth a look, but proceeding using F's source as the reference)\n")
  print(target_source_map[source_agrees == FALSE])
}

# use F's matched source as the authoritative source_id per target (arbitrary
# but consistent choice when they disagree; F is the leftmost designed primer)
target_source_map[, true_source_id := F]

fwrite(target_source_map[, .(target_id, true_source_id, source_agrees)],
       "target_id_to_source_map.csv")

cat("\nLookup table built for", nrow(target_source_map), "of your 380 targets.\n")
cat("(", 380 - nrow(target_source_map), "targets have no confident primer-to-source",
    "link at this threshold -- their primers may not BLAST cleanly back to the",
    "511-candidate file, worth noting but not necessarily a problem)\n\n")

## ---- 2. Redo the mismatch check on the top-100 results, correctly this time ---

results <- fread(file.path("top100_combined_results.csv"))

# expand matched_primer_targets (semicolon-separated) into long format so each
# read-target pairing gets its own row for checking
results_long <- results[n_primers_matched > 0]
results_long <- results_long[, .(target_id = str_trim(str_split(matched_primer_targets, ";")[[1]])),
                              by = .(qseqid, abundance, best_target_source_id,
                                     target_pident, n_primers_matched)]

results_long <- merge(results_long, target_source_map[, .(target_id, true_source_id)],
                       by = "target_id", all.x = TRUE)

# NOW a real, ID-scheme-compatible comparison: does the read's own best BLAST
# match (best_target_source_id) equal the TRUE source candidate its matched
# primer was built from?
results_long[, primer_source_agrees_with_read := !is.na(true_source_id) &
                !is.na(best_target_source_id) &
                true_source_id == best_target_source_id]

## collapse back to one row per read: does AT LEAST ONE matched primer's true
## source agree with this read's own best target match?
read_level_check <- results_long[, .(
  any_primer_source_agrees = any(primer_source_agrees_with_read, na.rm = TRUE),
  matched_target_ids       = paste(unique(target_id), collapse = "; "),
  matched_true_sources      = paste(unique(true_source_id), collapse = "; ")
), by = .(qseqid, abundance, best_target_source_id)]

setorder(read_level_check, -abundance)

fwrite(read_level_check, file.path("top100_corrected_mismatch_check.csv"))

cat("=== Corrected check: do matched primers trace back to the SAME source",
    "as the read's own best BLAST match? ===\n")
print(read_level_check[1:20])

cat("\n=== Summary ===\n")
print(table(read_level_check$any_primer_source_agrees, useNA = "ifany"))


problem_reads <- read_level_check[any_primer_source_agrees == FALSE]
print(problem_reads)

# check whether any of these involve one of the 25 known F/R-disagreement targets
disagreement_targets <- target_source_map[source_agrees == FALSE, target_id]
problem_reads[, involves_disagreement_target := sapply(strsplit(matched_target_ids, "; "),
                                                       function(x) any(x %in% disagreement_targets))]
print(problem_reads[, .(qseqid, abundance, best_target_source_id,
                        matched_target_ids, involves_disagreement_target)])
