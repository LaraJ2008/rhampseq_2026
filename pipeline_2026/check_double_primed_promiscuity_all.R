library(data.table)
library(stringr)

unexamined_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/"

## ---- 1. Load everything needed ------------------------------------------------

sampled_meta <- fread(file.path(unexamined_dir, "double_primed_sample_metadata.csv"))

target_hits <- fread(file.path(unexamined_dir, "double_primed_vs_targets.tsv"), header = FALSE,
                      col.names = c("qseqid","sseqid","pident","length",
                                    "mismatch","evalue","bitscore","qlen"))
target_hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]

# best whole-read match per read
best_whole_read_match <- target_hits[order(-bitscore)][, .SD[1], by = qseqid]
setnames(best_whole_read_match, c("sseqid","pident","length","bitscore"),
         c("read_best_source_id","read_pident","read_length","read_bitscore"))

# the target_id -> true source lookup built earlier in this investigation
source_map <- fread("target_id_to_source_map.csv")

## ---- 2. For each sampled read, get the source ID(s) its OWN primer(s) -------
## were designed to amplify -- matched_pair reads have one target (agree),
## chimera reads have two candidate targets (f_targets, r_targets)

sampled_meta[, f_target_id := str_split(f_targets, "; ")]
sampled_meta[, r_target_id := str_split(r_targets, "; ")]

get_expected_sources <- function(f_ids, r_ids) {
  all_ids <- unique(c(f_ids, r_ids))
  sources <- source_map[target_id %in% all_ids, true_source_id]
  paste(unique(sources), collapse = "; ")
}

sampled_meta[, expected_sources := mapply(get_expected_sources, f_target_id, r_target_id)]

## ---- 3. Join with the actual whole-read BLAST result and compare -----------

results <- merge(sampled_meta, best_whole_read_match, by = "qseqid", all.x = TRUE)

results[, matches_expected := !is.na(read_best_source_id) &
          str_detect(expected_sources, fixed(read_best_source_id))]

## coverage: what fraction of the read's own length does the alignment span?
## a real, legitimate match should cover most of the read, not just a short
## coincidental stretch
results[, read_coverage := read_length / qlen]

## Thresholds for calling something a legitimate promiscuity signal (adjust
## as needed): high identity AND good length/coverage -- matching the
## confidence level of the RH.E04EE9/RH.8B25A0 hits already confirmed
## (97-100% identity, ~150bp, essentially full-length)
MIN_PIDENT <- 90
MIN_COVERAGE <- 0.60

results[, classification := fcase(
  is.na(read_best_source_id), "No whole-read match to any target",

  matches_expected == TRUE & read_pident >= 90, "Matches own primer's target (expected)",
  matches_expected == TRUE & read_pident < 90, "Matches own primer's target but LOW identity (divergent variant)",

  matches_expected == FALSE & read_pident >= MIN_PIDENT & read_coverage >= MIN_COVERAGE,
    "Matches a DIFFERENT target -- LEGITIMATE promiscuity signal (high identity + coverage)",

  matches_expected == FALSE,
    "Matches a DIFFERENT target but WEAK/short match -- likely coincidental, not legitimate",

  default = "Unclassified"
)]

fwrite(results, file.path(unexamined_dir, "double_primed_promiscuity_check_full.csv"))

cat("=== Promiscuity check results (identity/coverage-filtered) ===\n\n")
cat("By pair_status category:\n")
print(table(results$pair_status, results$classification))

cat("\n=== Reads showing a LEGITIMATE promiscuity signal (high identity + coverage) ===\n")
print(results[classification == "Matches a DIFFERENT target -- LEGITIMATE promiscuity signal (high identity + coverage)",
              .(qseqid, sample, pair_status, expected_sources, read_best_source_id,
                read_pident, read_length, read_coverage)])

cat("\n=== Weak/coincidental different-target hits (NOT counted as legitimate) ===\n")
print(results[classification == "Matches a DIFFERENT target but WEAK/short match -- likely coincidental, not legitimate",
              .(qseqid, sample, pair_status, read_pident, read_length, read_coverage)])

cat("\n=== Overall summary ===\n")
print(table(results$classification))
