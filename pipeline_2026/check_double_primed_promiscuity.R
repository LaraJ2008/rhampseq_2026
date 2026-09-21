library(data.table)
library(stringr)

unexamined_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/"

## ---- 1. Load everything needed ------------------------------------------------

sampled_meta <- fread(file.path(unexamined_dir, "double_primed_sample_meta_lookup.csv"))

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

results[, classification := fcase(
  is.na(read_best_source_id), "No whole-read match to any target",
  matches_expected == TRUE & read_pident >= 90, "Matches own primer's target (expected)",
  matches_expected == TRUE & read_pident < 90, "Matches own primer's target but LOW identity (divergent variant)",
  matches_expected == FALSE, "Matches a DIFFERENT target than its own primer(s) -- promiscuity signal",
  default = "Unclassified"
)]

fwrite(results, file.path("double_primed_promiscuity_check.csv"))

cat("=== Promiscuity check results ===\n\n")
cat("By pair_status category:\n")
print(table(results$pair_status, results$classification))

cat("\n=== Reads showing a promiscuity signal ===\n")
print(results[classification == "Matches a DIFFERENT target than its own primer(s) -- promiscuity signal",
              .(qseqid, sample, pair_status, expected_sources, read_best_source_id,
                read_pident, read_length)])

cat("\n=== Overall summary ===\n")
print(table(results$classification))

promiscuity_reads <- results[classification == "Matches a DIFFERENT target than its own primer(s) -- promiscuity signal"]
print(promiscuity_reads)

fwrite(promiscuity_reads, file.path("promiscuous_primer_pairs.csv"))
