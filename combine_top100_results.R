library(data.table)
library(stringr)
library(Biostrings)

## ---- 0. CONFIG ---------------------------------------------------------------
out_dir <- "dada2_bypass_check"
DOESN'T TOTALLY WORK?

## ---- 1. Get abundance per read from the FASTA headers -----------------------
## vsearch --sizeout headers look like: >abc123;size=4567
top100_seqs <- readDNAStringSet(file.path(out_dir, "top100_abundant.fasta"))
header_dt <- data.table(
  qseqid    = str_remove(names(top100_seqs), ";size=\\d+$"),
  abundance = as.integer(str_extract(names(top100_seqs), "(?<=size=)\\d+"))
)
setorder(header_dt, -abundance)

cat("Top 100 reads loaded, abundance range:",
    min(header_dt$abundance), "-", max(header_dt$abundance), "\n\n")

## ---- 2. Best target match per read --------------------------------------------

target_hits <- fread(file.path(out_dir, "top100_vs_all_targets.tsv"), header = FALSE,
                      col.names = c("qseqid","sseqid","pident","length",
                                    "mismatch","evalue","bitscore","qlen"))
target_hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]

best_target <- target_hits[order(-bitscore)][, .SD[1], by = qseqid]
setnames(best_target, c("sseqid","pident","length","bitscore"),
         c("best_target_source_id","target_pident","target_length","target_bitscore"))

## ---- 3. Genuine primer matches per read ---------------------------------------
## Same filtering logic as the Caulobacter primer check earlier: require
## near-full-length, high-identity matches to distinguish real primer binding
## from short coincidental hits (word_size 7 finds plenty of noise).

primer_hits <- fread(file.path(out_dir, "top100_vs_all_primers.tsv"), header = FALSE,
                      col.names = c("qseqid","sseqid","pident","length",
                                    "mismatch","qstart","qend","sstart","send",
                                    "evalue","bitscore"))
primer_hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]

primer_hits_filtered <- primer_hits[length >= 18 & pident >= 90]

cat("Primer hits before filtering:", nrow(primer_hits), "\n")
cat("Primer hits after filtering (length>=18, pident>=90):", nrow(primer_hits_filtered), "\n\n")

primer_hits_filtered[, primer_target_id := str_remove(sseqid, "_[FR]$")]
primer_hits_filtered[, primer_orientation := str_extract(sseqid, "[FR]$")]

## collapse to one row per read: list of matched primers + their target_ids
primer_summary <- primer_hits_filtered[, .(
  matched_primers    = paste(unique(sseqid), collapse = "; "),
  matched_primer_targets = paste(unique(primer_target_id), collapse = "; "),
  n_primers_matched  = uniqueN(sseqid)
), by = qseqid]

## ---- 4. Combine everything into one results table -----------------------------

results <- merge(header_dt, best_target[, .(qseqid, best_target_source_id,
                                             target_pident, target_length, target_bitscore)],
                  by = "qseqid", all.x = TRUE)
results <- merge(results, primer_summary, by = "qseqid", all.x = TRUE)

results[is.na(n_primers_matched), n_primers_matched := 0]
results[is.na(matched_primers), matched_primers := ""]

setorder(results, -abundance)

fwrite(results, file.path(out_dir, "top100_combined_results.csv"))

cat("=== Top 20 most abundant reads: target + primer match summary ===\n")
print(results[1:20, .(qseqid, abundance, best_target_source_id,
                       target_pident, matched_primer_targets, n_primers_matched)])

## ---- 5. The "mapped read specific primers" list -- for your figure -----------
## Only the primers that GENUINELY matched at least one of the top 100 reads.

mapped_primer_ids <- unique(primer_hits_filtered$sseqid)
mapped_primer_target_ids <- unique(primer_hits_filtered$primer_target_id)

cat("\n=== Primers that genuinely matched at least one top-100 read ===\n")
cat("Individual primers (F/R):", length(mapped_primer_ids), "\n")
cat("Distinct target_ids represented:", length(mapped_primer_target_ids), "\n\n")
print(mapped_primer_target_ids)

writeLines(mapped_primer_ids, file.path(out_dir, "mapped_primers_for_figure.txt"))
cat("\nSaved to", file.path(out_dir, "mapped_primers_for_figure.txt"), "\n")

## ---- 6. Flag any read where the matched primer's target DISAGREES with the ----
## best BLAST target match -- worth a look, could indicate cross-reactivity
## or a primer binding somewhere other than its intended target.

results[, target_primer_mismatch := !is.na(best_target_source_id) &
           n_primers_matched > 0 &
           !str_detect(matched_primer_targets, fixed(str_extract(best_target_source_id, "^[^_]+")))]

cat("\n=== Reads where matched primer's target differs from best BLAST target match ===\n")
print(results[target_primer_mismatch == TRUE,
              .(qseqid, abundance, best_target_source_id, matched_primer_targets)])
