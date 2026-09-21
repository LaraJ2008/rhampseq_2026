library(data.table)
library(stringr)
library(Biostrings)

## =============================================================================
## Isolate the "no primer match at all" reads -- these were only ever counted
## by subtraction before (total_rare - everything classified), never saved as
## their own file. This extracts the real sequences.
## =============================================================================

unexamined_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/unexamined_fasta"  # <-- confirm this is correct

pair_status <- fread(file.path("matched_but_rare_pair_status.csv"))

## every (sample, qseqid) combination that DID match at least one primer
matched_keys <- pair_status[, paste(sample, qseqid)]

all_fastas <- list.files(unexamined_dir, pattern = "_unexamined\\.fasta$", full.names = TRUE)
sample_names <- str_remove(basename(all_fastas), "_unexamined\\.fasta$")

no_match_list <- list()

for (i in seq_along(all_fastas)) {
  s <- sample_names[i]
  seqs <- readDNAStringSet(all_fastas[i])
  ids <- str_remove(names(seqs), ";size=\\d+$")
  abundance <- as.integer(str_extract(names(seqs), "(?<=size=)\\d+"))
  seq_strings <- as.character(seqs)   # actual DNA content, not just the ID

  keys_this_sample <- paste(s, ids)
  is_unmatched <- !(keys_this_sample %in% matched_keys)

  if (sum(is_unmatched) > 0) {
    no_match_list[[s]] <- data.table(
      sample = s,
      qseqid = ids[is_unmatched],
      abundance = abundance[is_unmatched],
      sequence = seq_strings[is_unmatched]   # <-- the actual sequence, for
                                               # real cross-sample dereplication
    )
  }
}

no_match_dt <- rbindlist(no_match_list)
cat("Total 'no primer match' READS found:", nrow(no_match_dt), "\n")
cat("Distinct READ IDs:", uniqueN(no_match_dt$qseqid),
    "(always equals total -- read IDs are unique per physical cluster",
    "by definition, this does NOT measure sequence redundancy)\n\n")

## ---- The check that actually matters: true sequence-level redundancy -------
## Group by actual DNA content, summed across ALL samples -- reveals whether
## the same real sequence recurs (possibly across different samples, since
## dereplication was only ever done WITHIN each sample separately)

seq_level <- no_match_dt[, .(
  n_occurrences = .N,
  total_abundance = sum(abundance),
  samples_seen_in = uniqueN(sample)
), by = sequence]
setorder(seq_level, -total_abundance)

cat("=== TRUE distinct sequence variants:", nrow(seq_level), "===\n")
cat("(this is the number that matters -- compare to", nrow(no_match_dt),
    "total reads)\n\n")

cat("Top 20 most-recurring sequences (by total abundance across all samples):\n")
print(head(seq_level[, .(n_occurrences, total_abundance, samples_seen_in)], 20))

cat("\nSequences appearing in more than one sample:",
    sum(seq_level$samples_seen_in > 1), "\n")

fwrite(seq_level[, .(n_occurrences, total_abundance, samples_seen_in)],
       file.path(unexamined_dir, "no_primer_match_sequence_level_summary.csv"))
fwrite(no_match_dt[, .(sample, qseqid, abundance)],
       file.path(unexamined_dir, "no_primer_match_all.csv"))

cat("\nReads per sample (top 10 by count):\n")
print(head(no_match_dt[, .N, by = sample][order(-N)], 10))
cat("\n")

## ---- Stratified sample for downstream analysis: most-recurring sequences ---
## PLUS a random spread, built directly from already-loaded content (no need
## to re-read fasta files -- avoids repeating the path-matching issues from
## earlier tonight)

set.seed(789)

## one representative read per unique sequence (first occurrence), joined
## back to its recurrence stats
first_occurrence <- no_match_dt[, .SD[1], by = sequence]
first_occurrence <- merge(first_occurrence, seq_level, by = "sequence")

N_TOP_RECURRING <- 20
N_RANDOM <- 71

top_recurring <- first_occurrence[order(-total_abundance)][1:min(N_TOP_RECURRING, .N)]
remaining <- first_occurrence[!qseqid %in% top_recurring$qseqid]
random_draw <- remaining[sample(.N, min(N_RANDOM, .N))]

sampled <- rbind(top_recurring, random_draw)

cat("Sampled", nrow(sampled), "sequences:",
    nrow(top_recurring), "most-recurring +", nrow(random_draw), "random\n")
cat("\nSampled reads per source sample:\n")
print(table(sampled$sample))

## build the fasta directly from the sequence strings already in hand
no_match_seqs <- DNAStringSet(setNames(sampled$sequence, sampled$qseqid))
writeXStringSet(no_match_seqs, file.path(unexamined_dir, "no_primer_match_sample.fasta"))

fwrite(sampled[, .(sample, qseqid, abundance, total_abundance, samples_seen_in)],
       file.path(unexamined_dir, "no_primer_match_sample_metadata.csv"))

cat("\nSaved", length(no_match_seqs), "sequences to",
    file.path(unexamined_dir, "no_primer_match_sample.fasta"), "\n")
cat("(includes the", N_TOP_RECURRING, "most-recurring sequences by design --",
    "not a purely random sample -- see no_primer_match_sample_metadata.csv",
    "for which reads are which)\n")
