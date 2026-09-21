library(data.table)
library(stringr)
library(Biostrings)

unexamined_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/"
unexamined_fasta <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/unexamined_fasta/"

## ---- 1. Load pair-status classification (reload if not already in session) --

combined_pair_status <- fread(file.path(unexamined_dir, "matched_but_rare_pair_status.csv"))

cat("Full population sizes:\n")
print(table(combined_pair_status$pair_status))

## ---- 2. Random, stratified sampling ------------------------------------------
## Two separate populations tested here:
##  - "double-primed" (matched_pair, chimera): both ends anchored by a primer
##  - "single-primed" (F_only, R_only): only one end anchored -- much larger
##    population (~1.47M reads), more unconstrained internal sequence exposed,
##    but promiscuity conclusions here are slightly less certain since there's
##    no second-primer cross-confirmation of specificity like matched_pair has

N_PER_CATEGORY <- 150
set.seed(123)

all_categories <- c("matched_pair", "chimera", "F_only", "R_only")
selected <- combined_pair_status[pair_status %in% all_categories]

sampled <- selected[, .SD[sample(.N, min(.N, N_PER_CATEGORY))], by = pair_status]

cat("\nSampled reads per category:\n")
print(table(sampled$pair_status))
cat("\nSample spread across source samples:\n")
print(table(sampled$sample, sampled$pair_status))

## ---- 3. Pull actual sequences for the sampled reads --------------------------

all_unexamined_fastas <- list.files(unexamined_fasta, pattern = "_unexamined\\.fasta$", full.names = TRUE)
seqs_by_sample <- lapply(
  setNames(all_unexamined_fastas, str_remove(basename(all_unexamined_fastas), "_unexamined\\.fasta$")),
  readDNAStringSet
)

pull_seq <- function(sample, qseqid) {
  s <- seqs_by_sample[[sample]]
  if (is.null(s)) return(NULL)
  hit <- s[str_remove(names(s), ";size=\\d+$") == qseqid]
  if (length(hit) != 1) return(NULL)
  hit
}

seq_list <- Map(pull_seq, sampled$sample, sampled$qseqid)
valid <- sapply(seq_list, function(x) is(x, "XStringSet") && length(x) == 1)
cat("\nValid sequences pulled:", sum(valid), "out of", length(seq_list), "\n")

sampled_seqs <- do.call(c, unname(seq_list[valid]))

## keep a lookup table of which category/target-assignment each sequence had,
## since headers alone won't carry that after writing to fasta
sampled_meta <- sampled[valid, .(qseqid, sample, pair_status, abundance, f_targets, r_targets)]
fwrite(sampled_meta, file.path("double_primed_sample_metadata.csv"))

writeXStringSet(sampled_seqs, file.path("double_primed_sample.fasta"))
cat("\nSaved", length(sampled_seqs), "sequences to double_primed_sample.fasta\n")
cat("(filename kept for compatibility -- this now includes matched_pair, chimera,\n")
cat(" F_only, and R_only categories, not just the double-primed ones)\n")
cat("Next: BLAST this file (whole-read, nucleotide) against P450_final_targets_IDT_db\n")
cat("  blastn -query unexamined_check/double_primed_sample.fasta -db P450_final_targets_IDT_db -outfmt '6 qseqid sseqid pident length mismatch evalue bitscore qlen' -evalue 1e-10 -max_target_seqs 3 -out unexamined_check/double_primed_vs_targets.tsv\n")
