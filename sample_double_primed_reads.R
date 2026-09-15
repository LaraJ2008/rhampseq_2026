library(data.table)
library(stringr)
library(Biostrings)

unexamined_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/"
unexamined_fasta <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/unexamined_fasta/"
## ---- 1. Load pair-status classification (reload if not already in session) --

combined_pair_status <- fread(file.path(unexamined_dir, "matched_but_rare_pair_status.csv"))

cat("Full population sizes:\n")
print(table(combined_pair_status$pair_status))

## ---- 2. Random, stratified sampling from BOTH double-primed categories -----
## "double-primed" = both F and R hit something (matched_pair OR chimera),
## as opposed to F_only/R_only (single end only). Sample generously but
## tractably from each; random draw avoids the alphabetical-tiebreak bug
## from the earlier chimera-only sampling.

N_PER_CATEGORY <- 150
set.seed(123)

double_primed <- combined_pair_status[pair_status %in% c("matched_pair", "chimera")]

sampled <- double_primed[, .SD[sample(.N, min(.N, N_PER_CATEGORY))], by = pair_status]

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
fwrite(sampled_meta, file.path(unexamined_dir, "double_primed_sample_meta_lookup.csv"))

writeXStringSet(sampled_seqs, file.path(unexamined_dir, "double_primed_sample.fasta"))
cat("\nSaved", length(sampled_seqs), "sequences to double_primed_sample.fasta\n")
cat("Next: BLAST this file (whole-read, nucleotide) against P450_final_targets_IDT_db\n")
cat("  blastn -query unexamined_check/double_primed_sample.fasta -db P450_final_targets_IDT_db -outfmt '6 qseqid sseqid pident length mismatch evalue bitscore qlen' -evalue 1e-10 -max_target_seqs 3 -out unexamined_check/double_primed_vs_targets.tsv\n")

#double_primed_sample.fasta in unexamined_dir