## Run this AFTER reformat_rhampseq_primers_v2.R has produced `primers` in your
## R session (or reload rhampseq_primers_clean.csv if starting fresh).

library(data.table)
library(Biostrings)

primers <- fread("rhampseq_primers_clean.csv")  # skip if `primers` is already loaded

## Build one FASTA entry per primer (F and R separately), ID = target_id + orientation
all_primer_seqs <- c(
  setNames(primers$primer_F_clean, paste0(primers$target_id, "_F")),
  setNames(primers$primer_R_clean, paste0(primers$target_id, "_R"))
)

writeXStringSet(DNAStringSet(all_primer_seqs), filepath = "all_primers.fasta")

cat("Wrote", length(all_primer_seqs), "primer sequences to all_primers.fasta\n")
