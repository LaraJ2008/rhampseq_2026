library(data.table)
library(stringr)
library(Biostrings)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stopifnot(dir.exists(data_dir)) is a good tool
## ---- 0. CONFIG -----------------------------------------------------------
## Point this at wherever your unexamined_check files actually live.
data_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/unexamined_fasta"

stopifnot(dir.exists(data_dir))  # fails loudly here instead of silently
                                   # pulling 0 sequences later if the path is wrong

## ---- 1. Load the promiscuity check results and pull the unexplained chimeras -

promiscuity_results <- fread(file.path("double_primed_promiscuity_check.csv"))

unexplained_chimera <- promiscuity_results[
  pair_status == "chimera" & classification == "No whole-read match to any target"
]

cat("Unexplained chimera reads to re-test:", nrow(unexplained_chimera), "\n")

## ---- 2. Pull actual sequences from the per-sample unexamined fasta files -----

all_unexamined_fastas <- list.files(data_dir, pattern = "_unexamined\\.fasta$",
                                     full.names = TRUE)
cat("Per-sample fasta files found:", length(all_unexamined_fastas), "\n")

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

seq_list <- Map(pull_seq, unexplained_chimera$sample, unexplained_chimera$qseqid)
valid <- sapply(seq_list, function(x) is(x, "XStringSet") && length(x) == 1)
cat("Valid sequences pulled:", sum(valid), "out of", length(seq_list), "\n")

unexplained_seqs <- do.call(c, unname(seq_list[valid]))
writeXStringSet(unexplained_seqs, file.path(data_dir, "unexplained_chimera_for_uniprot.fasta"))

cat("\nSaved to", file.path(data_dir, "unexplained_chimera_for_uniprot.fasta"), "\n")
cat("Next: run download_uniprot_p450.sh (if not already done), then:\n\n")
cat("  blastx -query", file.path(data_dir, "unexplained_chimera_for_uniprot.fasta"),
    "-db uniprot_p450_db -outfmt '6 qseqid sseqid pident length evalue bitscore qstart qend qlen stitle' -max_target_seqs 3 -evalue 1e-5 -out",
    file.path(data_dir, "unexplained_chimera_vs_uniprot.tsv"), "\n")
