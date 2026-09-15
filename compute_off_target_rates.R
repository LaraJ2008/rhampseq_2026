library(data.table)
library(stringr)
library(Biostrings)

## ---- CONFIG: two separate locations -------------------------------------------
## de_rep_dir: where the (large) *_derep.fasta files actually live (e.g. D drive)
## out_dir:    where everything else lives -- blast .tsv files, the bash
##             script's off_target_summary.csv, and this script's final output
de_rep_dir <- "D:/NSF_FAST_poc_extracted/D_off_target_check/"   # <-- EDIT THIS
out_dir    <- "C_off_target_check"                   # <-- EDIT THIS if needed

## ---- Find all samples that were processed -----------------------------------
derep_files <- list.files(de_rep_dir, pattern = "_derep\\.fasta$", full.names = TRUE)
samples <- str_remove(basename(derep_files), "_derep\\.fasta$")

cat("Samples found:", length(samples), "\n")
print(samples)

## ---- Total read counts (recompute here for the summary table, or read from --
## the bash script's own summary CSV if it saved one -- using that directly
## if present, since it already has total_reads captured)

if (file.exists(file.path(out_dir, "off_target_summary.csv"))) {
  bash_summary <- fread(file.path(out_dir, "off_target_summary.csv"))
} else {
  bash_summary <- NULL
}

## ---- Process each sample -------------------------------------------------------
## IMPORTANT: run this entire lapply(...) call as ONE block, not line by line --
## `s` only has a value while lapply is actively iterating; pasting individual
## lines from inside the function body on their own will always error with
## "object 's' not found", since nothing ever assigns it in that context.

results <- rbindlist(lapply(samples, function(s) {
  
  derep_path <- file.path(de_rep_dir, paste0(s, "_derep.fasta"))
  blast_path <- file.path(out_dir, paste0(s, "_vs_primers.tsv"))
  
  # get abundance per unique sequence from the fasta headers (;size=N)
  seqs <- readDNAStringSet(derep_path)
  seq_dt <- data.table(
    qseqid    = str_remove(names(seqs), ";size=\\d+$"),
    abundance = as.integer(str_extract(names(seqs), "(?<=size=)\\d+"))
  )
  
  reads_examined <- sum(seq_dt$abundance)
  
  if (!file.exists(blast_path) || file.info(blast_path)$size == 0) {
    # no primer hits at all for this sample's examined sequences
    matched_ids <- character(0)
  } else {
    blast_hits <- fread(blast_path, header = FALSE,
                        col.names = c("qseqid","sseqid","pident","length",
                                      "mismatch","evalue","bitscore"))
    blast_hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]
    blast_filtered <- blast_hits[length >= 18 & pident >= 90]
    matched_ids <- unique(blast_filtered$qseqid)
  }
  
  seq_dt[, matched := qseqid %in% matched_ids]
  reads_matched <- sum(seq_dt[matched == TRUE, abundance])
  reads_unmatched <- reads_examined - reads_matched
  
  data.table(
    sample = s,
    reads_examined = reads_examined,
    reads_matched = reads_matched,
    reads_unmatched = reads_unmatched,
    off_target_rate_among_examined = round(100 * reads_unmatched / reads_examined, 2)
  )
}))

## ---- Merge in total read counts and compute coverage --------------------------

if (!is.null(bash_summary)) {
  results <- merge(results, bash_summary[, .(sample, total_reads)], by = "sample", all.x = TRUE)
  results[, pct_of_total_examined := round(100 * reads_examined / total_reads, 2)]
  setcolorder(results, c("sample", "total_reads", "reads_examined",
                         "pct_of_total_examined", "reads_matched",
                         "reads_unmatched", "off_target_rate_among_examined"))
}

setorder(results, sample)

print(results)

fwrite(results, file.path(out_dir, "final_off_target_rates.csv"))

cat("\n=== Overall (pooled across all samples, weighted by examined reads) ===\n")
cat("Off-target rate:",
    round(100 * sum(results$reads_unmatched) / sum(results$reads_examined), 2),
    "%\n")

