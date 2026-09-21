library(data.table)
library(stringr)
library(Biostrings)

unexamined_dir <- "D:/NSF_FAST_poc_extracted/D_low_abundance_analysis/unexamined_smplvsprimers_fasta/" 

derep_files <- list.files(unexamined_dir, pattern = "_unexamined\\.fasta$", full.names = TRUE)
samples <- str_remove(basename(derep_files), "_unexamined\\.fasta$")

cat("Samples with unexamined sequences found:", length(samples), "\n\n")

## ---- Process each sample -------------------------------------------------------

all_matched_seqs <- list()
all_pair_status <- list()

results <- rbindlist(lapply(samples, function(s) {

  fasta_path <- file.path(unexamined_dir, paste0(s, "_unexamined.fasta"))
  blast_path <- file.path(unexamined_dir, paste0(s, "_unexamined_vs_primers.tsv"))

  seqs <- readDNAStringSet(fasta_path)
  seq_dt <- data.table(
    qseqid    = str_remove(names(seqs), ";size=\\d+$"),
    abundance = as.integer(str_extract(names(seqs), "(?<=size=)\\d+"))
  )

  reads_unexamined <- sum(seq_dt$abundance)
  n_unique_unexamined <- nrow(seq_dt)

  if (!file.exists(blast_path) || file.info(blast_path)$size == 0) {
    matched_ids <- character(0)
    pair_status_dt <- data.table(qseqid = character(0), pair_status = character(0),
                                  f_targets = character(0), r_targets = character(0))
  } else {
    blast_hits <- fread(blast_path, header = FALSE,
                         col.names = c("qseqid","sseqid","pident","length",
                                       "mismatch","qstart","qend","sstart",
                                       "send","evalue","bitscore"))
    blast_hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]
    blast_filtered <- blast_hits[length >= 18 & pident >= 90]
    matched_ids <- unique(blast_filtered$qseqid)

    blast_filtered[, target_id := str_remove(sseqid, "_[FR]$")]
    blast_filtered[, orientation := str_extract(sseqid, "[FR]$")]

    ## classify each matched read: matched_pair (F and R hit the SAME target),
    ## chimera (F and R hit DIFFERENT targets), or F_only / R_only (single end)
    pair_status_dt <- blast_filtered[, {
      f_ids <- unique(target_id[orientation == "F"])
      r_ids <- unique(target_id[orientation == "R"])
      status <- if (length(f_ids) > 0 && length(r_ids) > 0) {
        if (length(intersect(f_ids, r_ids)) > 0) "matched_pair" else "chimera"
      } else if (length(f_ids) > 0) {
        "F_only"
      } else {
        "R_only"
      }
      list(pair_status = status,
           f_targets = paste(f_ids, collapse = "; "),
           r_targets = paste(r_ids, collapse = "; "))
    }, by = qseqid]
  }

  seq_dt[, matched_a_primer := qseqid %in% matched_ids]
  seq_dt <- merge(seq_dt, pair_status_dt, by = "qseqid", all.x = TRUE)

  # stash the actual sequences + classification of matched-but-rare reads
  if (length(matched_ids) > 0) {
    matched_seqs_this_sample <- seqs[str_remove(names(seqs), ";size=\\d+$") %in% matched_ids]
    all_matched_seqs[[s]] <<- matched_seqs_this_sample
    all_pair_status[[s]] <<- seq_dt[matched_a_primer == TRUE,
                                     .(sample = s, qseqid, abundance,
                                       pair_status, f_targets, r_targets)]
  }

  data.table(
    sample = s,
    n_unique_unexamined = n_unique_unexamined,
    reads_unexamined = reads_unexamined,
    n_unique_matched_a_primer = sum(seq_dt$matched_a_primer),
    reads_matched_a_primer = sum(seq_dt[matched_a_primer == TRUE, abundance]),
    pct_of_unexamined_that_match_a_primer = round(
      100 * sum(seq_dt[matched_a_primer == TRUE, abundance]) / reads_unexamined, 2),
    n_matched_pair = sum(seq_dt$pair_status == "matched_pair", na.rm = TRUE),
    n_chimera      = sum(seq_dt$pair_status == "chimera", na.rm = TRUE),
    n_single_end   = sum(seq_dt$pair_status %in% c("F_only","R_only"), na.rm = TRUE)
  )
}))

print(results)

fwrite(results, file.path(unexamined_dir, "unexamined_summary.csv"))

cat("\n=== Overall ===\n")
cat("Total unexamined reads across all samples:", sum(results$reads_unexamined), "\n")
cat("Of those, matching at least one primer:", sum(results$reads_matched_a_primer), "\n")
cat("Percent:", round(100 * sum(results$reads_matched_a_primer) / sum(results$reads_unexamined), 2), "%\n")

## ---- Save the actual matched-but-rare sequences for direct inspection --------
## This is the population directly relevant to the "primer promiscuity" question
## -- reads that a primer clearly recognizes, but that never reached abundance
## 5+ as an exact duplicate. Worth BLASTing these (whole-read, not just the
## primer ends) against the confirmed P450 set or nr to see whether their
## "internal" sequence matches the SAME P450 gene the primer was designed for,
## or something genuinely different.

if (length(all_matched_seqs) > 0) {

  ## defensive check: make sure every stored entry is actually a proper
  ## DNAStringSet before trying to combine them -- drops anything malformed
  ## (e.g. an unexpected NULL or wrong type from an edge-case sample) rather
  ## than letting one bad entry break the whole combine step
  entry_classes <- sapply(all_matched_seqs, function(x) class(x)[1])
  cat("Entry classes found:", paste(unique(entry_classes), collapse = ", "), "\n")

  valid_entries <- Filter(function(x) is(x, "XStringSet") && length(x) > 0,
                           all_matched_seqs)
  n_dropped <- length(all_matched_seqs) - length(valid_entries)
  if (n_dropped > 0) {
    cat("WARNING: dropped", n_dropped,
        "malformed/empty entries before combining sequences\n")
  }

  combined_matched <- do.call(c, unname(valid_entries))
  writeXStringSet(combined_matched, file.path(unexamined_dir, "matched_but_rare_sequences.fasta"))
  cat("\nSaved", length(combined_matched),
      "matched-but-rare sequences to matched_but_rare_sequences.fasta\n")

  combined_pair_status <- rbindlist(all_pair_status)
  fwrite(combined_pair_status, file.path(unexamined_dir, "matched_but_rare_pair_status.csv"))

  cat("\n=== Pair-status breakdown across all matched-but-rare reads ===\n")
  cat("(matched_pair = F and R hit the SAME target -- looks like a real amplicon)\n")
  cat("(chimera = F and R hit DIFFERENT targets -- possible primer-jumping artifact)\n")
  cat("(F_only/R_only = only one end matched, likely a longer amplicon that didn't\n")
  cat(" read through to the other primer)\n\n")
  print(table(combined_pair_status$pair_status))

  n_chimera <- sum(combined_pair_status$pair_status == "chimera")
  if (n_chimera > 0) {
    cat("\n=== Chimera reads: which target combinations are involved? ===\n")
    print(combined_pair_status[pair_status == "chimera",
                                .(sample, qseqid, abundance, f_targets, r_targets)])
  }

  cat("\nNext step to test the promiscuity hypothesis: BLAST this file's FULL\n",
      "sequences (not just primer regions) against your confirmed P450 set\n",
      "or NCBI nr, and compare each hit's identity to what the ORIGINAL primer\n",
      "target was designed from -- agreement supports 'just a rare real read',\n",
      "disagreement supports 'primer caught a different, related P450'.\n",
      "Focus this follow-up especially on the matched_pair reads -- chimeras are\n",
      "more likely a PCR/library artifact than genuine biological promiscuity.\n")
}

#~~~~~~~~~~~~~~~~~~~What's up with thousands of chimeras? Probably low validity matches~~~~~~~~~~

blast_all <- rbindlist(lapply(list.files(unexamined_dir, pattern = "_vs_primers\\.tsv$", full.names = TRUE),
                              function(f) {
                                dt <- fread(f, header = FALSE,
                                            col.names = c("qseqid","sseqid","pident","length","mismatch",
                                                          "qstart","qend","sstart","send","evalue","bitscore"))
                                dt[, qseqid := str_remove(qseqid, ";size=\\d+$")]
                                dt
                              }))
#~~~~~~~~~~STRINGENCY OF PRIMER MATCH BEAKDOWN; TWO OPTIONS~~~~~~~~~~~~~~~~~~~

#blast_filtered <- blast_all[length >= 18 & pident >= 90]
blast_filtered <- blast_all[length >= 22 & pident >= 95]
chimera_reads <- combined_pair_status[pair_status == "chimera", qseqid]
pair_reads    <- combined_pair_status[pair_status == "matched_pair", qseqid]

cat("=== Match quality: chimera-classified reads ===\n")
print(summary(blast_filtered[qseqid %in% chimera_reads, .(pident, length)]))

cat("\n=== Match quality: matched_pair-classified reads ===\n")
print(summary(blast_filtered[qseqid %in% pair_reads, .(pident, length)]))


#~~~~~~~~~~~~~~~~~~Chimeric pairs look prett robust, are they PCR issues or homology?~~~~~~~~~~~~~~

chimera_dt <- combined_pair_status[pair_status == "chimera"]

# dereplicate to distinct sequences within the chimera set, prioritize the
# most abundant ones as a practical, defensible subset to test
chimera_dt[, seq_key := qseqid]  # already unique reads; abundance tells us weight
top_chimera <- chimera_dt[order(-abundance)][1:30]  # top 30 most abundant

# pull their actual sequences
all_unexamined_seqs <- list.files(unexamined_dir, pattern = "_unexamined\\.fasta$", full.names = TRUE)
seqs_by_sample <- lapply(setNames(all_unexamined_seqs, str_remove(basename(all_unexamined_seqs), "_unexamined\\.fasta$")),
                         readDNAStringSet)

pull_seq <- function(sample, qseqid) {
  s <- seqs_by_sample[[sample]]
  s[str_remove(names(s), ";size=\\d+$") == qseqid]
}
chimera_seqs <- do.call(c, Map(pull_seq, top_chimera$sample, top_chimera$qseqid))
writeXStringSet(chimera_seqs, file.path(unexamined_dir, "top30_chimera_reads.fasta"))


#grab random set of 'chimera' reads to blast against p450s~~~~~~~~~~~~

set.seed(42)
chimera_dt_shuffled <- chimera_dt[sample(.N)]
top_chimera <- chimera_dt_shuffled[order(-abundance)][1:30]

table(top_chimera$sample)
table(top_chimera$abundance)

seq_list <- Map(pull_seq, top_chimera$sample, top_chimera$qseqid)
valid <- sapply(seq_list, function(x) is(x, "XStringSet") && length(x) == 1)
chimera_seqs <- do.call(c, unname(seq_list[valid]))
writeXStringSet(chimera_seqs, file.path(unexamined_dir, "top30_chimera_reads.fasta"))


