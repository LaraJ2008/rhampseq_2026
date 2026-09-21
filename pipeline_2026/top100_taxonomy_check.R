## =============================================================================
## Taxonomy check on the top-100 abundant reads, using local Swiss-Prot+taxdb
## =============================================================================
## Same approach that worked for identifying the divergent RH.E04EE9/RH.8B25A0
## sequence earlier -- blastx (translated nucleotide query vs protein db),
## this time requesting scientific names and taxids directly in the output,
## and run across all 100 reads at once since local BLAST handles this easily.
##
## CAVEAT: Swiss-Prot is small and biased toward well-studied model organisms.
## A best hit here means "closest known relative Swiss-Prot contains," not
## necessarily the read's true organism of origin -- same limitation we hit
## with the very first Caulobacter Swiss-Prot search early in this investigation.
## =============================================================================

## ---- Step 1: run in bash --------------------------------------------------
## (paste this whole block directly, no line breaks, to avoid paste-mangling)

# blastx -query dada2_bypass_check/top100_abundant.fasta -db swissprot_db_clean/swissprot -outfmt '6 qseqid sseqid pident length evalue bitscore stitle sscinames staxids' -max_target_seqs 3 -evalue 1e-5 -out dada2_bypass_check/top100_vs_swissprot_taxonomy.tsv

## ---- Step 2: parse and summarize in R --------------------------------------

library(data.table)
library(stringr)

hits <- fread("dada2_bypass_check/top100_vs_swissprot_taxonomy.tsv", header = FALSE,
              col.names = c("qseqid","sseqid","pident","length","evalue",
                            "bitscore","stitle","sscinames","staxids"),
              quote = "")   # protein titles can contain quote characters

hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]

## keep only the single best hit per read, by bitscore
best_hit <- hits[order(-bitscore)][, .SD[1], by = qseqid]

cat("Reads with at least one Swiss-Prot hit:", nrow(best_hit), "out of 100\n")
cat("Reads with NO hit at all (nothing in Swiss-Prot resembles them):",
    100 - nrow(best_hit), "\n\n")

## ---- Organism distribution across the top 100 -------------------------------

cat("=== Top organisms represented (best hit per read) ===\n")
print(best_hit[, .N, by = sscinames][order(-N)])

## ---- Identity distribution -- low identity = "real hit but distant relative" -

cat("\n=== Identity distribution of best hits ===\n")
print(summary(best_hit$pident))

## flag reads whose best Swiss-Prot hit is weak (low identity), meaning even
## the "closest" thing Swiss-Prot has isn't very close -- these are the reads
## most likely to represent something genuinely novel/under-characterized
weak_hits <- best_hit[pident < 50]
cat("\nReads with <50% identity to their best Swiss-Prot hit (likely novel/",
    "poorly-represented genes):", nrow(weak_hits), "\n")

fwrite(best_hit, "dada2_bypass_check/top100_taxonomy_summary.csv")
cat("\nSaved full results to dada2_bypass_check/top100_taxonomy_summary.csv\n")
