library(data.table)
library(stringr)

## =============================================================================
## Clean up NCBI ClusteredNR/nr blastx hit table (raw CSV export, no header)
## =============================================================================
## Column structure (reconstructed from the raw file -- NCBI's web CSV export
## does NOT include a header row or organism/description column, just raw
## accession numbers):
##   qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
##   sstart, send, evalue, bitscore, positive_pct, qframe, sframe
## =============================================================================

hits <- fread("AMU87N88014-Alignment-HitTable_top38targets.csv", header = FALSE,
              col.names = c("qseqid","sseqid","pident","length","mismatch",
                            "gapopen","qstart","qend","sstart","send",
                            "evalue","bitscore","positive_pct","qframe","sframe"))

# clean up query ID (strip the ;size= suffix if present, harmless if not)
hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]

cat("Total hit rows:", nrow(hits), "\n")
cat("Distinct queries (should be up to 38):", uniqueN(hits$qseqid), "\n\n")

## ---- 1. Best hit per query (by bitscore) -------------------------------------

best_hit <- hits[order(-bitscore)][, .SD[1], by = qseqid]
setorder(best_hit, -bitscore)

cat("=== Best hit per query ===\n")
print(best_hit[, .(qseqid, sseqid, pident, length, evalue, bitscore)])

fwrite(best_hit, "best_hit_per_target.csv")

## ---- 2. Top 5 hits per query (useful for eyeballing consistency/ambiguity) ---

top5 <- hits[order(qseqid, -bitscore)][, head(.SD, 5), by = qseqid]
fwrite(top5, "top5_hits_per_target.csv")

cat("\n=== Summary stats on best hits ===\n")
cat("Best-hit identity range:", round(min(best_hit$pident), 1), "-",
    round(max(best_hit$pident), 1), "%\n")
cat("Queries with best-hit identity < 50% (likely novel/poorly represented):",
    sum(best_hit$pident < 50), "\n")

cat("\n=== IMPORTANT ===\n")
cat("This file has NO organism/species information -- just raw accession\n")
cat("numbers (e.g. 'WP_298603505.1'). To get real taxonomy, the accessions\n")
cat("in best_hit$sseqid need to be looked up separately via NCBI. Next step:\n")
cat("save the unique accession list and run an esummary/efetch lookup:\n\n")

writeLines(unique(best_hit$sseqid), "accessions_for_taxonomy_lookup.txt")
cat("Saved", uniqueN(best_hit$sseqid), "unique accessions to",
    "accessions_for_taxonomy_lookup.txt\n")
