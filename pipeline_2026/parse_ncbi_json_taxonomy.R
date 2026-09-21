install.packages("jsonlite")
install.packages("Rtools
  ")
library(jsonlite)
library(data.table)
library(stringr)

## =============================================================================
## Parse NCBI blastx JSON export -- includes organism names directly, no
## separate taxonomy lookup needed (unlike the plain Hit Table CSV)
## =============================================================================

raw <- fromJSON("AMU87N88014-Alignment.json", simplifyVector = FALSE)
queries <- raw$BlastOutput2

cat("Queries with results:", length(queries), "(expect up to 38 -- any missing",
    "likely had zero hits at all)\n\n")

## ---- 1. Flatten into one long table: one row per (query, hit) --------------

all_hits <- rbindlist(lapply(queries, function(q) {
  search <- q$report$results$search
  qid <- search$query_title
  qlen <- search$query_len
  hits <- search$hits

  if (length(hits) == 0) {
    return(data.table(qseqid = qid, qlen = qlen, hit_num = NA, accession = NA,
                       title = "NO HITS FOUND", sciname = NA, taxid = NA,
                       pident = NA, align_len = NA, evalue = NA, bitscore = NA))
  }

  rbindlist(lapply(hits, function(h) {
    desc <- h$description[[1]]   # top/representative description for this hit
    hsp  <- h$hsps[[1]]          # best HSP (NCBI already orders these)

    data.table(
      qseqid    = qid,
      qlen      = qlen,
      hit_num   = h$num,
      accession = desc$accession,
      title     = desc$title,
      sciname   = desc$sciname,
      taxid     = desc$taxid,
      pident    = round(100 * hsp$identity / hsp$align_len, 2),
      align_len = hsp$align_len,
      evalue    = hsp$evalue,
      bitscore  = hsp$bit_score
    )
  }))
}))

cat("Total hit rows:", nrow(all_hits), "\n\n")

## ---- 2. Best hit per query ----------------------------------------------------

best_hit <- all_hits[order(-bitscore)][, .SD[1], by = qseqid]
setorder(best_hit, -bitscore)

cat("=== Best hit per query (organism included) ===\n")
print(best_hit[, .(qseqid, sciname, title, pident, evalue, bitscore)])

fwrite(best_hit, "best_hit_per_target_with_taxonomy.csv")

## ---- 3. Organism distribution across all 37 best hits -----------------------

cat("\n=== Organism distribution (best hit per query) ===\n")
print(best_hit[, .N, by = sciname][order(-N)])

## ---- 4. Cross-check: do the 4 known-Caulobacter targets actually show ------
## Caulobacter here? Good validation that this identification approach works.

caulobacter_sources <- readLines("caulobacter_matching_targets.txt")

best_hit[, is_known_caulobacter_source := qseqid %in% caulobacter_sources]

cat("\n=== Validation: known Caulobacter-source queries -- do their best",
    "NCBI hits actually say Caulobacter? ===\n")
print(best_hit[is_known_caulobacter_source == TRUE,
               .(qseqid, sciname, title, pident, evalue)])

## ---- 5. Top 5 hits per query (for eyeballing ambiguous/scattered results) ---

top5 <- all_hits[!is.na(hit_num)][order(qseqid, -bitscore)][, head(.SD, 5), by = qseqid]
fwrite(top5, "top5_hits_per_target_with_taxonomy.csv")

cat("\n=== Queries with NO hits at all ===\n")
print(all_hits[title == "NO HITS FOUND", .(qseqid)])

cat("\nSaved:\n")
cat("  best_hit_per_target_with_taxonomy.csv\n")
cat("  top5_hits_per_target_with_taxonomy.csv\n")
