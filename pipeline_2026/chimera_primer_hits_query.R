library(data.table)
library(stringr)

hits <- fread("unresolved_relaxed_primer_hits.tsv", header = FALSE,
              col.names = c("qseqid","sseqid","pident","length","mismatch",
                            "qstart","qend","sstart","send","evalue","bitscore"))
hits[, qseqid := str_remove(qseqid, ";size=\\d+$")]
hits[, target_id := str_remove(sseqid, "_[FR]$")]
hits[, orientation := str_extract(sseqid, "[FR]$")]

## best (highest bitscore) relaxed match per read, regardless of how weak
best_relaxed <- hits[order(-bitscore)][, .SD[1], by = qseqid]

cat("Reads with SOME primer signal, even relaxed:", nrow(best_relaxed),
    "out of 91\n")
cat("Reads with NO primer signal at all, even relaxed:", 91 - nrow(best_relaxed), "\n\n")

cat("=== Match quality distribution of these 'best available' relaxed hits ===\n")
print(summary(best_relaxed[, .(pident, length, mismatch)]))

## show the full picture per read: identity, length, mismatches, which primer
cat("\n=== Full per-read best relaxed match ===\n")
print(best_relaxed[order(-pident), .(qseqid, target_id, orientation,
                                      pident, length, mismatch, evalue)])

## Which primer/target shows up most often as the best (even weak) explanation?
cat("\n=== Which primers most often explain these misprimes? ===\n")
print(best_relaxed[, .N, by = target_id][order(-N)])

fwrite(best_relaxed, "unresolved_relaxed_best_hit.csv")
