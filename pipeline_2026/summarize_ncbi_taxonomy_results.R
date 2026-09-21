library(data.table)
library(stringr)

## Works directly from top5_hits_per_target_with_taxonomy.csv -- no need for
## a separate best_hit_per_target_with_taxonomy.csv file.

top5 <- fread("top5_hits_per_target_with_taxonomy.csv")

## ---- 1. Best hit per query (hit_num == 1, already NCBI-sorted, but confirm --
## via max bitscore just to be safe) ------------------------------------------

best_hit <- top5[order(-bitscore)][, .SD[1], by = qseqid]
setorder(best_hit, -bitscore)

cat("Queries represented:", uniqueN(top5$qseqid), "\n\n")
cat("=== Best hit per query ===\n")
print(best_hit[, .(qseqid, sciname, title, pident, evalue, bitscore)])

fwrite(best_hit, "best_hit_per_target_derived.csv")

## ---- 2. Organism distribution across all best hits ---------------------------

cat("\n=== Organism distribution (best hit per query) ===\n")
print(best_hit[, .N, by = sciname][order(-N)])

## ---- 3. Flag anything in the Caulobacteraceae family, not just genus --------
## Caulobacter, Phenylobacterium, and a few other genera all sit in this same
## family -- worth seeing the broader picture, not just exact-genus matches

caulobacteraceae_genera <- c("Caulobacter", "Phenylobacterium", "Brevundimonas",
                              "Asticcacaulis", "Caulobacterales")

best_hit[, is_caulobacteraceae_relative := str_detect(sciname,
                        paste(caulobacteraceae_genera, collapse = "|"))]

cat("\n=== Targets whose best hit is Caulobacter OR a close relative",
    "(same family, Caulobacteraceae) ===\n")
print(best_hit[is_caulobacteraceae_relative == TRUE,
               .(qseqid, sciname, title, pident, evalue)])

## ---- 4. Cross-check against the confirmed 31 Caulobacter-genome-matching ---
## sources, if that file is available

if (file.exists("caulobacter_matching_targets.txt")) {
  caulobacter_sources <- readLines("caulobacter_matching_targets.txt")
  best_hit[, is_known_caulobacter_source := qseqid %in% caulobacter_sources]

  cat("\n=== Of your", length(caulobacter_sources), "already-confirmed",
      "Caulobacter sources, how many appear among these 37 queries,",
      "and what does NCBI call them? ===\n")
  print(best_hit[is_known_caulobacter_source == TRUE,
                 .(qseqid, sciname, title, pident, evalue)])

  cat("\n=== Any Caulobacteraceae-relative hit NOT already in your confirmed",
      "31 -- i.e. a related genus your genome-based check would have missed ===\n")
  print(best_hit[is_caulobacteraceae_relative == TRUE & is_known_caulobacter_source == FALSE,
               .(qseqid, sciname, title, pident, evalue)])
} else {
  cat("\n(caulobacter_matching_targets.txt not found in this directory --",
      "skipping cross-check against the confirmed-31 list)\n")
}

print(top5[qseqid == "A00405_539_HLLNCDSX3_4_2144_25943_22498_1_291_+_b63a675850008016b91f6d9617d49ae1"])
