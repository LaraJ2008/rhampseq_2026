library(data.table)

library(data.table)
library(stringr)

## ---- 1. Load demux counts (this file is fine as-is) --------------------------
## Note: cutadapt paired demux names each output by the FORWARD (-g) primer's
## matched adapter name -- there's no separate "_R" row per target, the "_F"
## name represents the whole read PAIR. This is expected, not missing data.

demux_raw <- fread("demux_counts_21_targets.tsv", header = FALSE,
                    col.names = c("qseqid", "demux_reads"))
demux_raw[, target_id := str_remove(qseqid, "_F$")]
demux <- demux_raw[, .(target_id, demux_reads)]

## ---- 2. Load pileup counts, robust to a merged target_id/region column ------
## If confirmed_targets.bed had a delimiter problem (e.g. leftover CRLF from
## an earlier WSL/Windows edit breaking the tab-split in the bash while-loop),
## target_id and genomic_region may have merged together with no separator,
## producing a 2-column file instead of 3. Detect and handle both cases.

pileup_lines <- readLines("pileup_counts_21_targets.tsv")
n_fields <- str_count(pileup_lines[1], "\t") + 1

if (n_fields == 3) {
  pileup_raw <- fread("pileup_counts_21_targets.tsv", header = FALSE,
                       col.names = c("target_id", "genomic_region", "pileup_reads"))
} else {
  cat("NOTE: pileup file has", n_fields, "tab-separated fields, not the expected 3 --",
      "parsing the merged target_id/region column with a regex fallback.\n\n")

  raw2col <- fread("pileup_counts_21_targets.tsv", header = FALSE,
                    col.names = c("merged", "pileup_reads"))

  # target_id pattern: "RH." + hex-ish characters; genomic_region pattern:
  # an optional "NZ_" prefix, then accession letters/digits, a period+version,
  # a colon, then start-end coordinates
  extracted <- str_match(raw2col$merged,
                          "^(RH\\.[A-Z0-9]+)((?:NZ_)?[A-Z]+[0-9.]+:[0-9]+-[0-9]+)$")

  pileup_raw <- data.table(
    target_id      = extracted[, 2],
    genomic_region = extracted[, 3],
    pileup_reads   = raw2col$pileup_reads
  )

  n_failed <- sum(is.na(pileup_raw$target_id))
  if (n_failed > 0) {
    cat("WARNING:", n_failed, "rows didn't match the expected pattern -- inspect these:\n")
    print(raw2col[is.na(pileup_raw$target_id)])
  }
}

## ---- 2b. Collapse duplicate CP/NZ_CP accession pairs ------------------------
## NCBI's esearch/efetch pulled both a GenBank-style (e.g. CP033875.1) and a
## RefSeq-style (NZ_CP033875.1) record for several of these genomes -- same
## underlying sequence, counted twice under two accessions. Normalize by
## stripping "NZ_" before aggregating, so real coverage per target isn't
## artificially split/duplicated across both copies.

pileup_raw[, accession_normalized := str_remove(genomic_region, "^NZ_")]

pileup <- pileup_raw[, .(pileup_reads = sum(pileup_reads, na.rm = TRUE)),
                      by = .(target_id, accession_normalized)]
pileup <- pileup[, .(pileup_reads = sum(pileup_reads)), by = target_id]
# ^ second aggregation collapses multiple genomic regions per target_id
#   (e.g. if a target's source sequence matched several distinct genome
#   locations) into one total pileup count per target

## ---- 3. Join into one table --------------------------------------------------

results <- merge(demux, pileup, by = "target_id", all = TRUE)
results[is.na(demux_reads), demux_reads := 0]
results[is.na(pileup_reads), pileup_reads := 0]

## ---- 4. Apply the interpretation logic directly as a column -----------------

results[, interpretation := fcase(
  demux_reads == 0 & pileup_reads == 0,
    "Genuine amplification failure",

  pileup_reads > 0 & demux_reads == 0,
    "Reads exist but primer didn't recognize them (mismatch/demux issue)",

  demux_reads > 0 & pileup_reads > 0,
    "Caulobacter reads genuinely present -- original 'zero hits' was a classification problem",

  default = "Unexpected pattern -- inspect manually"
)]

setorder(results, -pileup_reads, -demux_reads)

## ---- 5. Print and save --------------------------------------------------------

print(results)

fwrite(results, "combined_results_21_targets.csv")

cat("\n=== Summary ===\n")
print(table(results$interpretation))

results[interpretation == "Unexpected pattern -- inspect manually"]
