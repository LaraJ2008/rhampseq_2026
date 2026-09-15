## =============================================================================
## From filtered BLAST hits: find TRUE primer-pair matches against Caulobacter
## =============================================================================
## A single primer hit means little on its own -- short partial matches occur
## by chance across large genomes. Real evidence of a Caulobacter-amplifiable
## primer pair requires BOTH the F and R primer from the SAME target_id to
## hit the SAME genome record, in the correct relative orientation, within a
## realistic amplicon-sized distance of each other.
## =============================================================================

library(data.table)
library(stringr)

hits <- fread("primers_vs_caulobacter_filtered.tsv", header = FALSE,
              col.names = c("qseqid","sseqid","pident","length","mismatch",
                            "gapopen","qstart","qend","sstart","send",
                            "evalue","bitscore"))

cat("Filtered hits loaded:", nrow(hits), "\n")

## ---- split qseqid into target_id + orientation (F/R) -----------------------
hits[, target_id := str_remove(qseqid, "_[FR]$")]
hits[, orientation := str_extract(qseqid, "[FR]$")]

## ---- BLAST reports sstart > send for reverse-strand hits; normalize --------
hits[, strand := ifelse(sstart < send, "+", "-")]
hits[, pos_min := pmin(sstart, send)]
hits[, pos_max := pmax(sstart, send)]

f_hits <- hits[orientation == "F"]
r_hits <- hits[orientation == "R"]

## ---- pair up F and R hits sharing the same target_id AND genome record ----
paired <- merge(f_hits, r_hits, by = c("target_id","sseqid"),
                 suffixes = c("_F","_R"), allow.cartesian = TRUE)

## ---- realistic amplicon check: opposite strands, sensible distance apart --
## rhAmpSeq amplicons in this panel are short (observed ORF lengths earlier
## topped out around 276bp) -- use a generous 50-1000bp window to avoid being
## overly strict, but strand orientation must be correct (F on one strand,
## R on the other, pointing toward each other).
MAX_AMPLICON <- 1000
MIN_AMPLICON <- 20

paired[, amplicon_dist := ifelse(pos_min_F < pos_min_R,
                                   pos_min_R - pos_max_F,
                                   pos_min_F - pos_max_R)]

true_pairs <- paired[
  strand_F != strand_R &
  amplicon_dist > MIN_AMPLICON & amplicon_dist < MAX_AMPLICON
]

cat("\n=== Candidate TRUE primer-pair matches against Caulobacter genomes ===\n")
cat("N =", nrow(true_pairs), "\n\n")

if (nrow(true_pairs) > 0) {
  print(true_pairs[, .(target_id, sseqid, pident_F, pident_R, length_F, length_R,
                        strand_F, strand_R, amplicon_dist)])
  fwrite(true_pairs, "confirmed_caulobacter_primer_pairs.csv")
  cat("\nSaved to confirmed_caulobacter_primer_pairs.csv\n")
  cat("\n*** These are targets whose primers WOULD plausibly amplify real\n",
      "Caulobacter genomic sequence -- worth cross-checking against whether\n",
      "these specific target_ids ended up in your final 415-target panel or\n",
      "were dropped during multiplex QC. If they're in the panel but showing\n",
      "zero reads in your fastq data, that points toward a true wet-lab\n",
      "amplification failure rather than a design-stage dropout.\n")
} else {
  cat("\n*** No F/R primer pair from any target lands on the same Caulobacter\n",
      "genome in correct orientation at a realistic amplicon distance.\n",
      "Individual short partial matches may still show up in the filtered\n",
      "hits table (this is expected/normal), but none of them combine into a\n",
      "real amplifiable pair. This is genuine evidence that your primers, as\n",
      "designed, do not target Caulobacter sequence -- consistent with\n",
      "primers having been designed from non-Caulobacter source contigs in\n",
      "the first place, not a cross-reactivity or off-target concern.\n")
}


involved_ids <- unique(c(paste0(true_pairs$target_id, "_F"),
                         paste0(true_pairs$target_id, "_R")))
writeLines(involved_ids, "involved_primer_ids.txt")
