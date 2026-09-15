## =============================================================================
## Get exact genomic coordinates for the 21 confirmed Caulobacter-derived targets
## =============================================================================
## We already have two pieces that need joining:
##  1. primers_vs_caulobacter_candidates_filtered.tsv -- confirms which
##     target_id's F/R primers match which SOURCE candidate sequence
##     (sseqid here = a source read ID from the pre-filter 511-candidate pool)
##  2. targets_vs_caulobacter_genomes.tsv -- confirms where that SAME source
##     candidate sequence (as qseqid there) lands on an actual Caulobacter
##     genome (sseqid = genome accession, sstart/send = coordinates)
##
## Joining on the shared source-read ID gives us genome coordinates for each
## of the 21 confirmed targets -- a BED region we can directly check for real
## read pileup, independent of whether primer-based demultiplexing works.
## =============================================================================

library(data.table)
library(stringr)

## ---- 1. Load and find true pairs (target + matched source) -----------------

pv_c <- fread("primers_vs_caulobacter_candidates_filtered.tsv", header = FALSE,
              col.names = c("qseqid","sseqid","pident","length","mismatch",
                            "qstart","qend","sstart","send","evalue","bitscore"))

pv_c[, target_id := str_remove(qseqid, "_[FR]$")]
pv_c[, orientation := str_extract(qseqid, "[FR]$")]

# true pair = same target_id + same source sseqid hit by BOTH F and R
pair_counts <- pv_c[, .(n_orientations = uniqueN(orientation)), by = .(target_id, sseqid)]
true_pairs <- pair_counts[n_orientations == 2]

cat("Confirmed true pairs:", nrow(true_pairs), "\n")
print(true_pairs$target_id)

## ---- 2. Join to genome coordinates from the target-vs-genome BLAST ---------

genome_hits <- fread("targets_vs_caulobacter_genomes.tsv", header = FALSE,
                      col.names = c("qseqid","sseqid","pident","length",
                                    "evalue","bitscore","sstart","send"))
# here qseqid = the source candidate read ID (same ID space as pv_c$sseqid)

coords <- merge(true_pairs, genome_hits,
                 by.x = "sseqid", by.y = "qseqid",
                 suffixes = c("_primer_match","_genome"))

coords[, bed_start := pmin(sstart, send) - 1]  # BED is 0-based, half-open
coords[, bed_end   := pmax(sstart, send)]
coords[, genome_accession := sseqid_genome]

cat("\nTargets with resolved genome coordinates:", uniqueN(coords$target_id), "\n")
print(coords[, .(target_id, genome_accession, bed_start, bed_end, pident_genome)])

## ---- 3. Write BED file for samtools region-based read counting -------------

# if a target has multiple hit regions (e.g. from repeated/near-identical
# genome copies), keep all of them -- worth checking coverage at every one
fwrite(coords[, .(genome_accession, bed_start, bed_end, target_id)],
       "confirmed_targets.bed", sep = "\t", col.names = FALSE)

cat("\nSaved", nrow(coords), "regions to confirmed_targets.bed\n")

## also save the flat list of 21 target IDs for the primer-demux step
writeLines(unique(true_pairs$target_id), "confirmed_21_target_ids.txt")
cat("Saved target ID list to confirmed_21_target_ids.txt\n")
