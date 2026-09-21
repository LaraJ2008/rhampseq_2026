## =============================================================================
## From-scratch rhAmpSeq multiplex P450 panel analysis
## =============================================================================
## Context: primers were designed purely from MG-RAST P450-annotated contigs,
## with NO taxonomic consideration. So we don't yet know which of the 415
## targets are Caulobacter-derived -- that has to be established first (Part 1)
## before we can ask why Caulobacter reads are missing from the targeted data
## (Parts 2-4). Part 5 tests whether Caulobacter is actually an outlier or
## just the visible tip of a broader panel dropout pattern.
##
## Requires (one-time):
##   BiocManager::install(c("Rsubread","Rsamtools","GenomicAlignments",
##                           "Biostrings","ShortRead"))
##   install.packages(c("data.table","rentrez","ggplot2"))
##   cutadapt on PATH (pip install cutadapt --break-system-packages, or conda)
## =============================================================================

library(Rsubread)
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(ShortRead)
library(data.table)
library(rentrez)

## ---- 0. CONFIG ----------------------------------------------------------------

cfg <- list(
  # Finalized panel metadata you already have: one row per target, both
  # sub-panels (380 + 35) combined, with primer sequences and the
  # amplicon/target sequence pulled from the source MG-RAST contig.
  panel_csv        = "p450_panel_415_targets.csv",
  # required columns: target_id, primer_F, primer_R, target_seq
  # optional if you have it: source_contig_id, mgrast_job_id, mgrast_annotation

  fastq_dir        = "raw_fastq",
  fastq_R1_pattern = "_R1.fastq.gz",
  fastq_R2_pattern = "_R2.fastq.gz",
  paired_end       = TRUE,

  out_dir          = "rhampseq_diagnosis",
  cutadapt_bin     = "cutadapt",

  min_pident       = 95,   # nucleotide-level, rhAmp is high fidelity
  min_qcov         = 90
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
panel <- fread(cfg$panel_csv)
stopifnot(all(c("target_id","primer_F","primer_R","target_seq") %in% names(panel)))
cat("Panel size:", nrow(panel), "targets\n")

## =============================================================================
## PART 1 -- Assign taxonomy to each of the 415 targets
## =============================================================================
## You don't have this yet, since design was taxonomy-blind. Two options:
##  (a) PREFERRED if you still have MG-RAST job access: pull the taxonomy call
##      MG-RAST already made for each source contig/gene via its API (matR
##      package) -- this is consistent with how the panel was built and avoids
##      a second, possibly-disagreeing classification step.
##  (b) FALLBACK used here: BLAST each target_seq against nt via NCBI's remote
##      blastn and pull taxonomy for the best hit via rentrez. Fine for ~415
##      sequences; slow but only needs to be run once, so results are cached.

taxonomy_cache_file <- file.path(cfg$out_dir, "target_taxonomy.csv")

if (!file.exists(taxonomy_cache_file)) {

  target_fasta <- file.path(cfg$out_dir, "targets.fasta")
  writeXStringSet(
    DNAStringSet(setNames(panel$target_seq, panel$target_id)),
    filepath = target_fasta
  )

  blast_tsv <- file.path(cfg$out_dir, "targets_vs_nt.tsv")
  # -remote avoids needing a local nt copy; for 415 short sequences this is
  # workable but can take a while and is rate-limited. Batch if it stalls.
  cmd <- sprintf(
    "blastn -query %s -db nt -remote -outfmt '6 qseqid sseqid pident length qlen evalue bitscore staxid' -max_target_seqs 3 > %s",
    target_fasta, blast_tsv
  )
  system(cmd)

  blast_res <- fread(blast_tsv, col.names = c("target_id","sseqid","pident",
                        "length","qlen","evalue","bitscore","staxid"))
  blast_res[, qcov := 100 * length / qlen]
  blast_res <- blast_res[pident >= cfg$min_pident & qcov >= cfg$min_qcov]
  best_hit <- blast_res[order(-bitscore)][, .SD[1], by = target_id]

  # Resolve taxid -> genus/species name via rentrez (batched, be gentle on NCBI)
  uniq_taxids <- unique(best_hit$staxid)
  tax_lookup <- rbindlist(lapply(uniq_taxids, function(tid) {
    Sys.sleep(0.34)  # stay under NCBI's rate limit without an API key
    s <- tryCatch(entrez_fetch(db = "taxonomy", id = tid, rettype = "xml",
                                parsed = TRUE), error = function(e) NULL)
    if (is.null(s)) return(data.table(staxid = tid, genus = NA, species = NA))
    rec <- XML::xmlToList(s)$Taxon
    data.table(staxid = tid,
               genus  = rec$LineageEx[sapply(rec$LineageEx,
                          function(x) x$Rank == "genus")][[1]]$ScientificName %||% NA,
               species = rec$ScientificName %||% NA)
  }), fill = TRUE)

  target_taxonomy <- merge(best_hit[, .(target_id, pident, qcov, staxid)],
                            tax_lookup, by = "staxid", all.x = TRUE)
  target_taxonomy <- merge(panel[, .(target_id)], target_taxonomy,
                            by = "target_id", all.x = TRUE)
  fwrite(target_taxonomy, taxonomy_cache_file)
} else {
  target_taxonomy <- fread(taxonomy_cache_file)
}

n_caulo_targets <- target_taxonomy[genus == "Caulobacter", .N]
cat("\nTargets taxonomically assigned to Caulobacter:", n_caulo_targets, "of", nrow(panel), "\n")

if (n_caulo_targets == 0) {
  cat("\n*** No target in the finalized panel traces back to Caulobacter. ***\n",
      "This means the dropout happened BEFORE any wet-lab step -- either in\n",
      "MG-RAST's contig binning/annotation, or in the panel-narrowing step from\n",
      "the full P450 candidate list down to 380+35. Check design-stage records\n",
      "before doing any further fastq analysis: this fully explains a zero-hit\n",
      "result and no primer troubleshooting is needed.\n")
}

## =============================================================================
## PART 2 -- Primer-based demultiplexing of raw reads (which targets actually
## got amplified, by design)
## =============================================================================
## cutadapt with linked adapters (forward primer 5', reverse primer's
## revcomp expected downstream) lets us assign each read pair to its intended
## target and get a genuine "did this specific primer pair produce product"
## count, distinct from post-hoc sequence similarity.

r1_files <- list.files(cfg$fastq_dir, pattern = cfg$fastq_R1_pattern, full.names = TRUE)
r2_files <- sub(cfg$fastq_R1_pattern, cfg$fastq_R2_pattern, r1_files)

demux_dir <- file.path(cfg$out_dir, "demux")
dir.create(demux_dir, showWarnings = FALSE)

# Build a cutadapt multi-primer fasta pair (name=target_id) once
primers_fwd_fa <- file.path(cfg$out_dir, "primers_fwd.fasta")
primers_rev_fa <- file.path(cfg$out_dir, "primers_rev.fasta")
writeXStringSet(DNAStringSet(setNames(panel$primer_F, panel$target_id)), primers_fwd_fa)
writeXStringSet(DNAStringSet(setNames(panel$primer_R, panel$target_id)), primers_rev_fa)

demux_counts <- list()
for (i in seq_along(r1_files)) {
  sample_name <- sub(cfg$fastq_R1_pattern, "", basename(r1_files[i]))
  out_r1 <- file.path(demux_dir, paste0(sample_name, "_{name}_R1.fastq.gz"))
  out_r2 <- file.path(demux_dir, paste0(sample_name, "_{name}_R2.fastq.gz"))
  log_file <- file.path(demux_dir, paste0(sample_name, "_cutadapt.json"))

  # -g / -G = 5' adapters (forward/reverse primers), allow a couple of
  # mismatches (-e), require full primer length overlap (--no-indels keeps
  # this a clean primer-identity check rather than a fuzzy alignment)
  cmd <- sprintf(
    "%s -g file:%s -G file:%s -e 0.1 --no-indels --pair-adapters -o %s -p %s --json=%s %s %s",
    cfg$cutadapt_bin, primers_fwd_fa, primers_rev_fa, out_r1, out_r2,
    log_file, r1_files[i], r2_files[i]
  )
  system(cmd)

  # tally trimmed read counts per target_id from the output filenames
  out_files <- list.files(demux_dir, pattern = paste0("^", sample_name, "_.*_R1\\.fastq\\.gz$"),
                           full.names = TRUE)
  counts <- rbindlist(lapply(out_files, function(f) {
    tid <- sub(paste0("^", sample_name, "_"), "", sub("_R1\\.fastq\\.gz$", "", basename(f)))
    n <- length(ShortRead::readFastq(f))
    data.table(sample = sample_name, target_id = tid, demux_reads = n)
  }))
  demux_counts[[sample_name]] <- counts
}

demux_summary <- rbindlist(demux_counts)
fwrite(demux_summary, file.path(cfg$out_dir, "demux_per_target_counts.csv"))

## =============================================================================
## PART 3 -- Direct reference mapping (primer-agnostic check)
## =============================================================================
## This catches reads that match a target's sequence even if cutadapt didn't
## assign them there via primer match -- e.g. multiplex cross-talk, primer
## dimers, or a read whose primer region has enough mismatch to fail
## demultiplexing but whose amplicon body still matches the intended target.

index_dir <- file.path(cfg$out_dir, "ref_index")
dir.create(index_dir, showWarnings = FALSE)
target_fasta <- file.path(cfg$out_dir, "targets.fasta")
index_prefix <- file.path(index_dir, "targets_idx")
if (!file.exists(paste0(index_prefix, ".reads"))) {
  buildindex(basename = index_prefix, reference = target_fasta)
}

map_dir <- file.path(cfg$out_dir, "direct_map")
dir.create(map_dir, showWarnings = FALSE)
bam_out <- file.path(map_dir, sub(cfg$fastq_R1_pattern, ".bam", basename(r1_files)))

align(index = index_prefix, readfile1 = r1_files, readfile2 = r2_files,
      input_format = "gzFASTQ", output_file = bam_out,
      unique = FALSE, nBestLocations = 5, nthreads = 4)

direct_map_counts <- rbindlist(lapply(bam_out, function(b) {
  param <- ScanBamParam(what = c("qname","rname"),
                         flag = scanBamFlag(isUnmappedQuery = FALSE))
  aln <- as.data.table(scanBam(b, param = param)[[1]])
  aln[, sample := sub("\\.bam$", "", basename(b))]
  aln[, .N, by = .(sample, target_id = rname)]
}))
fwrite(direct_map_counts, file.path(cfg$out_dir, "direct_mapping_per_target_counts.csv"))

## =============================================================================
## PART 4 -- Combine everything per target, flag Caulobacter specifically
## =============================================================================

combined <- merge(demux_summary, direct_map_counts,
                   by = c("sample","target_id"), all = TRUE)
combined[is.na(demux_reads), demux_reads := 0]
combined[is.na(N), N := 0]
setnames(combined, "N", "direct_map_reads")
combined <- merge(combined, target_taxonomy[, .(target_id, genus)],
                   by = "target_id", all.x = TRUE)

fwrite(combined, file.path(cfg$out_dir, "combined_per_target_per_sample.csv"))

cat("\n=== Caulobacter-derived targets: primer-demux vs direct-mapping read counts ===\n")
print(combined[genus == "Caulobacter"])

cat("
Interpretation of the Caulobacter rows above:
  - both demux_reads and direct_map_reads are ~0
      -> no evidence of product at all; real amplification failure. Now do a
         nucleotide-level primer-vs-template mismatch check for these
         specific target_ids (Part 6 note below).
  - demux_reads ~0 but direct_map_reads > 0
      -> product exists in the data but isn't being correctly identified by
         its intended primer pair -- e.g. primer-region mismatches large
         enough to fail cutadapt's primer match but not enough to prevent
         some other mechanism producing amplicon-matching sequence, or
         cross-amplification from a different primer pair in the pool that
         happens to produce a similar product. Worth a manual look at a few
         reads.
  - both > 0
      -> Caulobacter is NOT actually missing; whatever produced the original
         'zero hits' conclusion (the by-hand top-100 BLAST) was the problem,
         not the panel or the sequencing.
\n")

## =============================================================================
## PART 5 -- Is Caulobacter actually an outlier, or part of a broader dropout?
## =============================================================================
## Before concluding anything Caulobacter-specific, check whether many
## targets across many genera show zero reads -- if dropout is widespread,
## the Caulobacter case is just the most visible instance of a systemic
## issue (e.g. panel-wide GC bias, uneven primer Tm, pooling imbalance)
## rather than something particular to Caulobacter primers or template.

per_target_total <- combined[, .(total_demux = sum(demux_reads),
                                  total_direct = sum(direct_map_reads)),
                              by = target_id]
per_target_total <- merge(per_target_total, target_taxonomy[, .(target_id, genus)],
                           by = "target_id", all.x = TRUE)

dropout_by_genus <- per_target_total[, .(
  n_targets = .N,
  n_zero    = sum(total_demux == 0),
  pct_zero  = 100 * sum(total_demux == 0) / .N
), by = genus][order(-pct_zero)]

fwrite(dropout_by_genus, file.path(cfg$out_dir, "dropout_rate_by_genus.csv"))
cat("\n=== Dropout rate by genus (targets with zero demuxed reads) ===\n")
print(dropout_by_genus)

cat("\nIf Caulobacter's dropout rate sits well above most other genera, that's\n",
    "a genuine Caulobacter-specific effect worth chasing (GC content, primer\n",
    "design filtering, rhAmp chemistry at those loci). If it's in line with a\n",
    "broadly high dropout rate across the whole panel, focus on panel-wide QC\n",
    "(primer Tm/GC distribution, pooling balance, PCR conditions) instead.\n")

## =============================================================================
## PART 6 -- (only if Part 4 showed a true dropout) primer-vs-template
## nucleotide mismatch check for the Caulobacter target_ids
## =============================================================================

check_primer_mismatches <- function(target_id, panel) {
  row <- panel[target_id == target_id]
  fwd <- DNAString(row$primer_F)
  rev <- DNAString(row$primer_R)
  amp <- DNAString(row$target_seq)

  fwd_hit <- matchPattern(fwd, amp, max.mismatch = 5, with.indels = TRUE)
  rev_hit <- matchPattern(reverseComplement(rev), amp, max.mismatch = 5, with.indels = TRUE)

  list(target_id = target_id,
       fwd_matches_own_template = length(fwd_hit) > 0,
       rev_matches_own_template = length(rev_hit) > 0)
}

# Sanity check: do the designed primers even match their OWN source
# sequence perfectly? (Should be trivially yes -- if not, something is off
# in how target_seq/primer sequences were paired in your panel table.)
caulo_ids <- target_taxonomy[genus == "Caulobacter", target_id]
if (length(caulo_ids) > 0) {
  self_check <- rbindlist(lapply(caulo_ids, check_primer_mismatches, panel = panel))
  fwrite(self_check, file.path(cfg$out_dir, "caulobacter_primer_self_check.csv"))
  print(self_check)
}
