#!/bin/bash
## =============================================================================
## Skip DADA2: get top 100 most abundant raw sequences directly, no denoising
## =============================================================================
## Rationale: DADA2's error model and ASV inference are appropriate for
## single-locus amplicon data. Applying one uniform truncLen (chosen by eye
## from a quality plot) across a panel spanning ~150-290bp amplicons risks
## truncating longer targets incorrectly and conflating short-target reads.
## Exact-sequence dereplication (vsearch --derep_fulllength) makes no such
## assumption -- it just counts identical sequences as-is. Less powerful for
## detecting single-base sequencing errors, but transparent and defensible
## for a first-pass abundance/QC check across a heterogeneous panel.
## =============================================================================

set -euo pipefail

## ---- 0. CONFIG ---------------------------------------------------------------
FASTQ_DIR="/mnt/d/NSF_FAST_poc_extracted/p450_seq/fastq"   # adjust if needed
OUT_DIR="dada2_bypass_check"
TOP_N=100

mkdir -p "$OUT_DIR"

## ---- 1. Install vsearch, if not already present -------------------------------
if ! command -v vsearch &> /dev/null; then
  sudo apt install -y vsearch
fi

## ---- 2. Pool all R1 reads across all samples into one file --------------------
## Using R1 only for a quick global abundance ranking -- keeps this simple and
## avoids read-merging complications for the longer amplicons where R1/R2
## may not fully overlap. Every read still carries its primer sequence at the
## start, which is fine/expected -- we WANT primers present for the
## primer-matching step later.

echo "=== Pooling all R1 fastq files (excluding Undetermined) ==="
for f in "$FASTQ_DIR"/*_R1_*.fastq.gz; do
  [[ "$f" == *Undetermined* ]] && continue
  echo "  including: $(basename "$f")"
done
zcat $(for f in "$FASTQ_DIR"/*_R1_*.fastq.gz; do
         [[ "$f" == *Undetermined* ]] && continue
         echo "$f"
       done) > "$OUT_DIR/all_samples_R1.fastq"
echo "Total reads pooled: $(($(wc -l < "$OUT_DIR/all_samples_R1.fastq") / 4))"

## ---- 3. Convert to FASTA, then dereplicate: collapse identical sequences, ------
## count occurrences
## NOTE: vsearch --derep_fulllength requires FASTA input, not FASTQ -- quality
## scores aren't used in exact-match dereplication anyway, so nothing is lost
## by dropping them here.

echo "=== Converting pooled reads to FASTA ==="
seqkit fq2fa "$OUT_DIR/all_samples_R1.fastq" -o "$OUT_DIR/all_samples_R1.fasta"

echo "=== Dereplicating (exact match, no denoising) ==="
vsearch --derep_fulllength "$OUT_DIR/all_samples_R1.fasta" \
        --sizeout \
        --output "$OUT_DIR/derep_all.fasta" \
        --minuniquesize 1

echo "Unique sequences found: $(grep -c '^>' "$OUT_DIR/derep_all.fasta")"

## ---- 4. Sort by abundance, take top 100 ---------------------------------------
## vsearch's --sizeout already sorts by decreasing abundance by default when
## combined with --derep_fulllength, but --sortbysize makes it explicit/certain

vsearch --sortbysize "$OUT_DIR/derep_all.fasta" \
        --output "$OUT_DIR/derep_sorted.fasta"

seqkit head -n $TOP_N "$OUT_DIR/derep_sorted.fasta" > "$OUT_DIR/top100_abundant.fasta"

echo "Top $TOP_N sequences extracted:"
grep '^>' "$OUT_DIR/top100_abundant.fasta" | head -5
echo "..."

## ---- 5. BLAST top 100 against panel source/candidate target sequences ---------

echo ""
echo "=== BLASTing top 100 against panel target sequences ==="

if [ ! -f "P450_final_targets_IDT_db.nsq" ]; then
  makeblastdb -in P450_final_targets_IDT.fasta -dbtype nucl -out P450_final_targets_IDT_db
fi

blastn -query "$OUT_DIR/top100_abundant.fasta" -db P450_final_targets_IDT_db \
  -outfmt '6 qseqid sseqid pident length mismatch evalue bitscore qlen' \
  -evalue 1e-10 -max_target_seqs 3 \
  -out "$OUT_DIR/top100_vs_all_targets.tsv"

echo "Target matches found: $(wc -l < "$OUT_DIR/top100_vs_all_targets.tsv")"

## ---- 6. BLAST top 100 against ALL primers (deliberately, to check every one) --

echo ""
echo "=== BLASTing top 100 against panel primers ==="
echo "NOTE: all_primers.fasta currently contains only the 380 PRIMARY POOL"
echo "primers (760 sequences: 380 F + 380 R). The 35 SecondaryPool and 8"
echo "Singles primers are NOT included unless all_primers.fasta has been"
echo "regenerated with that data merged in."

if [ ! -f "all_primers_db.nsq" ]; then
  makeblastdb -in all_primers.fasta -dbtype nucl -out all_primers_db
fi

blastn -task blastn-short -query "$OUT_DIR/top100_abundant.fasta" -db all_primers_db \
  -word_size 7 -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore' \
  -out "$OUT_DIR/top100_vs_all_primers.tsv"

echo "Primer matches found: $(wc -l < "$OUT_DIR/top100_vs_all_primers.tsv")"

echo ""
echo "=== DONE ==="
echo "Files ready for R analysis in $OUT_DIR/:"
echo "  - top100_abundant.fasta          (the 100 sequences themselves, with abundance in header)"
echo "  - top100_vs_all_targets.tsv      (which panel target each read best matches)"
echo "  - top100_vs_all_primers.tsv      (every primer that matches each read, not just the intended one)"
