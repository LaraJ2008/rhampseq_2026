#!/bin/bash
## =============================================================================
## Per-sample off-target rate: what fraction of each sample's reads match
## NONE of the 380 panel primers?
## =============================================================================
## Unlike the top-100 check, this needs each sample dereplicated SEPARATELY --
## the earlier pooled-samples file lost per-sample identity. Using
## --minuniquesize 5 keeps this tractable on samples with millions of reads,
## while still capturing the large majority of total read mass (rare
## one-off sequencing-error variants contribute little to the total either way).
## =============================================================================

set -euo pipefail

FASTQ_DIR="/mnt/d/NSF_FAST_poc_extracted/p450_seq/fastq"
OUT_DIR="off_target_check"
MIN_SIZE=5

mkdir -p "$OUT_DIR"

if [ ! -f "all_primers_db.nsq" ]; then
  makeblastdb -in all_primers.fasta -dbtype nucl -out all_primers_db
fi

echo "sample,total_reads" > "$OUT_DIR/off_target_summary.csv"

for f in "$FASTQ_DIR"/*_R1_*.fastq.gz; do
  [[ "$f" == *Undetermined* ]] && continue

  sample=$(basename "$f" | sed 's/_S[0-9]*_L001_R1_001\.fastq\.gz$//')
  echo ""
  echo "=== Processing sample: $sample ==="

  total_reads=$(( $(zcat "$f" | wc -l) / 4 ))
  echo "Total reads: $total_reads"
  echo "${sample},${total_reads}" >> "$OUT_DIR/off_target_summary.csv"

  zcat "$f" | seqkit fq2fa -o "$OUT_DIR/${sample}.fasta"

  vsearch --derep_fulllength "$OUT_DIR/${sample}.fasta" \
          --sizeout \
          --minuniquesize $MIN_SIZE \
          --output "$OUT_DIR/${sample}_derep.fasta" 2>/dev/null

  n_unique=$(grep -c '^>' "$OUT_DIR/${sample}_derep.fasta" || echo 0)
  echo "Unique sequences kept (size>=$MIN_SIZE): $n_unique"

  if [ "$n_unique" -eq 0 ]; then
    echo "  no sequences passed the size cutoff, skipping BLAST for this sample"
    continue
  fi

  blastn -task blastn-short -query "$OUT_DIR/${sample}_derep.fasta" -db all_primers_db \
    -word_size 7 -evalue 1e-5 \
    -outfmt '6 qseqid sseqid pident length mismatch evalue bitscore' \
    -out "$OUT_DIR/${sample}_vs_primers.tsv"

  echo "  BLAST complete: $(wc -l < "$OUT_DIR/${sample}_vs_primers.tsv") raw hit rows"
done

echo ""
echo "=== DONE. Per-sample fasta/derep/blast files are in $OUT_DIR/ ==="
echo "Run the companion R script to compute final per-sample off-target rates."
