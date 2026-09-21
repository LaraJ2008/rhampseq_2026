#!/bin/bash
## =============================================================================
## What are the "unexamined" reads? Isolate sequences that occurred <5 times
## (excluded by --minuniquesize 5 in the original off-target check) and test
## whether they still match our primers.
## =============================================================================
## Rationale: a real, rare-but-genuine P450 variant -- e.g. a related gene or
## allele that a primer happens to bind despite different "internal" sequence
## -- would still show a clean primer match at the read's ends, even though
## it never reached abundance 5+ as an EXACT duplicate. This directly tests
## that hypothesis, separate from whether it's just sequencing-error noise.
## =============================================================================

set -euo pipefail

OUT_DIR="off_target_check"
UNEXAMINED_DIR="unexamined_check"
mkdir -p "$UNEXAMINED_DIR"

if [ ! -f "all_primers_db.nsq" ]; then
  makeblastdb -in all_primers.fasta -dbtype nucl -out all_primers_db
fi

for fasta in "$OUT_DIR"/*.fasta; do
  # skip the derep outputs themselves, we want the ORIGINAL per-sample fasta
  [[ "$fasta" == *_derep.fasta ]] && continue

  sample=$(basename "$fasta" .fasta)
  echo ""
  echo "=== Isolating unexamined (size 1-4) sequences: $sample ==="

  vsearch --derep_fulllength "$fasta" \
          --sizeout \
          --maxuniquesize 4 \
          --output "$UNEXAMINED_DIR/${sample}_unexamined.fasta" 2>/dev/null

  n_unexamined_uniques=$(grep -c '^>' "$UNEXAMINED_DIR/${sample}_unexamined.fasta" || echo 0)
  echo "Unique sequences with size 1-4: $n_unexamined_uniques"

  if [ "$n_unexamined_uniques" -eq 0 ]; then
    echo "  none found, skipping BLAST"
    continue
  fi

  blastn -task blastn-short -query "$UNEXAMINED_DIR/${sample}_unexamined.fasta" -db all_primers_db \
    -word_size 7 -evalue 1e-5 \
    -outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore' \
    -out "$UNEXAMINED_DIR/${sample}_unexamined_vs_primers.tsv"

  echo "  BLAST complete: $(wc -l < "$UNEXAMINED_DIR/${sample}_unexamined_vs_primers.tsv") raw hit rows"
done

echo ""
echo "=== DONE. Run the companion R script to summarize matched vs unmatched"
echo "    among the previously-excluded low-abundance reads. ==="
