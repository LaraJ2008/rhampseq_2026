#!/bin/bash
## =============================================================================
## Final targeted check: do the 21 confirmed Caulobacter-derived targets show
## ANY evidence in the real targeted fastq data?
## =============================================================================
## Two independent tests, run in parallel:
##   A) Primer-based demux -- did the EXACT designed primer pair for these 21
##      targets recognize and trim any real reads? (tests wet-lab amplification
##      specifically at these loci)
##   B) Direct genomic pileup -- do ANY reads, regardless of primer matching,
##      align to the exact genomic coordinates these 21 targets came from?
##      (tests whether the sequence exists in the data at all, independent of
##      primer-level demultiplexing -- catches cases where amplification
##      happened but got mis-assigned)
##
## Run from ~/p450_stringent_analysis_2026 after running
## get_confirmed_target_regions.R (which produces confirmed_targets.bed and
## confirmed_21_target_ids.txt)
## =============================================================================

set -euo pipefail

## ---- CONFIG: point these at your real fastq files ---------------------------
#FASTQ_R1="CHANGE_ME_R1.fastq.gz"
#FASTQ_R2="CHANGE_ME_R2.fastq.gz"
#SAMPLE_NAME="S9"   # or whatever labels your sample

FASTQ_R1="/mnt/d/NSF_FAST_poc_extracted/p450_seq/fastq/as3-1_S1_L001_R1_001.fastq.gz"
FASTQ_R2="/mnt/d/NSF_FAST_poc_extracted/p450_seq/fastq/as3-1_S1_L001_R2_001.fastq.gz"
SAMPLE_NAME="as3-1"

mkdir -p final_check

## -----------------------------------------------------------------------------
## PART A -- primer-based demux, just these 21 target pairs
## -----------------------------------------------------------------------------

echo "=== Extracting primer sequences for the 21 confirmed targets ==="
grep -f <(sed 's/$/_F/' confirmed_21_target_ids.txt) all_primers.fasta -A1 --no-group-separator > final_check/confirmed_primers_F.fasta
grep -f <(sed 's/$/_R/' confirmed_21_target_ids.txt) all_primers.fasta -A1 --no-group-separator > final_check/confirmed_primers_R.fasta

echo "F primers extracted: $(grep -c '^>' final_check/confirmed_primers_F.fasta)"
echo "R primers extracted: $(grep -c '^>' final_check/confirmed_primers_R.fasta)"

echo "=== Running cutadapt demux (21 target pairs only) ==="
cutadapt -g file:final_check/confirmed_primers_F.fasta -G file:final_check/confirmed_primers_R.fasta -e 0.1 --no-indels --pair-adapters -o "final_check/${SAMPLE_NAME}_{name}_R1.fastq.gz" -p "final_check/${SAMPLE_NAME}_{name}_R2.fastq.gz" --json="final_check/${SAMPLE_NAME}_cutadapt.json" "$FASTQ_R1" "$FASTQ_R2"

echo ""
echo "=== Per-target demux read counts ==="
for f in final_check/${SAMPLE_NAME}_RH.*_R1.fastq.gz; do
[ -e "$f" ] || continue
tid=$(basename "$f" | sed "s/^${SAMPLE_NAME}_//; s/_R1\.fastq\.gz$//")
n=$(zcat "$f" | wc -l)
n=$((n / 4))
echo -e "${tid}\t${n}"
done | tee final_check/demux_counts_21_targets.tsv

## -----------------------------------------------------------------------------
## PART B -- direct genomic pileup at the 21 confirmed loci, primer-agnostic
## -----------------------------------------------------------------------------

echo ""
echo "=== Aligning full fastq to Caulobacter genomes (bwa mem) ==="
bwa mem -t 4 caulobacter_genomes.fasta "$FASTQ_R1" "$FASTQ_R2" \
| samtools sort -@ 4 -o final_check/reads_vs_caulobacter.sorted.bam -
  samtools index final_check/reads_vs_caulobacter.sorted.bam

echo ""
echo "=== Read counts at each of the 21 confirmed genomic loci (primer-agnostic) ==="
while IFS=$'\t' read -r chrom start end target_id; do
count=$(samtools view -c final_check/reads_vs_caulobacter.sorted.bam "${chrom}:${start}-${end}")
echo -e "${target_id}\t${chrom}:${start}-${end}\t${count}"
done < confirmed_targets.bed | tee final_check/pileup_counts_21_targets.tsv

echo ""
echo "=== DONE. Compare final_check/demux_counts_21_targets.tsv (primer-specific)"
echo "    against final_check/pileup_counts_21_targets.tsv (primer-agnostic)"
echo ""
echo "Interpretation:"
echo "  - Both near-zero for a target      -> genuine amplification failure"
echo "  - Pileup > 0 but demux ~0          -> reads exist but primer wasn't"
echo "                                          recognizing them (mismatch/"
echo "                                          demux parameter issue)"
echo "  - Both > 0                          -> Caulobacter reads genuinely"
echo "                                          present; original 'zero hits'"
echo "                                          conclusion was a classification"
echo "                                          problem, not a data problem"