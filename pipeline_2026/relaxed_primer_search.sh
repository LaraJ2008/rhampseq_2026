#!/bin/bash
set -euo pipefail

## Relaxed primer search: same database, but we deliberately do NOT filter
## by the strict length>=18/pident>=90 threshold used everywhere else tonight.
## -evalue 1000 keeps the search maximally sensitive; filtering happens in R,
## where we can look at the FULL distribution of match quality rather than a
## hard cutoff.

blastn -task blastn-short -query /mnt/d/NSF_FAST_poc_extracted/unexplained_chimera_for_uniprot.fasta -db all_primers_db -word_size 6 -evalue 1000 -outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore' -out unexamined_check/unresolved_relaxed_primer_hits.tsv

wc -l unexamined_check/unresolved_relaxed_primer_hits.tsv
