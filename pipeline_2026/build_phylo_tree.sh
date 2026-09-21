#!/bin/bash
set -euo pipefail

sudo apt install -y mafft fasttree

## ---- 1. Use the combined fasta already produced by R Part 1 -----------------
## (targets + chimeras + single-primed reads, tagged in tree_input_metadata.csv)

cp tree_input_all_groups.fasta tree_targets.fasta
grep -c '^>' tree_targets.fasta

## ---- 2. Align with MAFFT ------------------------------------------------------

mafft --auto tree_targets.fasta > tree_targets_aligned.fasta

## ---- 3. Build tree with FastTree ---------------------------------------------
## -nt for nucleotide; FastTree reports local support values (similar in
## spirit to bootstrap, much faster to compute) at internal nodes

fasttree -nt -gtr tree_targets_aligned.fasta > tree_targets.nwk

echo "Tree written to tree_targets.nwk -- move to R for visualization"
