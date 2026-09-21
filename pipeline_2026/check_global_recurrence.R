library(Biostrings)
library(data.table)
library(stringr)

## Uses derep_sorted.fasta -- the GLOBAL pooled-across-all-samples dereplicated
## set (from the top-100/accumulation-curve work) -- to check total abundance
## of each of the 91 unresolved sequences across the ENTIRE dataset, not just
## the one sample each happened to be sampled from.

global_pool <- readDNAStringSet("dada2_bypass_check/derep_sorted.fasta")
global_ids <- str_remove(names(global_pool), ";size=\\d+$")
global_abundance <- as.integer(str_extract(names(global_pool), "(?<=size=)\\d+"))

unresolved <- readDNAStringSet("/mnt/d/NSF_FAST_poc_extracted/unexplained_chimera_for_uniprot.fasta")
unresolved_ids <- str_remove(names(unresolved), ";size=\\d+$")

lookup <- data.table(qseqid = unresolved_ids)
lookup[, global_abundance := global_abundance[match(qseqid, global_ids)]]

cat("=== Global recurrence of the 91 unresolved sequences ===\n")
cat("(this is total copies across ALL pooled samples, not just their own sample)\n\n")
print(lookup[order(-global_abundance)])

cat("\nSummary:\n")
print(summary(lookup$global_abundance))

fwrite(lookup, "unresolved_reads_global_abundance.csv")

cat("\nIf several of these show high global abundance (well above the 1-4",
    "seen in their own sample), that's evidence of a real, recurring sequence\n")
cat("appearing across many samples -- not isolated noise.\n")
