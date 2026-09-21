
library(Biostrings)

results <- fread("dada2_bypass_check/top100_combined_results.csv")
#split primer pairs and search each individually through the caulobacter hits
top100_target_long <- results[, .(target_id = str_trim(str_split(matched_primer_targets, ";")[[1]])),
                              by = qseqid]
top100_unique_targets <- unique(top100_target_long$target_id)

confirmed_caulobacter_21 <- readLines("confirmed_21_target_ids.txt")

overlap <- intersect(top100_unique_targets, confirmed_caulobacter_21)

cat("Distinct primer targets represented in top 100:", length(top100_unique_targets), "\n")
cat("Of those, confirmed Caulobacter targets:", length(overlap), "\n")
print(overlap)

if (length(overlap) > 0) {
  cat("\nReads in the top 100 using a confirmed Caulobacter primer pair:\n")
  print(top100_target_long[target_id %in% overlap])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~Data for above section: 4 primer pairs, 7 of the top 38 pairs/100 reads.
# Distinct primer targets represented in top 100: 38 
# > cat("Of those, confirmed Caulobacter targets:", length(overlap), "\n")
# Of those, confirmed Caulobacter targets: 4 
# > print(overlap)
# [1] "RH.44DCEC943BED412Z0Z" "RH.0C6EF86AB29A422Z0Z" "RH.517A9FBA2B494D5Z0Z" "RH.04A0205BEA6C49EZ0Z"
# > if (length(overlap) > 0) {
#   +   cat("\nReads in the top 100 using a confirmed Caulobacter primer pair:\n")
#   +   print(top100_target_long[target_id %in% overlap])
#   + }
# 
# Reads in the top 100 using a confirmed Caulobacter primer pair:
#   qseqid             target_id
# <char>                <char>
#   1: M03520:605:000000000-LVNKY:1:1101:15415:1367 RH.44DCEC943BED412Z0Z
# 2: M03520:605:000000000-LVNKY:1:1101:10216:5547 RH.0C6EF86AB29A422Z0Z
# 3: M03520:605:000000000-LVNKY:1:1101:14854:2015 RH.517A9FBA2B494D5Z0Z
# 4: M03520:605:000000000-LVNKY:1:1101:17700:3522 RH.04A0205BEA6C49EZ0Z
# 5: M03520:605:000000000-LVNKY:1:1101:17582:2929 RH.0C6EF86AB29A422Z0Z
# 6: M03520:605:000000000-LVNKY:1:1101:21262:3305 RH.517A9FBA2B494D5Z0Z
# 7: M03520:605:000000000-LVNKY:1:1101:18037:2609 RH.517A9FBA2B494D5Z0Z
# > 

top100_source_map <- source_map[target_id %in% top100_unique_targets]
cat("Targets with a resolved source sequence:", nrow(top100_source_map), "out of", length(top100_unique_targets), "\n")

candidates <- readDNAStringSet("P450_final_targets_IDT.fasta")
candidate_ids <- names(candidates)

top100_seqs_for_blast <- candidates[candidate_ids %in% top100_source_map$true_source_id]

writeXStringSet(top100_seqs_for_blast, "dada2_bypass_check/top100_targets_for_ncbi_web.fasta")
cat("Wrote", length(top100_seqs_for_blast), "sequences for manual NCBI BLAST submission\n")

#Wrote 37 sequences for manual NCBI BLAST submission
