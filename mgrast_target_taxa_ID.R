## =============================================================================
## Reconstructing the P450 x Caulobacter candidate pool directly from MG-RAST
## =============================================================================
## Rationale: you don't have the IDT scientist's pre-panel candidate list, but
## MG-RAST itself already made both the functional call (cytochrome P450) and
## the organism call (Caulobacter) on your original contigs/genes. So instead
## of re-deriving taxonomy from scratch, we ask MG-RAST directly:
##   "of everything you annotated as cytochrome P450, which of those did you
##    ALSO annotate as Caulobacter?"
## That reconstructs the pre-panel candidate pool. Then we check whether any
## of those candidate sequences ended up as one of your final 415 targets.
##
## API docs:  https://help.mg-rast.org/api.html
## Base URL:  https://api.mg-rast.org/1/
## Requires: httr2, jsonlite, data.table, Biostrings
## =============================================================================

library(httr2)
library(jsonlite)
library(data.table)
library(Biostrings)

##--------00 MG-RAST DATAMINING VIA BASH -------------------------------------
#Used terminal through R to get analysis details
$ curl -X GET -H "auth: uep4GzARH35ZBjif2UJyS9JfN" "https://api.mg-rast.org/1/metagenome/mgm4962884.3"
{"pipeline_id":"34824161-0d3b-4ed9-b183-5c6918ec1e43",
  "created":"2022-03-18 15:44:17",
  "status":"private",
  "project":["mgp102073","https://api.mg-rast.org/project/mgp102073"],"job_id":530517,
  "owner":"mgu78486",
  "url":"https://api.mg-rast.org/metagenome/mgm4962884.3?verbosity=minimal",
  "pipeline_version":"4.0.3",
  "library":null,
  "sequence_type":"WGS",
  "id":"mgm4962884.3",
  "submission":"7e142e62-68df-45c4-a8e9-29e99cb8aa79",
  "md5_checksum":"ce3aad8168ab829c2ebe4dc365e0d087",
  "sample":null,
  "name":"DeLong_2_S10_L004_R",
  "version":1,
  "pipeline_parameters":{
      "dereplicate":"yes",
      "m5nr_sims_version":"7",
      "aa_pid":"90", ****************This is a very stringent match requirement to mgRAST template
      "assembled":"no", *****************Why not? 
      "file_type":"fastq",
      "min_qual":"15",
      "priority":"never",
      "max_lqb":"5",
      "m5rna_sims_version":"7",
      "fgs_type":"454", ***************This is pyroseq not illumina
      "bowtie":"yes",
      "m5nr_annotation_version":"1",
      "prefix_length":"50",
      "screen_indexes":"h_sapiens",
      "m5rna_annotation_version":"1",
      "dynamic_trim":"yes",
      "rna_pid":"97"}}


curl -m 60 -L -o local_p450_diagnosis/PF00067.hmm.gz \
"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00067?annotation=hmm"
gunzip -f local_p450_diagnosis/PF00067.hmm.gz


nohup hmmsearch --tblout local_p450_diagnosis/p450_hits.tbl \
--domtblout local_p450_diagnosis/p450_domhits.tbl \
-E 1e-3 \
local_p450_diagnosis/PF00067.hmm \
putative_coding_features_rRNA_filtered_S9.faa \
> local_p450_diagnosis/hmmsearch_full.out 2>&1 &

## ---- 0. CONFIG -----------------------------------------------------------

cfg <- list(
  auth_key   = "uep4GzARH35ZBjif2UJyS9JfN",   # from mg-rast.org my profile, show webkey
  mgm_ids    = c("mgm4962884.3", "mgm4962885.3"),         # your metagenome job ID(s); add all
                                                           # relevant jobs if data spans several
  api_base   = "https://api.mg-rast.org/1",

  # annotation-source choices matter here: SEED/Subsystems is the classic
  # functional source MG-RAST uses for pathway/subsystem-style function
  # names like "Cytochrome P450"; RefSeq is generally the more complete
  # organism-level source. Adjust if your job was annotated against
  # different default sources.
  function_source = "Subsystems",
  function_filter = "Cytochrome P450",
  organism_source = "RefSeq",
  organism_filter = "Caulobacter",
  evalue_cutoff    = 10,   # matches typical MG-RAST default; tighten if needed

  panel_csv  = "p450_panel_415_targets.csv",   # target_id, target_seq, ...
  out_dir    = "mgrast_diagnosis"
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Reachability / auth sanity check ----------------------------------

test_req <- request(cfg$api_base) |> req_url_path_append("metagenome", cfg$mgm_ids[1])
test_resp <- tryCatch(req_perform(test_req), error = function(e) e)

if (inherits(test_resp, "error")) {
  stop("Could not reach MG-RAST API or the job ID/auth_key is wrong. ",
       "Error: ", conditionMessage(test_resp),
       "\nTry the same URL in a browser first: ",
       cfg$api_base, "/metagenome/", cfg$mgm_ids[1])
}
cat("MG-RAST API reachable, job", cfg$mgm_ids[1], "found.\n")

## ---- 2. Helper: paged annotation/sequence query ---------------------------
## The annotation/sequence resource returns matched reads/features as a
## tab-delimited block (md5, sequence, semicolon-list of annotations, seq id)
## rather than JSON -- handle both possible response shapes defensively.

query_annotation <- function(mgm_id, type, source, filter, evalue, auth_key,
                              api_base, limit = 1000) {
  req <- request(api_base) |>
    req_url_path_append("annotation", "sequence", mgm_id) |>
    req_url_query(type = type, source = source, filter = filter,
                   evalue = evalue, auth_key = auth_key, limit = limit)

  resp <- req_perform(req)
  body <- resp_body_string(resp)

  # Response is typically tab-delimited: md5 <tab> sequence <tab> annotations <tab> seq_id
  if (grepl("^\\{", trimws(body))) {
    # got JSON (e.g. an error or a wrapped object) -- surface it
    parsed <- fromJSON(body)
    if (!is.null(parsed$ERROR)) {
      warning("MG-RAST API error for ", mgm_id, " / ", type, "=", filter, ": ",
              parsed$ERROR)
      return(data.table())
    }
  }

  lines <- strsplit(body, "\n")[[1]]
  lines <- lines[nzchar(lines)]
  parts <- tstrsplit(lines, "\t", fixed = TRUE)
  if (length(parts) < 4) return(data.table())

  dt <- data.table(md5 = parts[[1]], sequence = parts[[2]],
                    annotation = parts[[3]], seq_id = parts[[4]],
                    mgm_id = mgm_id, query_filter = filter)
  dt
}

## ---- 3. Pull P450-annotated features and Caulobacter-annotated features ---

p450_hits <- rbindlist(lapply(cfg$mgm_ids, function(id) {
  cat("Querying P450 function annotations for", id, "...\n")
  query_annotation(id, type = "function", source = cfg$function_source,
                    filter = cfg$function_filter, evalue = cfg$evalue_cutoff,
                    auth_key = cfg$auth_key, api_base = cfg$api_base)
}))
fwrite(p450_hits, file.path(cfg$out_dir, "mgrast_p450_function_hits.csv"))
cat("P450-annotated features pulled from MG-RAST:", nrow(p450_hits), "\n")

caulo_hits <- rbindlist(lapply(cfg$mgm_ids, function(id) {
  cat("Querying Caulobacter organism annotations for", id, "...\n")
  query_annotation(id, type = "organism", source = cfg$organism_source,
                    filter = cfg$organism_filter, evalue = cfg$evalue_cutoff,
                    auth_key = cfg$auth_key, api_base = cfg$api_base)
}))
fwrite(caulo_hits, file.path(cfg$out_dir, "mgrast_caulobacter_organism_hits.csv"))
cat("Caulobacter-annotated features pulled from MG-RAST:", nrow(caulo_hits), "\n")

## ---- 4. Intersect: candidate pool of P450 x Caulobacter features ---------
## Join on md5 (the M5nr identity checksum) which is consistent across
## annotation types for the same underlying sequence/feature.

candidate_pool <- merge(
  p450_hits[, .(md5, sequence, mgm_id, p450_annotation = annotation)],
  caulo_hits[, .(md5, caulo_annotation = annotation)],
  by = "md5"
)
fwrite(candidate_pool, file.path(cfg$out_dir, "candidate_p450_caulobacter_pool.csv"))

cat("\n=== Candidate pool: features called BOTH P450 AND Caulobacter ===\n")
cat("N =", nrow(candidate_pool), "\n")

if (nrow(candidate_pool) == 0) {
  cat("\n*** MG-RAST itself never jointly annotated any feature as both\n",
      "cytochrome P450 AND Caulobacter. This pushes the explanation upstream\n",
      "of panel design entirely -- into MG-RAST's own annotation/assembly\n",
      "step. Possible reasons: Caulobacter contigs assembled but their P450\n",
      "genes fell below the evalue/similarity cutoff for a function call, a\n",
      "chimeric-assembly issue merged Caulobacter sequence with another\n",
      "organism's annotation, or the two annotation SOURCES used here\n",
      "(Subsystems / RefSeq) simply don't have good Caulobacter P450\n",
      "coverage -- worth re-running this query with source='SEED' or\n",
      "source='KEGG' for function, and source='SILVA'/'GenBank' for organism,\n",
      "since MG-RAST computes annotations against multiple source databases\n",
      "independently and results can differ meaningfully between them.\n")
}

## ---- 5. Check whether any candidate sequence made it into the final panel -

if (nrow(candidate_pool) > 0 && file.exists(cfg$panel_csv)) {
  panel <- fread(cfg$panel_csv)

  # Exact/near-exact match: does any panel target_seq contain (or match) a
  # candidate sequence? Panel targets are presumably shorter amplicon-region
  # subsequences of the original MG-RAST contig, so containment is the
  # right test, not full-length equality.
  candidate_seqs <- DNAStringSet(candidate_pool$sequence)
  names(candidate_seqs) <- candidate_pool$md5

  panel_seqs <- DNAStringSet(panel$target_seq)
  names(panel_seqs) <- panel$target_id

  match_results <- rbindlist(lapply(seq_along(candidate_seqs), function(i) {
    hits <- vcountPattern(as.character(candidate_seqs[[i]]), panel_seqs,
                           max.mismatch = 5)
    matched_targets <- panel$target_id[hits > 0]
    if (length(matched_targets) == 0) {
      data.table(candidate_md5 = names(candidate_seqs)[i], matched_target_id = NA)
    } else {
      data.table(candidate_md5 = names(candidate_seqs)[i],
                 matched_target_id = matched_targets)
    }
  }))

  fwrite(match_results, file.path(cfg$out_dir, "candidate_vs_panel_match.csv"))
  n_survived <- sum(!is.na(match_results$matched_target_id))

  cat("\n=== Did any Caulobacter P450 candidate survive into the final panel? ===\n")
  cat("Candidates found in MG-RAST:", nrow(candidate_pool), "\n")
  cat("Candidates matched to a final target_id:", n_survived, "\n")

  if (n_survived == 0) {
    cat("\n*** Caulobacter P450 candidates existed in MG-RAST's own annotation,\n",
        "but NONE of them ended up as one of the 415 finalized targets. This\n",
        "points to the panel-narrowing/multiplex-QC step (run by IDT) as the\n",
        "dropout point -- worth asking IDT directly whether their multiplex\n",
        "compatibility screen (GC%, secondary structure, primer-dimer risk)\n",
        "rejected these specific candidate sequences. Caulobacter's high GC\n",
        "content (~67%) is a plausible mechanistic reason.\n")
  } else {
    cat("\n*** At least one Caulobacter-derived P450 target did make it into\n",
        "the final panel. Now go back to the rhampseq_panel_diagnosis.R script\n",
        "and check demux/direct-mapping read counts specifically for these\n",
        "matched_target_id values -- this tells you whether it's a true wet-lab\n",
        "amplification failure for those specific loci.\n")
    print(match_results[!is.na(matched_target_id)])
  }
} else if (nrow(candidate_pool) > 0) {
  cat("\nCandidate pool exists in MG-RAST but panel_csv wasn't found -- point\n",
      "cfg$panel_csv at your finalized target list and re-run this last step\n",
      "to check for survival into the panel.\n")
}
