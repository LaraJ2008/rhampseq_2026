install.packages(c("rentrez", "XML"))
library(data.table)
library(rentrez)
library(XML)
library(stringr)

## =============================================================================
## Get full taxonomic lineage (family/order/class) for each best-hit taxid
## =============================================================================
## Same lookup approach used early in this investigation for MG-RAST/Caulobacter
## taxonomy resolution. Unlike a live BLAST search, this is a lightweight
## per-ID XML fetch against NCBI's taxonomy database -- small batch size here
## should be quick and reliable.

best_hit <- fread("best_hit_per_target_derived.csv")

unique_taxids <- unique(best_hit$taxid)
cat("Unique taxids to look up:", length(unique_taxids), "\n")

## ---- Fetch full lineage for each taxid ----------------------------------------

get_lineage <- function(taxid) {
  Sys.sleep(0.34)  # stay under NCBI's rate limit without an API key
  res <- tryCatch(
    entrez_fetch(db = "taxonomy", id = taxid, rettype = "xml", parsed = TRUE),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(data.table(taxid = taxid, family = NA, order = NA, class = NA,
                       phylum = NA, full_sciname = NA))
  }

  rec <- tryCatch(xmlToList(res)$Taxon, error = function(e) NULL)
  if (is.null(rec)) {
    return(data.table(taxid = taxid, family = NA, order = NA, class = NA,
                       phylum = NA, full_sciname = NA))
  }

  lineage <- rec$LineageEx
  get_rank <- function(rank_name) {
    hit <- Filter(function(x) !is.null(x$Rank) && x$Rank == rank_name, lineage)
    if (length(hit) > 0) hit[[1]]$ScientificName else NA
  }

  data.table(
    taxid = taxid,
    family = get_rank("family"),
    order  = get_rank("order"),
    class  = get_rank("class"),
    phylum = get_rank("phylum"),
    full_sciname = rec$ScientificName
  )
}

lineage_table <- rbindlist(lapply(unique_taxids, get_lineage), fill = TRUE)

fwrite(lineage_table, "taxid_lineage_lookup.csv")
cat("Saved lineage lookup for", nrow(lineage_table), "taxids\n\n")

## ---- Merge back into best_hit table -------------------------------------------

best_hit_with_lineage <- merge(best_hit, lineage_table, by = "taxid", all.x = TRUE)
fwrite(best_hit_with_lineage, "best_hit_per_target_with_lineage.csv")

## ---- Specifically look at the "uncultured" entries ---------------------------

uncultured <- best_hit_with_lineage[str_detect(sciname, "(?i)uncultured")]

cat("=== 'Uncultured' entries -- higher-level taxonomy ===\n")
print(uncultured[, .(qseqid, sciname, family, order, class, phylum)])

cat("\n=== Family-level distribution across ALL best hits",
    "(including cultured) ===\n")
print(best_hit_with_lineage[, .N, by = family][order(-N)])

cat("\n=== Order-level distribution ===\n")
print(best_hit_with_lineage[, .N, by = order][order(-N)])

