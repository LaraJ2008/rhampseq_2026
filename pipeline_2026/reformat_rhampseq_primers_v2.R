## =============================================================================
## Reformat rhAmpSeq primer export (v2) into a tidy per-target primer table
## =============================================================================
## Matches actual CSV structure:
##   reverse_name, reverse_sequence, forward_name , forward_sequence
## Note: some header/value fields carry stray leading/trailing whitespace
## from the Excel export -- trimmed explicitly below.
##
## Raw sequence format example:
##   /rhSeq-r/ACGCTGGTAGTGCGTCrGGGTG/GT3/
## Three parts:
##   /rhSeq-r/    -- chemistry-class tag (reverse vs forward), not sequence
##   ...GTCrGGGTG -- real sequence; lowercase 'r' is a prefix marker on the
##                   base immediately following it ("rG" = ribo-G, the single
##                   RNA base central to rhAmp's cleavage chemistry). Dropping
##                   just the lowercase 'r' and keeping the base after it
##                   gives a clean DNA sequence for alignment.
##   /GT3/        -- trailing IDT blocking-moiety/modification code. Preserved
##                   in its own column this time (not discarded) in case
##                   moiety code correlates with anything downstream (e.g.
##                   which primers survived multiplex compatibility QC).
## =============================================================================

library(data.table)
library(stringr)

## ---- 0. CONFIG -------------------------------------------------------------

out_csv <- "rhampseq_primers_clean.csv"

## ---- 1. Read and clean up whitespace ----------------------------------------
raw_csv <- "rhamp_seq_primer_list_full.csv"   # the filename, as a string
raw <- read.csv(raw_csv, header = TRUE, stringsAsFactors = FALSE)
names(raw) <- str_trim(names(raw))
raw <- as.data.table(raw)


# trim whitespace from every character column's VALUES too
char_cols <- names(raw)[sapply(raw, is.character)]
raw[, (char_cols) := lapply(.SD, str_trim), .SDcols = char_cols]

stopifnot(all(c("reverse_name","reverse_sequence",
                "forward_name","forward_sequence") %in% names(raw)))

## ---- 2. Helper: strip tags, extract clean DNA seq + RNA-base index + moiety -

parse_primer_seq <- function(raw_seq) {
  # blocking moiety: whatever's inside the trailing /GT<digits>/ tag
  moiety <- str_extract(raw_seq, "(?<=/)GT\\d+(?=/$)")

  # strip leading /rhSeq-r/ or /rhSeq-f/ tag
  s <- str_remove(raw_seq, "^/rhSeq-[rf]/")
  # strip trailing /GT<digits>/ tag
  s <- str_remove(s, "/GT\\d+/$")

  # position of the lowercase 'r' (RNA base marker), 1-based from 5' end,
  # measured AFTER tag-stripping but BEFORE removing the 'r' itself
  rna_pos_5p <- str_locate(s, "r")[, "start"]

  # drop just the lowercase 'r' marker, keep the base that follows
  s_clean <- str_remove(s, "r")

  list(clean_seq   = s_clean,
       blocking_moiety = moiety,
       rna_base_pos_5prime = rna_pos_5p,
       rna_base_pos_from_3prime = nchar(s_clean) - rna_pos_5p + 1)
}

## ---- 3. Extract target_id (strip trailing .R / .F from names) --------------

raw[, target_id := str_remove(reverse_name, "\\.R$")]
stopifnot(all(raw$target_id == str_remove(raw$forward_name, "\\.F$")))

## ---- 4. Apply parsing to both primer columns --------------------------------

r_parsed <- rbindlist(lapply(raw$reverse_sequence, parse_primer_seq))
f_parsed <- rbindlist(lapply(raw$forward_sequence, parse_primer_seq))

primers <- data.table(
  target_id              = raw$target_id,

  primer_R_raw            = raw$reverse_sequence,
  primer_R_clean          = r_parsed$clean_seq,
  primer_R_blocking_moiety= r_parsed$blocking_moiety,
  primer_R_rna_pos_5p     = r_parsed$rna_base_pos_5prime,
  primer_R_rna_pos_3p     = r_parsed$rna_base_pos_from_3prime,

  primer_F_raw            = raw$forward_sequence,
  primer_F_clean          = f_parsed$clean_seq,
  primer_F_blocking_moiety= f_parsed$blocking_moiety,
  primer_F_rna_pos_5p     = f_parsed$rna_base_pos_5prime,
  primer_F_rna_pos_3p     = f_parsed$rna_base_pos_from_3prime
)

## ---- 5. Sanity checks --------------------------------------------------------

cat("Rows parsed:", nrow(primers), "\n\n")

n_missing_rna <- primers[is.na(primer_R_rna_pos_5p) | is.na(primer_F_rna_pos_5p), .N]
cat("Rows with a missing/failed RNA-base detection (should be 0):", n_missing_rna, "\n")
if (n_missing_rna > 0) print(primers[is.na(primer_R_rna_pos_5p) | is.na(primer_F_rna_pos_5p)])

n_missing_moiety <- primers[is.na(primer_R_blocking_moiety) | is.na(primer_F_blocking_moiety), .N]
cat("Rows with a missing blocking-moiety tag (should be 0):", n_missing_moiety, "\n\n")

cat("Blocking moiety code distribution (R primers):\n")
print(table(primers$primer_R_blocking_moiety, useNA = "ifany"))
cat("\nBlocking moiety code distribution (F primers):\n")
print(table(primers$primer_F_blocking_moiety, useNA = "ifany"))

cat("\nSample of cleaned output:\n")
print(head(primers[, .(target_id, primer_R_clean, primer_R_blocking_moiety,
                        primer_R_rna_pos_3p, primer_F_clean,
                        primer_F_blocking_moiety, primer_F_rna_pos_3p)]))

fwrite(primers, out_csv)
cat("\nSaved tidy primer table to", out_csv, "\n")
