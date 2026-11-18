#!/usr/bin/env Rscript
# Paper 1 - Phase 1 step 1a: harmonise every sample across the 13 accessions into one metadata
# table (accession, GSM, patient, compartment, disease, platform, assay), plus an accession
# summary and a pairwise sample-overlap matrix that exposes duplicated patients between releases.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")

dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---- disease vocabulary ------------------------------------------------------------------------
# Each cohort spells the same diagnosis differently, so collapse them onto one canonical label.

disease_map <- c(
  "mn" = "MN", "mgn" = "MN", "membranous nephropathy" = "MN",
  "membranous glomerulonephropathy" = "MN", "idiopathic membranous nephropathy" = "MN",
  "fsgs" = "FSGS", "focal segmental glomerulosclerosis" = "FSGS",
  "focal and segmental glomerulosclerosis" = "FSGS", "focal segmental glomerular sclerosis" = "FSGS",
  "mcd" = "MCD", "minimal change disease" = "MCD",
  "fsgs_mcd" = "FSGS_MCD",
  "other" = "Other", "other nephrotic syndrome" = "Other",
  "iga" = "IgAN", "iga nephropathy" = "IgAN",
  "sle" = "SLE", "systemic lupus erythematosus" = "SLE",
  "rpgn" = "RPGN", "anca" = "ANCA", "aav" = "ANCA", "anca-associated vasculitis" = "ANCA",
  "dn" = "DN", "diabetic nephropathy" = "DN",
  "ht" = "HT", "hypertensive nephropathy" = "HT",
  "tmd" = "TMD", "thin membrande disease" = "TMD", "thin membrane disease" = "TMD",
  "cd" = "CD",
  "tn" = "Control-TN", "tumor nephrectomy" = "Control-TN",
  "ld" = "Control-LD", "living donor" = "Control-LD", "livingdonor" = "Control-LD",
  "living_donor" = "Control-LD", "healthy living transplant donor" = "Control-LD",
  "healthy control" = "Control-LD", "nor" = "Control-LD"
)

canon_dis <- function(x) {
  k <- tolower(gsub("[[:space:]_-]+$", "", trimws(x)))
  out <- disease_map[k]
  out[is.na(out)] <- "UNKNOWN"
  unname(out)
}

# ---- generic GEO series-matrix sample-metadata reader ------------------------------------------

read_sm <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  lines <- lines[startsWith(lines, "!Sample_")]
  keys <- sub("^!", "", sub("\t.*$", "", lines))
  vals <- strsplit(sub("^[^\t]*\t", "", lines), "\t", fixed = TRUE)
  vals <- lapply(vals, gsub, pattern = '^"|"$', replacement = "")
  out <- list()
  for (k in unique(keys)) {
    idx <- which(keys == k)
    m <- do.call(rbind, vals[idx])
    out[[k]] <- if (ncol(m) == 1) m[1, ] else apply(m, 2, function(z) paste(z[nzchar(z)], collapse = " | "))
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

sm_disease <- function(sm) {
  chr_cols <- grep("characteristics|source_name", names(sm), value = TRUE)
  blob <- do.call(paste, c(sm[, chr_cols, drop = FALSE], list(sep = " | ")))
  m <- regmatches(blob, regexpr("(disease state|diagnosis|disease|phenotype):\\s*[^|]+", blob))
  sub("^[^:]*:\\s*", "", m)
}

# ---- array cohorts: the CEL file name encodes compartment, disease and patient -----------------
# NEPTUNE-numbered releases:  B1_TI_MN_94   /  210_Other  /  B1_Glom_FSGS_1
# ERCB-numbered releases:     H5-Glom-LD_478  /  H7-Tub-FSGS_MCD905
# 'AAV' (ANCA-associated vasculitis) and 'FSGS_MCD' (a mixed diagnostic label) also occur.

disease_re <- "(FSGS_MCD|Living_donor|LivingDonor|MGN|MN|FSGS|MCD|Other|Nor|LD|AAV|SLE|IgA|RPGN|ANCA|DN|HT|TMD|TN|CD)"

extract_dis <- function(t) {
  m <- regmatches(t, regexpr(disease_re, t))
  if (length(m)) m else NA_character_
}

extract_num <- function(t) sub("^.*?([0-9]+)[^0-9]*$", "\\1", t)

extract_bat <- function(t) {
  m <- regmatches(t, regexpr("^[BH][0-9]+", t))
  if (length(m)) m else NA_character_
}

parse_array <- function(acc, comp, platform) {
  fs <- list.files(file.path(DATA, acc, "raw"), pattern = "\\.CEL\\.gz$", full.names = TRUE)
  t <- vapply(fs, function(f) sub("\\.CEL\\.gz$", "", sub("^GSM[0-9]+_", "", basename(f))), "")
  dis <- vapply(t, extract_dis, "")
  num <- vapply(t, extract_num, "")
  bat <- vapply(t, extract_bat, "")
  ercb <- grepl("^[BH][0-9]+-", t)
  cmp <- ifelse(ercb, ifelse(grepl("-Tub-", t), "TI", "Glom"), comp)
  # living donors are numbered 1..n independently of the patients, so raw donor IDs collide with
  # real patient IDs (1_FSGS vs 1_Living_donor). Namespace them to keep the two apart.
  pid <- ifelse(ercb, paste(bat, dis, num, sep = "-"), num)
  ctrl <- !ercb & canon_dis(dis) == "Control-LD"
  pid[ctrl] <- paste0("LD-", num[ctrl])
  data.frame(accession = acc, gsm = sub("_.*$", "", basename(fs)), sample_title = t,
             patient_id = pid,
             compartment = cmp, disease = canon_dis(dis), disease_raw = dis,
             batch = bat, assay = "microarray", platform = platform, stringsAsFactors = FALSE)
}

# ---- RNA-seq cohorts ---------------------------------------------------------------------------

parse_rnaseq <- function(acc, comp, counts_gz) {
  cols <- strsplit(readLines(gzfile(counts_gz), n = 1), "\t")[[1]][-1]
  if (acc == "GSE216841") {
    dis <- ifelse(grepl("^MEMBR", cols), "MN", ifelse(grepl("^MCD", cols), "MCD", "Control-LD"))
    return(data.frame(accession = acc, gsm = cols, sample_title = cols, patient_id = cols,
                      compartment = comp, disease = dis, disease_raw = dis,
                      batch = NA_character_, assay = "RNA-seq", platform = "GPL18573",
                      stringsAsFactors = FALSE))
  }
  sm <- read_sm(file.path(DATA, acc, paste0(acc, "_series_matrix.txt.gz")))
  dmap <- setNames(sm_disease(sm), sub(" .*$", "", sm[["Sample_title"]]))
  num <- sub("^Sample_([0-9]+).*$", "\\1", cols)
  # the eight living donors are numbered 1..8 and would otherwise collide with patients 1..8
  ctrl <- grepl("LivingDonor", cols)
  num[ctrl] <- paste0("LD-", num[ctrl])
  data.frame(accession = acc, gsm = cols, sample_title = cols, patient_id = num,
             compartment = comp, disease = canon_dis(unname(dmap[cols])),
             disease_raw = unname(dmap[cols]), batch = NA_character_, assay = "RNA-seq",
             platform = "GPL18573", stringsAsFactors = FALSE)
}

# ---- Agilent PBMC cohort (GSE73953): only the series matrix carries the disease -----------------

parse_agilent <- function(acc) {
  sm <- read_sm(file.path(DATA, acc, paste0(acc, "_series_matrix.txt.gz")))
  data.frame(accession = acc, gsm = sub("_.*$", "", sm[["Sample_title"]]),
             sample_title = sm[["Sample_title"]],
             patient_id = paste0("PBMC-", seq_len(nrow(sm))),
             compartment = "PBMC", disease = canon_dis(sm_disease(sm)),
             disease_raw = sm_disease(sm), batch = NA_character_, assay = "microarray",
             platform = "GPL4133", stringsAsFactors = FALSE)
}

# ---- assemble ----------------------------------------------------------------------------------

md <- rbind(
  parse_array("GSE133288", "TI", "GPL19983"),
  parse_array("GSE200818", "TI", "GPL19983"),
  parse_array("GSE108112", "TI", "GPL19983"),
  parse_array("GSE200828", "Glom", "GPL19983"),
  parse_array("GSE108109", "Glom", "GPL19983"),
  parse_array("GSE104948", "Glom", "GPL22945/24120"),
  parse_array("GSE104954", "TI", "GPL22945/24120"),
  parse_array("GSE99339", "Glom", "GPL19109/19184"),
  parse_array("GSE99325", "TI", "GPL19109/19184"),
  parse_agilent("GSE73953"),
  parse_rnaseq("GSE197307", "Glom", file.path(DATA, "GSE197307", "processed", "GSE197307_counts_ensg.txt.gz")),
  parse_rnaseq("GSE182380", "TI", file.path(DATA, "GSE182380", "processed", "GSE182380_counts_symbol.txt.gz")),
  parse_rnaseq("GSE216841", "Glom", file.path(DATA, "GSE216841", "processed", "GSE216841_counts_ensg.txt.gz"))
)

md$cohort <- ifelse(grepl("^GSE(104948|104954|99339|99325)$", md$accession), "ERCB-numbered",
                    ifelse(md$accession == "GSE216841", "LCM-FFPE",
                           ifelse(md$accession == "GSE73953", "PBMC", "NEPTUNE-numbered")))

write.csv(md, file.path(tab_dir, "sample_metadata.csv"), row.names = FALSE)

cat("sample_metadata.csv:", nrow(md), "samples across", length(unique(md$accession)), "accessions\n")

# ---- per-accession summary ---------------------------------------------------------------------

acc_sum <- do.call(rbind, lapply(split(md, md$accession), function(d) {
  data.frame(accession = d$accession[1], cohort = d$cohort[1], assay = d$assay[1],
             compartment = paste(sort(unique(d$compartment)), collapse = "+"),
             platform = d$platform[1], n_samples = nrow(d),
             n_MN_samples = sum(d$disease == "MN"),
             n_MN_patients = length(unique(d$patient_id[d$disease == "MN"])),
             n_UNKNOWN = sum(d$disease == "UNKNOWN"),
             diseases = paste(names(sort(table(d$disease), decreasing = TRUE)), collapse = ","),
             stringsAsFactors = FALSE)
}))

write.csv(acc_sum, file.path(tab_dir, "accession_summary.csv"), row.names = FALSE)

cat("\n================ ACCESSION SUMMARY ================\n")

print(acc_sum[, c("accession", "cohort", "assay", "compartment", "n_samples", "n_MN_samples",
                  "n_MN_patients", "n_UNKNOWN", "diseases")], row.names = FALSE)

# ---- pairwise patient overlap among MN samples --------------------------------------------------
# Within each cohort family the patient identifier is comparable, so shared IDs reveal releases that
# profile the same people (possibly a different compartment or platform).

mn <- md[md$disease == "MN", ]

fam <- unique(mn$accession)

ov <- sapply(fam, function(a) sapply(fam, function(b)
  length(intersect(mn$patient_id[mn$accession == a], mn$patient_id[mn$accession == b]))))

write.csv(ov, file.path(tab_dir, "mn_patient_overlap_matrix.csv"))

cat("\n====== MN PATIENT-ID OVERLAP BETWEEN RELEASES ======\n")

print(ov)

# ---- the number that matters: unique MN patients per cohort family ------------------------------

cat("\n================ UNIQUE MN PATIENTS (de-duplicated) ================\n")

for (cf in c("NEPTUNE-numbered", "ERCB-numbered", "LCM-FFPE")) {
  d <- mn[mn$cohort == cf, ]
  cat(sprintf("%-18s MN samples=%-5d unique MN patients=%-5d\n", cf, nrow(d), length(unique(d$patient_id))))
}

cat("\nPBMC (GSE73953):", sum(md$accession == "GSE73953" & md$disease == "MN"), "MN samples (peripheral, optional)\n")

cat("\nUNKNOWN-disease samples:", sum(md$disease == "UNKNOWN"), "\n")

writeLines(c(
  paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 01a_build_sample_metadata.R"),
  paste0("  samples: ", nrow(md), "  accessions: ", length(unique(md$accession))),
  paste0("  UNKNOWN disease labels: ", sum(md$disease == "UNKNOWN")),
  paste0("  tables: ", tab_dir)
), file.path(PROJ, "logs", "build_sample_metadata.log"))

cat("\nDone.\n")
