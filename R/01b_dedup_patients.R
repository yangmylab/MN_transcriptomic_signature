#!/usr/bin/env Rscript
# Paper 1 - Phase 1 step 1b: turn the harmonised sample table into a patient-level map, quantify how
# much the 13 releases overlap, and pick one non-redundant accession per cohort/compartment/assay.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")

md <- read.csv(file.path(tab_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)

mn <- md[md$disease == "MN", ]

cat("MN samples:", nrow(mn), "  accessions:", length(unique(mn$accession)), "\n")

# ---- one row per unique MN patient per cohort family -------------------------------------------

pm <- do.call(rbind, lapply(split(mn, paste(mn$cohort, mn$patient_id, sep = "|")), function(d) {
  data.frame(cohort = d$cohort[1], patient_id = d$patient_id[1], n_samples = nrow(d),
             accessions = paste(sort(unique(d$accession)), collapse = ","),
             compartments = paste(sort(unique(d$compartment)), collapse = ","),
             assays = paste(sort(unique(d$assay)), collapse = ","),
             both_compartments = all(c("Glom", "TI") %in% d$compartment),
             stringsAsFactors = FALSE)
}))

write.csv(pm, file.path(tab_dir, "patient_map.csv"), row.names = FALSE)

cat("patient_map.csv:", nrow(pm), "unique MN patients\n")

# ---- redundancy per cohort family --------------------------------------------------------------

ds <- do.call(rbind, lapply(split(mn, mn$cohort), function(d) {
  n_pat <- length(unique(d$patient_id))
  data.frame(cohort = d$cohort[1], MN_samples = nrow(d), unique_MN_patients = n_pat,
             redundancy = round(nrow(d) / n_pat, 2),
             both_compartments = sum(sapply(split(d, d$patient_id), function(z)
               all(c("Glom", "TI") %in% z$compartment))),
             stringsAsFactors = FALSE)
}))

write.csv(ds, file.path(tab_dir, "dedup_summary.csv"), row.names = FALSE)

cat("\n=========== REDUNDANCY PER COHORT FAMILY ===========\n")

print(ds, row.names = FALSE)

# ---- subset relations between releases (MN patient sets) ----------------------------------------

set_of <- function(a) unique(mn$patient_id[mn$accession == a])

accs <- unique(mn$accession)

rel <- do.call(rbind, lapply(accs, function(a) do.call(rbind, lapply(accs, function(b) {
  if (a == b) return(NULL)
  sa <- set_of(a); sb <- set_of(b)
  data.frame(accession = a, other = b, n_a = length(sa), shared = length(intersect(sa, sb)),
             pct_of_a_in_b = round(100 * length(intersect(sa, sb)) / length(sa), 1),
             a_in_b = all(sa %in% sb), stringsAsFactors = FALSE)
}))))

write.csv(rel, file.path(tab_dir, "release_subset_relations.csv"), row.names = FALSE)

# ---- redundancy across EVERY disease class, not only MN ------------------------------------------
# The MN-only tables above cannot support a claim that overlap was assessed "across every disease
# class" (the manuscript says so). This table reports, per cohort family, how many samples map to
# how many (disease, patient) keys.

ov_all <- do.call(rbind, lapply(split(md, md$cohort), function(d) {
  d <- d[d$disease != "UNKNOWN", ]
  k <- paste(d$disease, d$patient_id)
  data.frame(cohort = d$cohort[1], samples = nrow(d), unique_patients = length(unique(k)),
             redundancy = round(nrow(d) / length(unique(k)), 2),
             diseases = length(unique(d$disease)), stringsAsFactors = FALSE)
}))

write.csv(ov_all, file.path(tab_dir, "overlap_all_diseases.csv"), row.names = FALSE)

cat("\n====== RELEASES WHOSE MN PATIENTS ARE ENTIRELY CONTAINED IN ANOTHER ======\n")

print(rel[rel$a_in_b, c("accession", "other", "n_a", "shared", "pct_of_a_in_b")], row.names = FALSE)

# ---- recommended non-redundant accession per cohort / compartment / assay -----------------------

grp <- split(mn, paste(mn$cohort, mn$compartment, mn$assay, sep = "|"))

rec <- do.call(rbind, lapply(grp, function(d) {
  n <- sapply(split(d$patient_id, d$accession), function(x) length(unique(x)))
  data.frame(cohort = d$cohort[1], compartment = d$compartment[1], assay = d$assay[1],
             keep = names(n)[which.max(n)], n_patients_kept = max(n),
             redundant = paste(names(n)[-which.max(n)], collapse = ","),
             stringsAsFactors = FALSE)
}))

write.csv(rec, file.path(tab_dir, "recommended_releases.csv"), row.names = FALSE)

cat("\n====== RECOMMENDED (non-redundant) RELEASE PER COHORT/COMPARTMENT/ASSAY ======\n")

print(rec, row.names = FALSE)

writeLines(c(
  paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 01b_dedup_patients.R"),
  paste0("  unique MN patients: ", nrow(pm)),
  paste0("  redundancy: ", paste(ds$cohort, ds$redundancy, sep = "=", collapse = "  "))
), file.path(PROJ, "logs", "dedup_patients.log"))

cat("\nDone.\n")
