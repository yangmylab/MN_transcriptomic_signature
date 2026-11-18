#!/usr/bin/env Rscript
# Paper 1 - Phase 1 step 1c: composition of the proposed non-redundant analysis tiers.
# For releases that overlap, keep one sample per patient (priority order below) so that no patient
# is counted twice inside a tier.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")

md <- read.csv(file.path(tab_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)

md$patient_key <- paste(md$cohort, md$patient_id)

# priority: within a tier, earlier accessions win for a duplicated patient
priority <- c("GSE133288", "GSE200818", "GSE200828", "GSE104948", "GSE104954", "GSE99339",
              "GSE99325", "GSE108112", "GSE108109", "GSE197307", "GSE182380", "GSE216841", "GSE73953")

md$prio <- match(md$accession, priority)

tiers <- list(
  "Discovery glom (microarray, NEPTUNE)"    = c("GSE200828", "GSE108109"),
  "Discovery TI (microarray, NEPTUNE)"      = c("GSE133288", "GSE200818", "GSE108112"),
  "Validation glom (microarray, ERCB)"      = c("GSE104948", "GSE99339"),
  "Validation TI (microarray, ERCB)"        = c("GSE104954", "GSE99325"),
  "Cross-platform glom (RNA-seq, NEPTUNE)"  = "GSE197307",
  "Cross-platform TI (RNA-seq, NEPTUNE)"    = "GSE182380",
  "PLA2R sub-analysis (LCM-FFPE)"           = "GSE216841",
  "Peripheral cross-check (PBMC)"           = "GSE73953"
)

rows <- list()

for (nm in names(tiers)) {
  d <- md[md$accession %in% tiers[[nm]], ]
  d <- d[order(d$patient_key, d$prio), ]
  d <- d[!duplicated(d$patient_key), ]
  comp <- as.data.frame.matrix(table(d$disease, d$compartment))
  comp$disease <- rownames(comp)
  rows[[nm]] <- data.frame(tier = nm, releases = paste(tiers[[nm]], collapse = "+"),
                           n_samples = nrow(d), n_patients = length(unique(d$patient_key)),
                           MN = sum(d$disease == "MN"),
                           control_LD = sum(d$disease == "Control-LD"),
                           FSGS = sum(d$disease == "FSGS"), MCD = sum(d$disease == "MCD"),
                           Other = sum(d$disease == "Other"),
                           stringsAsFactors = FALSE)
  cat("\n=====", nm, "=====\n")
  print(comp, row.names = FALSE)
}

comp_tab <- do.call(rbind, rows)

write.csv(comp_tab, file.path(tab_dir, "proposed_tiers.csv"), row.names = FALSE)

cat("\n================ PROPOSED NON-REDUNDANT TIERS ================\n")

print(comp_tab, row.names = FALSE)

# how many patients overlap between the NEPTUNE microarray discovery set and each validation layer?

disc <- md[md$accession %in% c("GSE200828", "GSE108109", "GSE133288", "GSE200818", "GSE108112"), ]

cat("\n====== PATIENT OVERLAP: NEPTUNE MICROARRAY DISCOVERY vs EACH VALIDATION LAYER ======\n")

for (acc in c("GSE104948", "GSE104954", "GSE197307", "GSE182380", "GSE216841", "GSE73953")) {
  v <- md[md$accession == acc & md$disease == "MN", ]
  d <- disc[disc$disease == "MN", ]
  cat(sprintf("%-10s  MN patients=%-3d  shared with discovery=%-3d  (%.0f%%)\n",
              acc, length(unique(v$patient_key)),
              length(intersect(v$patient_key, d$patient_key)),
              100 * length(intersect(v$patient_key, d$patient_key)) / length(unique(v$patient_key))))
}

writeLines(c(
  paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 01c_design_options.R"),
  paste0("  tiers: ", length(tiers)),
  paste0("  total non-redundant MN (disc+val+cross): ",
         sum(comp_tab$MN[1:6]))
), file.path(PROJ, "logs", "design_options.log"))

cat("\nDone.\n")
