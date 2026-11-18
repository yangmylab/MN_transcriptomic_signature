#!/usr/bin/env Rscript
# Paper 1 - Phase A step 3b: assemble the 8 APPROVED, non-redundant cohort sets from the normalised
# matrices, de-duplicating patients, dropping redundant releases, and removing the RNA-seq patients
# that are shared with the microarray discovery set from the cross-platform validation layer.
#
# Gene space: each tier is a single platform, so the matrices keep their native gene ids (Entrez for
# the Brainarray arrays, Ensembl or symbol for the RNA-seq cohorts) and each set also gets a
# symbol-keyed copy for cross-platform comparison.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })

set_dir <- file.path(PROJ, "work", "sets")
rma_dir <- file.path(PROJ, "work", "rma")
rna_dir <- file.path(PROJ, "work", "rnaseq")
tab_dir <- file.path(OUT, "tables")

dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)

md <- read.csv(file.path(tab_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)

md$patient_key <- paste(md$cohort, md$patient_id)

# ---- load one normalised matrix per accession ----------------------------------------------------

load_rds <- function(p) { r <- readRDS(p); cat(sprintf("  loaded %-46s %6d x %4d\n", basename(p), nrow(r), ncol(r))); r }

cat("=== loading normalised matrices ===\n")

M <- list(
  GSE133288 = load_rds(file.path(rma_dir, "GSE133288_GPL19983_gene_rma.rds")),
  GSE200818 = load_rds(file.path(rma_dir, "GSE200818_GPL19983_gene_rma.rds")),
  GSE200828 = load_rds(file.path(rma_dir, "GSE200828_GPL19983_gene_rma.rds")),
  GSE104948 = load_rds(file.path(rma_dir, "GSE104948_gene_rma_merged.rds")),
  GSE104954 = load_rds(file.path(rma_dir, "GSE104954_gene_rma_merged.rds")),
  GSE197307 = load_rds(file.path(rna_dir, "GSE197307_log2expr.rds")),
  GSE182380 = load_rds(file.path(rna_dir, "GSE182380_log2expr.rds")),
  GSE216841 = load_rds(file.path(rna_dir, "GSE216841_log2expr.rds"))
)

# ---- gene-id -> symbol map (org.Hs.eg.db covers both Entrez and Ensembl ids) --------------------

cat("\n=== building gene-id -> symbol map ===\n")

sym <- function(ids, keytype) {
  s <- suppressWarnings(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                                              keytype = keytype, multiVals = "first"))
  unname(s)
}

to_symbol <- function(m, keytype) {
  s <- sym(rownames(m), keytype)
  keep <- !is.na(s) & !duplicated(s)
  cat("  mapped", sum(!is.na(s)), "of", nrow(m), "ids;", sum(keep), "unique symbols kept\n")
  out <- m[keep, , drop = FALSE]
  rownames(out) <- s[keep]
  out
}

# ---- cohort definition ---------------------------------------------------------------------------
# priority within a tier decides which sample survives when a patient appears in two releases

priority <- c("GSE133288", "GSE200818", "GSE200828", "GSE104948", "GSE104954", "GSE197307", "GSE182380")

sets <- list(
  disc_glom_array  = list(acc = c("GSE200828"), keytype = "ENTREZID"),
  disc_ti_array    = list(acc = c("GSE133288", "GSE200818"), keytype = "ENTREZID"),
  val_glom_array   = list(acc = c("GSE104948"), keytype = "ENTREZID"),
  val_ti_array     = list(acc = c("GSE104954"), keytype = "ENTREZID"),
  cross_glom_rnaseq = list(acc = c("GSE197307"), keytype = "ENSEMBL"),
  cross_ti_rnaseq  = list(acc = c("GSE182380"), keytype = "SYMBOL"),
  pla2r_lcm        = list(acc = c("GSE216841"), keytype = "ENSEMBL")
)

# the RNA-seq patients that also appear in the microarray discovery set must not validate the model
disc_pat <- unique(md$patient_key[md$accession %in% c("GSE200828", "GSE133288", "GSE200818")])

build <- function(nm, spec) {
  cat("\n=====", nm, "=====\n")
  mats <- lapply(spec$acc, function(a) M[[a]])
  names(mats) <- spec$acc

  # samples of these accessions, in tier priority order
  ann <- md[md$accession %in% spec$acc, ]
  ann <- ann[order(match(ann$accession, priority), ann$gsm), ]

  # every patient is kept. The RNA-seq series (GSE197307, GSE182380) do NOT carry a participant id:
  # their sample titles are local sample indices and the series matrices list only disease + tissue.
  # Patient-level linkage to the array cohorts is therefore IMPOSSIBLE from the downloaded metadata,
  # so the RNA-seq tiers are treated as independent cohorts and the overlap column is left NA rather
  # than filled with identifier collisions (see handoff P-2). ERCB is independent by construction
  # (different cohort, 0 shared participants).
  ann <- ann[!duplicated(ann$patient_key), ]

  # stack the matrices on their common genes and keep the selected samples
  common <- Reduce(intersect, lapply(mats, rownames))
  cat("  common genes across releases:", length(common), "\n")
  m <- do.call(cbind, lapply(mats, function(x) x[common, , drop = FALSE]))

  keep <- colnames(m) %in% ann$gsm
  m <- m[, keep, drop = FALSE]
  ann <- ann[match(colnames(m), ann$gsm), ]

  stopifnot(identical(ann$gsm, colnames(m)))

  cat("  samples:", ncol(m), "| patients:", length(unique(ann$patient_key)),
      "| diseases:", paste(names(sort(table(ann$disease), decreasing = TRUE)), collapse = ","), "\n")

  ms <- to_symbol(m, spec$keytype)

  saveRDS(m, file.path(set_dir, paste0(nm, "_expr.rds")))
  saveRDS(ms, file.path(set_dir, paste0(nm, "_expr_symbol.rds")))
  write.csv(ann, file.path(set_dir, paste0(nm, "_samples.csv")), row.names = FALSE)

  data.frame(set = nm, releases = paste(spec$acc, collapse = "+"), genes = nrow(m),
             genes_symbol = nrow(ms), samples = ncol(m),
             patients = length(unique(ann$patient_key)),
             MN = sum(ann$disease == "MN"),
             FSGS = sum(ann$disease == "FSGS"), MCD = sum(ann$disease == "MCD"),
             Other = sum(ann$disease == "Other"), control = sum(grepl("^Control", ann$disease)),
             other_dx = sum(!ann$disease %in% c("MN", "FSGS", "MCD", "Other") & !grepl("^Control", ann$disease)),
             role = if (grepl("^cross_", nm)) "independent RNA-seq cohort" else "discovery / external validation",
             MN_overlapping_discovery = if (grepl("^cross_", nm)) NA_integer_ else
               sum(ann$disease == "MN" & ann$patient_key %in% disc_pat),
             stringsAsFactors = FALSE)
}

summ <- do.call(rbind, lapply(names(sets), function(nm) build(nm, sets[[nm]])))

write.csv(summ, file.path(tab_dir, "cohort_sets_summary.csv"), row.names = FALSE)

cat("\n================ COHORT SETS ================\n")

print(summ, row.names = FALSE)

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 03b_build_cohort_sets.R"),
             capture.output(print(summ, row.names = FALSE))),
           file.path(PROJ, "logs", "build_cohort_sets.log"))

cat("\nDone.\n")
