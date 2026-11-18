#!/usr/bin/env Rscript
# Paper 1 - scRNA arm verifier.
# Split in two: (A) infrastructure assertions that must hold even before the data arrive, and
# (B) data assertions that run only once scrna/samples.csv lists samples.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

scrna_dir <- file.path(PROJ, "scrna")
tab_dir <- file.path(OUT, "tables")

n_pass <- 0L
n_fail <- 0L

check <- function(ok, label, detail = "") {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat(sprintf("[PASS] %-70s %s\n", label, detail))
  } else {
    n_fail <<- n_fail + 1L
    cat(sprintf("[FAIL] %-70s %s\n", label, detail))
  }
}

cat("=========== scRNA ARM VERIFICATION ===========\n\n")

# --- A. infrastructure -------------------------------------------------------

docs <- c("DATA_REQUEST.md", "DATA_REQUEST.docx", "INTEGRATION_PLAN.md", "INTEGRATION_PLAN.docx",
          "manifest.tsv")

for (d in docs) check(file.exists(file.path(scrna_dir, d)), paste("document present:", d))

scripts <- c("07a_scrna_manifest.R", "07b_scrna_qc.R", "07c_scrna_integrate_annotate.R",
             "07d_scrna_pseudobulk.R", "07e_bulk_scrna_mapping.R", "07f_deconvolution.R",
             "07z_install_scrna_pkgs.R")

for (s in scripts) {
  f <- file.path(PROJ, "R", s)
  ok <- file.exists(f) && !inherits(try(parse(f), silent = TRUE), "try-error")
  check(ok, paste("script parses:", s))
}

man <- read.delim(file.path(scrna_dir, "manifest.tsv"), stringsAsFactors = FALSE)

check(nrow(man) == 8, "manifest lists 8 datasets", nrow(man))
check(sum(man$priority == "tier1") == 3, "three tier-1 datasets", sum(man$priority == "tier1"))

ready <- read.csv(file.path(scrna_dir, "readiness.csv"), stringsAsFactors = FALSE)

check(nrow(ready) == nrow(man), "readiness covers every manifest row", nrow(ready))
check(all(c("status", "detected_samples") %in% names(ready)), "readiness reports status and counts")

# --- B. data-dependent -------------------------------------------------------

samples_file <- file.path(scrna_dir, "samples.csv")

if (!file.exists(samples_file) || !nrow(read.csv(samples_file))) {

  cat("\n[SKIP] no samples yet - tier-1 status:",
      paste(ready$status[ready$priority == "tier1"], collapse = "/"), "\n")
  cat("       place the archives from scrna/DATA_REQUEST.md and re-run 07a\n")

} else {

  samples <- read.csv(samples_file, stringsAsFactors = FALSE)

  check(nrow(samples) > 0, "samples inventoried", nrow(samples))
  check(all(file.exists(samples$path)), "every sample path exists")

  # 10x datasets are read in place: 07a records the untarred directory + a filename prefix and
  # 07b reads the triplet with ReadMtx. No per-sample link folders exist (this volume is exfat,
  # P-34), so assert the prefix actually resolves to one barcodes/features/matrix triple.
  tenx <- samples[samples$format == "10x", ]

  check("prefix" %in% names(samples), "samples.csv carries the 10x prefix column")

  trio_ok <- vapply(seq_len(nrow(tenx)), function(i) {
    f <- list.files(tenx$path[i], pattern = paste0("^", tenx$prefix[i], "_.*\\.gz$"), full.names = TRUE)
    sum(grepl("_matrix\\.mtx\\.gz$", f)) == 1 &&
      sum(grepl("_barcodes\\.tsv\\.gz$", f)) == 1 &&
      sum(grepl("_(features|genes)\\.tsv\\.gz$", f)) == 1
  }, logical(1))

  check(all(trio_ok), "every 10x prefix resolves to one barcode/feature/matrix triple",
        sprintf("%d/%d", sum(trio_ok), length(trio_ok)))

  check(!any(dir.exists(file.path(unique(tenx$path), "samples"))),
        "no per-sample link folders under extracted (exfat cannot hold links)")

  check(all(file.exists(file.path(PROJ, man$archive[man$priority == "tier1"]))),
        "tier-1 archives resolve from scrna/manifest.tsv")
  check(sum(!is.na(samples$disease_group)) > 0, "disease labels attached",
        sum(!is.na(samples$disease_group)))

  qc <- read.csv(file.path(tab_dir, "scrna_qc.csv"), stringsAsFactors = FALSE)
  check(nrow(qc) == nrow(samples), "QC table covers every sample", nrow(qc))
  check(sum(qc$cells_kept) > 1000, "cells retained", sum(qc$cells_kept))

  check(file.exists(file.path(PROJ, "work", "scrna", "atlas.rds")), "atlas object written")
  check(file.exists(file.path(PROJ, "work", "scrna", "atlas_cells.csv")), "cell metadata written")

  counts <- read.csv(file.path(tab_dir, "scrna_celltype_counts.csv"), stringsAsFactors = FALSE)
  check(nrow(counts) > 0, "cell-type counts written", nrow(counts))
  check(all(c("Podocyte", "Endothelium", "Myeloid") %in% unique(counts$celltype)),
        "expected kidney lineages recovered",
        paste(intersect(c("Podocyte", "Endothelium", "Myeloid"),
                        unique(counts$celltype)), collapse = ","))

  check(file.exists(file.path(tab_dir, "scrna_pseudobulk_de.csv")), "pseudobulk DE written")
  check(file.exists(file.path(tab_dir, "signature_matrix.csv")), "signature matrix written")

  sig <- read.csv(file.path(tab_dir, "signature_matrix.csv"), row.names = 1, check.names = FALSE)

  check(nrow(sig) > 5000, "signature matrix has enough genes", nrow(sig))
  check(ncol(sig) >= 8, "signature matrix has enough cell types", ncol(sig))

  check(file.exists(file.path(tab_dir, "bulk_deg_celltype_attribution.csv")),
        "bulk DEG attribution written")
  check(file.exists(file.path(tab_dir, "bulk_vs_pseudobulk_concordance.csv")),
        "bulk-vs-pseudobulk concordance written")
  check(file.exists(file.path(tab_dir, "deconvolution_fractions.csv")),
        "deconvolution fractions written")

  fr <- read.csv(file.path(tab_dir, "deconvolution_fractions.csv"), stringsAsFactors = FALSE)
  ct_cols <- setdiff(names(fr), c("sample", "set", "disease", "group", "rms_residual", "max_weight"))
  check(all(abs(rowSums(fr[, ct_cols]) - 1) < 1e-6), "deconvolved proportions sum to 1")

  # Regression guards for the two deconvolution faults found on 2026-10-19 (handoff P-28/P-29):
  # an unscaled platform mismatch collapsed every fit onto a single cell type, and the unfiltered
  # all-gene signature made every sample 66-88 % "Neutrophil".
  check("max_weight" %in% names(fr), "deconvolution records a fit-quality column")
  check(mean(fr$max_weight > 0.99) < 0.02, "no degenerate single-cell-type fits",
        sprintf("%.1f%% corner solutions", 100 * mean(fr$max_weight > 0.99)))

  check(file.exists(file.path(tab_dir, "deconvolution_fractions_allgenes.csv")),
        "all-gene sensitivity deconvolution written")
  check(file.exists(file.path(tab_dir, "deconvolution_signature_sensitivity.csv")),
        "signature sensitivity comparison written")

  glom <- fr[fr$set == "disc_glom_array", ]

  check(median(glom$Podocyte) > median(glom$Neutrophil),
        "glomerular deconvolution is podocyte-dominant, not neutrophil-dominated",
        sprintf("podocyte %.3f vs neutrophil %.3f", median(glom$Podocyte), median(glom$Neutrophil)))
}

cat("\n=============================================\n")
cat(sprintf("passed: %d | failed: %d\n", n_pass, n_fail))

if (n_fail == 0) {
  cat("ALL CHECKS PASSED\n")
  quit(status = 0)
} else {
  cat("CHECKS FAILED\n")
  quit(status = 1)
}
