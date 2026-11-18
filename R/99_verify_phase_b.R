#!/usr/bin/env Rscript
# Paper 1 - Phase B verifier. Every assertion is a claim the Phase B analysis makes; the script
# exits 0 only when all of them hold.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")
deg_dir <- file.path(tab_dir, "deg")
meta_dir <- file.path(tab_dir, "meta")
enr_dir <- file.path(tab_dir, "enrichment")
fig_enr <- file.path(OUT, "enrichment")
fig_deg <- file.path(OUT, "deg")

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

cat("=========== PHASE B VERIFICATION ===========\n\n")

# --- QC policy ---------------------------------------------------------------
pol <- read.csv(file.path(tab_dir, "qc_outlier_policy.csv"), stringsAsFactors = FALSE)

check(nrow(pol) == 7, "QC policy covers all 7 cohort sets", nrow(pol))
check(all(pol$decision == "retain"), "QC policy retains every sample")
check(sum(pol$flagged_outliers) == 200, "QC policy documents 200 flagged samples", sum(pol$flagged_outliers))

flags <- read.csv(file.path(tab_dir, "qc_flagged_samples_retained.csv"), stringsAsFactors = FALSE)

check(nrow(flags) == 200, "flagged-sample appendix lists 200 rows", nrow(flags))

# --- step 6: ComBat ----------------------------------------------------------
combat <- read.csv(file.path(tab_dir, "combat_summary.csv"), stringsAsFactors = FALSE)

check(identical(combat$set, c("disc_ti_array", "disc_glom_array")), "combat_summary covers both discovery sets")

check(file.exists(file.path(PROJ, "work", "combat", "disc_ti_array_expr_combat.rds")),
      "ComBat-adjusted Entrez matrix exists")
check(file.exists(file.path(PROJ, "work", "combat", "disc_ti_array_expr_combat_symbol.rds")),
      "ComBat-adjusted symbol matrix exists")

adj <- readRDS(file.path(PROJ, "work", "combat", "disc_ti_array_expr_combat_symbol.rds"))
raw <- readRDS(file.path(PROJ, "work", "sets", "disc_ti_array_expr_symbol.rds"))

check(identical(dim(adj), dim(raw)), "ComBat matrix keeps dimensions", paste(dim(adj), collapse = "x"))
check(identical(colnames(adj), colnames(raw)), "ComBat matrix keeps sample order")
check(!identical(adj, raw), "ComBat actually changed the matrix")

b <- combat$batch_r2_pc1to5_before[1]
a <- combat$batch_r2_pc1to5_after[1]

check(a < b, "batch R2 falls after ComBat", paste(round(b, 4), "->", round(a, 4)))
check(a < 0.01, "batch R2 after ComBat is negligible", round(a, 4))
check(abs(combat$disease_silhouette_before[1] - combat$disease_silhouette_after[1]) < 0.05,
      "disease signal is preserved (silhouette unchanged)")
check(combat$action[2] == "no ComBat: single series, no technical batch to remove",
      "single-series glom set is documented, not silently skipped")

# --- step 7: limma -----------------------------------------------------------
deg_sum <- read.csv(file.path(tab_dir, "deg_summary.csv"), stringsAsFactors = FALSE)

check(nrow(deg_sum) == 32, "deg_summary has 31 cohort x contrast rows (+header)", nrow(deg_sum))
check(length(unique(deg_sum$set)) == 7, "all 7 cohorts contributed",
      length(unique(deg_sum$set)))
check(all(c("Glom", "TI") %in% deg_sum$compartment), "both compartments represented")
check(all(c("MN_vs_Control", "MN_vs_FSGS", "MN_vs_MCD", "MN_vs_Other", "MN_vs_nonMN") %in%
            deg_sum$contrast), "all five contrasts present")
check(sum(deg_sum$deg_up) + sum(deg_sum$deg_down) == sum(deg_sum$deg_total),
      "up + down equals total in every row")

files <- list.files(deg_dir, pattern = "\\.csv$")

check(length(files) == 32, "32 per-contrast DEG tables", length(files))

one <- read.csv(file.path(deg_dir, "disc_glom_array__MN_vs_Control.csv"), stringsAsFactors = FALSE)

check(all(c("gene", "logFC", "SE", "t", "P.Value", "adj.P.Val", "deg", "direction") %in%
            names(one)), "DEG tables carry logFC, SE, t, P, BH P and a flag")
check(!anyNA(one$SE), "limma SE is present for every gene")
check(all(abs(one$logFC[one$deg]) > 1 & one$adj.P.Val[one$deg] < 0.05),
      "every flagged DEG meets |log2FC| > 1 and BH < 0.05")

sens <- read.csv(file.path(tab_dir, "deg_combat_sensitivity.csv"), stringsAsFactors = FALSE)

check(all(sens$lfc_pearson > 0.95), "ComBat leaves logFC essentially unchanged",
      paste(range(sens$lfc_pearson), collapse = "-"))

# --- step 8: meta-analysis ---------------------------------------------------
meta <- read.csv(file.path(tab_dir, "meta_summary.csv"), stringsAsFactors = FALSE)

check(nrow(meta) == 10, "10 primary metas (5 contrasts x 2 compartments)", nrow(meta))
check(all(meta$meta_deg_up + meta$meta_deg_down == meta$meta_deg_total),
      "meta-DEG up + down equals total")

for (i in seq_len(nrow(meta))) {

  f <- file.path(meta_dir, paste0("meta_", meta$compartment[i], "_", meta$contrast[i], ".csv"))
  d <- read.csv(f, stringsAsFactors = FALSE)
  sig <- d[d$meta_deg, ]

  check(nrow(sig) == meta$meta_deg_total[i],
        paste0("meta_", meta$compartment[i], "_", meta$contrast[i], ": table matches summary"),
        paste(nrow(sig), "=", meta$meta_deg_total[i]))

  check(all(!is.na(d$pval)), paste0("  no non-estimable pooled P (", meta$compartment[i], " ",
                                    meta$contrast[i], ")"))

  if (nrow(sig)) {
    check(all(abs(sig$logFC) > 1 & sig$p.adj < 0.05 & sig$consistent_direction),
          paste0("  every meta-DEG meets threshold + direction consistency (",
                 meta$compartment[i], " ", meta$contrast[i], ")"))
  }

  if (all(!is.na(d$I2))) {
    check(all(d$I2 >= 0 & d$I2 <= 100), paste0("  I2 within [0,100] (",
                                               meta$compartment[i], " ", meta$contrast[i], ")"))
  }
}

val <- read.csv(file.path(tab_dir, "meta_validation.csv"), stringsAsFactors = FALSE)

check(nrow(val) == 10, "ERCB validation computed for every primary meta", nrow(val))
check(all(val$validation_set %in% c("val_glom_array", "val_ti_array")),
      "validation uses only the independent ERCB cohorts")
check(all(val$compartment[val$validation_set == "val_glom_array"] == "Glom") &&
        all(val$compartment[val$validation_set == "val_ti_array"] == "TI"),
      "validation cohort matches compartment")

sensm <- read.csv(file.path(tab_dir, "meta_sensitivity_allcohorts.csv"), stringsAsFactors = FALSE)

check(nrow(sensm) == 10, "all-cohorts sensitivity computed for every meta", nrow(sensm))
check(all(sensm$k_max == 3 | (sensm$contrast == "MN_vs_Control" & sensm$compartment == "TI")),
      "sensitivity metas reach k = 3 (TI control has no RNA-seq donors)")

# --- step 9: enrichment ------------------------------------------------------
ora <- read.csv(file.path(tab_dir, "enrichment_summary_ora.csv"), stringsAsFactors = FALSE)
gsea <- read.csv(file.path(tab_dir, "enrichment_summary_gsea.csv"), stringsAsFactors = FALSE)
sig_e <- read.csv(file.path(tab_dir, "enrichment_expected_signals.csv"), stringsAsFactors = FALSE)

check(nrow(ora) >= 20, "ORA ran for every eligible contrast x collection", nrow(ora))
check(all(c("GO_BP", "GO_CC", "GO_MF", "Reactome", "KEGG") %in% ora$collection),
      "GO (BP/CC/MF), Reactome and KEGG all represented")
check(nrow(gsea) == 20, "GSEA ran for 10 metas x 2 collections", nrow(gsea))
check(all(gsea$collection %in% c("Hallmark", "C7")), "GSEA covers Hallmark and C7")
check(max(gsea$n_significant) > 0, "at least one GSEA collection has significant sets",
      max(gsea$n_significant))
check(nrow(sig_e) >= 10, "expected-signal audit produced results", nrow(sig_e))
check(any(sig_e$significant), "at least one expected signal is significant",
      sum(sig_e$significant))
check(!any(is.na(sig_e$best_term)), "expected-signal audit names a real term for every row")

enr_files <- list.files(enr_dir, pattern = "\\.csv$")

check(length(enr_files) >= 35, "enrichment tables written", length(enr_files))

hall <- read.csv(file.path(enr_dir, "gsea_Glom_MN_vs_Control_Hallmark.csv"), stringsAsFactors = FALSE)

check(all(c("pathway", "NES", "padj", "leadingEdge") %in% names(hall)), "GSEA tables carry NES, padj, leading edge")
check(all(abs(hall$NES) <= 5), "NES values are on a plausible normalised scale",
      paste(round(range(hall$NES), 2), collapse = " .. "))

# --- figures -----------------------------------------------------------------
pairs <- c("04a_disc_ti_array_pca_before", "04a_disc_ti_array_pca_after",
           "04a_disc_ti_array_pca_before_after")

for (stem in pairs) {
  d <- file.path(OUT, "qc")
  check(file.exists(file.path(d, paste0(stem, ".pdf"))) &&
          file.exists(file.path(d, paste0(stem, ".jpeg"))),
        paste("figure pair present:", stem))
}

n_pdf_deg <- length(list.files(fig_deg, pattern = "\\.pdf$"))
n_jpg_deg <- length(list.files(fig_deg, pattern = "\\.jpeg$"))

check(n_pdf_deg > 0 && n_pdf_deg == n_jpg_deg, "every DEG/meta figure exists as PDF + JPEG",
      paste(n_pdf_deg, "/", n_jpg_deg))

n_pdf_enr <- length(list.files(fig_enr, pattern = "\\.pdf$"))
n_jpg_enr <- length(list.files(fig_enr, pattern = "\\.jpeg$"))

check(n_pdf_enr > 0 && n_pdf_enr == n_jpg_enr, "every enrichment figure exists as PDF + JPEG",
      paste(n_pdf_enr, "/", n_jpg_enr))

cat("\n=============================================\n")
cat(sprintf("passed: %d | failed: %d\n", n_pass, n_fail))

if (n_fail == 0) {
  cat("ALL CHECKS PASSED\n")
  quit(status = 0)
} else {
  cat("CHECKS FAILED\n")
  quit(status = 1)
}
