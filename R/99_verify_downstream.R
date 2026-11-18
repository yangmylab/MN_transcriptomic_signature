#!/usr/bin/env Rscript
# Paper 1 - downstream verifier (separability gate + exploratory subtyping).
# Each assertion is a claim these two analyses make; the script exits 0 only if all hold.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")
gate_fig <- file.path(OUT, "models_gate")
sub_fig <- file.path(OUT, "subtypes")

n_pass <- 0L
n_fail <- 0L

check <- function(ok, label, detail = "") {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat(sprintf("[PASS] %-72s %s\n", label, detail))
  } else {
    n_fail <<- n_fail + 1L
    cat(sprintf("[FAIL] %-72s %s\n", label, detail))
  }
}

cat("=========== DOWNSTREAM VERIFICATION ===========\n\n")

# --- gate: nested CV and external transfer ----------------------------------
nested <- read.csv(file.path(tab_dir, "separability_nestedcv.csv"), stringsAsFactors = FALSE)
ext <- read.csv(file.path(tab_dir, "separability_external.csv"), stringsAsFactors = FALSE)

check(nrow(nested) == 6, "nested CV covers 6 feature sizes", nrow(nested))
check(all(nested$nested_cv_auc_mean > 0.8), "nested CV AUC > 0.8 for every K",
      paste(round(range(nested$nested_cv_auc_mean), 3), collapse = " .. "))
check(nrow(ext) == 12, "external transfer covers 6 K x 2 cohorts", nrow(ext))

ercb <- ext[ext$cohort == "ERCB (external)", ]
# Two corrections are encoded here, and both must stay:
#   F07/P-58 (2026-10-20) removed the cohort-internal z-scoring (transductive: the ERCB AUC leaned on
#     the validation cohort's own distribution).
#   F62 (2026-11-18) then showed that training-set mu/sd alone is NOT sufficient across array
#     platforms (discovery HuGene-2.1-ST vs ERCB HG-U133): the typical ERCB gene sat at |z| ~ 4
#     instead of ~0.6, the linear predictor saturated and the "transfer" measured platform shift -
#     per-diagnosis AUCs of 0.05-0.19 for RPGN/FSGS_MCD and a K-dependent AUC of 0.64-0.77.
# 05a/05b now quantile-map each validation cohort onto the discovery reference distribution
# (rank-preserving, label-blind) BEFORE applying training statistics. The external AUC is 0.861-0.922
# across feature sizes and 05b's independent re-fit agrees to within 0.01 (reconciliation line in
# logs/gate_robustness.out). The bound therefore returns to 0.8: it now tests a valid transfer.
check(all(ercb$auc > 0.8), "ERCB AUC > 0.8 for every K (harmonised transfer, F07 + F62)",
      paste(round(range(ercb$auc), 3), collapse = " .. "))
check(all(ercb$ci_low > 0.5), "every ERCB AUC CI excludes chance", min(ercb$ci_low))
check(all(grepl("independent cohort", ext$cohort[ext$cohort != "ERCB (external)"])),
      "RNA-seq cohort is labelled independent (no patient linkage)")

sec <- read.csv(file.path(tab_dir, "separability_secondary.csv"), stringsAsFactors = FALSE)
# MN vs ALL non-MN including healthy donors. Under the harmonised transfer (F62) this reaches 0.90, so
# the gate is not driven by a disease-vs-normal axis alone; the assertion keeps the row present and
# requires it to clear chance by a wide margin.
check(nrow(sec) == 1 && sec$auc > 0.75,
      "secondary endpoint (incl. donors) reported and clearly above chance (0.90 harmonised, F62)",
      sec$auc)

conf <- read.csv(file.path(tab_dir, "separability_confusion.csv"), stringsAsFactors = FALSE)
check(all(c("truth", "predicted", "K") %in% names(conf)), "multiclass confusion stored")
check(nrow(conf[conf$K == "all", ]) == 170, "ERCB non-control samples all classified",
      nrow(conf[conf$K == "all", ]))
check(all(conf$truth %in% c("MN", "FSGS", "MCD", "Other")), "multiclass uses the 4-class vocabulary")

# --- gate: robustness -------------------------------------------------------
rob <- read.csv(file.path(tab_dir, "gate_robustness.csv"), stringsAsFactors = FALSE)
val <- function(lbl) rob$value[rob$check == lbl]

check(nrow(rob) == 7, "all seven robustness checks reported", nrow(rob))
check(abs(val("permutation null: median AUC") - 0.5) < 0.05,
      "permutation null centred on chance", val("permutation null: median AUC"))
check(val("empirical p") <= 0.05, "observed AUC beats every shuffle", val("empirical p"))
check(val("ERCB QC-flagged samples removed") >= val("reference (discovery -> ERCB)"),
      "signal survives removal of QC-flagged ERCB samples",
      paste(val("reference (discovery -> ERCB)"), "->", val("ERCB QC-flagged samples removed")))
check(val("reverse direction (ERCB -> discovery)") > 0.7, "reverse-direction transfer holds",
      val("reverse direction (ERCB -> discovery)"))

comp <- read.csv(file.path(tab_dir, "gate_auc_by_comparator.csv"), stringsAsFactors = FALSE)
check(nrow(comp) == 9, "per-diagnosis AUCs cover all ERCB comparators", nrow(comp))
check(all(comp$auc > 0.6), "MN is discriminated from every single comparator", min(comp$auc))

# Gate stability across fold partitions (05c, added 2026-11-18). F62 showed the *unharmonised* gate
# swung 0.638-0.679 across twelve partition draws with inverted comparators; the harmonised gate must
# be demonstrably stable, or the reported 0.922 is one lucky draw.
stab <- read.csv(file.path(tab_dir, "05c_gate_seed_stability.csv"), stringsAsFactors = FALSE)
check(nrow(stab) == 20, "gate stability: 10 partition draws at 2 feature counts", nrow(stab))
check(all(table(stab$K) >= 8), "gate stability: at least 8 draws per feature count",
      paste(names(table(stab$K)), table(stab$K), sep = "=", collapse = " "))
check(all(stab$ercb_auc > 0.8), "gate stability: every draw clears 0.8",
      paste(range(stab$ercb_auc), collapse = " .. "))
rng <- tapply(stab$ercb_auc, stab$K, function(x) diff(range(x)))
check(all(rng < 0.05), "gate stability: AUC range < 0.05 within each feature count",
      paste(round(rng, 4), collapse = " / "))
check(all(stab$ci_low > 0.5), "gate stability: every DeLong CI excludes chance", min(stab$ci_low))

# --- subtyping --------------------------------------------------------------
stab <- read.csv(file.path(tab_dir, "subtype_stability.csv"), stringsAsFactors = FALSE)

check(nrow(stab) == 10, "stability computed for k = 2..6 on both platforms", nrow(stab))
check(any(stab$stable[grepl("array", stab$platform)]),
      "array cohort has at least one stable partition",
      paste(stab$k[stab$stable & grepl("array", stab$platform)], collapse = ","))
check(!any(stab$stable[grepl("rnaseq", stab$platform)]),
      "RNA-seq cohort has NO stable partition (negative result, must stay documented)")

conc <- read.csv(file.path(tab_dir, "subtype_signature_concordance.csv"), stringsAsFactors = FALSE)
cv <- function(lbl) conc$value[conc$metric == lbl]

check(abs(cv("Jaccard of defining genes") - 0.0068) < 0.01,
      "cross-platform signature overlap is reported as near-zero",
      cv("Jaccard of defining genes"))
check(cv("max cross-platform centroid correlation") < 0.7,
      "cross-platform centroid correlation is reported as moderate at best",
      cv("max cross-platform centroid correlation"))

qc <- read.csv(file.path(tab_dir, "subtype_qc_check.csv"), stringsAsFactors = FALSE)
check(min(qc$p_median_corr) < 0.05, "array subtype split is flagged as quality-associated",
      signif(min(qc$p_median_corr), 3))

cl <- read.csv(file.path(tab_dir, "subtype_ercb_assignment.csv"), stringsAsFactors = FALSE)
check(nrow(cl) == 21, "all ERCB MN samples assigned", nrow(cl))
check(!any(is.na(cl$array_framework)), "no ERCB sample left unassigned")

# --- figures ----------------------------------------------------------------
for (d in c(gate_fig, sub_fig)) {
  pdfs <- list.files(d, pattern = "\\.pdf$")
  jpegs <- list.files(d, pattern = "\\.jpeg$")
  check(length(pdfs) > 0 && length(pdfs) == length(jpegs),
        paste("every figure in", basename(d), "exists as PDF + JPEG"),
        paste(length(pdfs), "/", length(jpegs)))
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
