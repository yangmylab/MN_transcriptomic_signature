#!/usr/bin/env Rscript
# Track B verifier - the panel benchmark's own claims (2026-11-18).
# Each assertion is something 10c asserts about the attrition ladder; the script exits 0 only if all hold.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")
meta_dir <- file.path(PROJ, "metadata")
fig_dir <- file.path(OUT, "models_trackb")

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

cat("=== Track B verifier ===\n\n")

# --- register ---------------------------------------------------------------
reg <- read.csv(file.path(meta_dir, "trackB_panel_register.csv"), stringsAsFactors = FALSE)

check(nrow(reg) >= 10, "panel register holds at least 10 published panels (B15 bar)", nrow(reg))
# Identifier rule: a panel must be traceable. PMID preferred; a DOI alone is acceptable (P06's PMID
# was not retrievable from the source page, its DOI was). The gap is recorded in the row's note.
check(all(nzchar(reg$pmid) | nzchar(reg$doi)), "every panel carries a PMID or a DOI", "")
check(all(nzchar(reg$doi)), "every panel carries a DOI", "")
check(all(nzchar(reg$genes)), "every panel carries an explicit gene list", "")
check(all(nzchar(reg$reported_auc_text)), "every panel carries its reported performance text", "")
check(!any(duplicated(reg$panel_id)), "panel ids are unique", paste(reg$panel_id, collapse = " "))

# provenance discipline: the 'as published' description must state dedup/standardisation status
check(all(reg$deduplication_stated %in% c("yes", "no")),
      "de-duplication status recorded for every panel (yes/no)", "")
check(all(reg$standardisation_stated %in% c("yes", "no")),
      "standardisation status recorded for every panel (yes/no)", "")

# inclusion rule: a multivariable panel, i.e. at least two named genes
check(all(reg$n_genes_reported >= 2), "every included panel has at least 2 genes (panel, not single marker)",
      paste(range(reg$n_genes_reported), collapse = " .. "))

# exclusions must be logged with a reason, so the register cannot be padded or trimmed silently
exc_file <- file.path(meta_dir, "trackB_exclusions.csv")
check(file.exists(exc_file), "exclusion log exists", "")
if (file.exists(exc_file)) {
  exc <- read.csv(exc_file, stringsAsFactors = FALSE)
  check(all(nzchar(exc$reason_excluded)), "every excluded candidate has a stated reason", nrow(exc))
  check(!any(exc$first_author_year == ""), "every excluded candidate is dated", "")
}

# --- attrition ladder -------------------------------------------------------
att <- read.csv(file.path(tab_dir, "10c_attrition.csv"), stringsAsFactors = FALSE)

rungs_expect <- c("R1", "R2", "R3", "R4", "R5", "R5r", "R6")
check(setequal(unique(att$rung), rungs_expect), "all seven rungs present",
      paste(sort(unique(att$rung)), collapse = ","))
check(nrow(att) == nrow(reg) * length(rungs_expect),
      "ladder is complete: panels x rungs, no row missing",
      paste(nrow(att), "rows for", nrow(reg), "panels"))
check(!any(duplicated(att[c("panel_id", "rung")])), "no duplicated panel x rung", "")

auc_rows <- att[att$metric == "auc", ]
check(all(auc_rows$auc >= 0 & auc_rows$auc <= 1), "every AUC lies in [0,1]",
      paste(range(auc_rows$auc), collapse = " .. "))
check(all(!is.na(auc_rows$auc)), "no missing AUC in any computed rung", "")
check(all(auc_rows$ci_low <= auc_rows$auc & auc_rows$auc <= auc_rows$ci_high),
      "every inferred CI brackets its AUC", "")

p_rows <- att[att$metric == "empirical_p", ]
check(all(p_rows$auc > 0 & p_rows$auc <= 1), "every refit-null p-value lies in (0,1]",
      paste(range(p_rows$auc), collapse = " .. "))

# genes actually used can never exceed the published list, and at least one must be measurable
check(all(att$n_genes_used <= att$n_genes_published), "no panel uses more genes than it published", "")
check(all(att$n_genes_used >= 1), "every panel has at least one measurable gene", "")

# the note on R2 must state the de-duplication outcome, and it must be truthful about n
check(all(grepl("collapsed", att$note[att$rung == "R2"])),
      "R2 records what de-duplication did to n", "")

# --- summary ----------------------------------------------------------------
summ <- read.csv(file.path(tab_dir, "10c_panel_summary.csv"), stringsAsFactors = FALSE)

check(nrow(summ) == nrow(reg), "summary covers every panel once", nrow(summ))
check(all(c("R1", "R3", "R5", "perm_p") %in% names(summ)), "summary carries R1/R3/R5 and the null p", "")
check(all(abs(summ$drop_R1_R5 - (summ$R1 - summ$R5)) < 1e-6),
      "drop_R1_R5 is arithmetically consistent with the rungs", "")

# --- panel-size frontier (10d) ----------------------------------------------
fr <- read.csv(file.path(tab_dir, "10d_panel_size_frontier.csv"), stringsAsFactors = FALSE)

check(nrow(fr) == 9, "panel-size frontier covers 8 K values plus the locked reference", nrow(fr))
tk <- fr[is.na(fr$locked_panel), ]
check(nrow(tk) == 8 && !any(duplicated(tk$K)), "top-K curve has one row per K", paste(sort(tk$K), collapse = ","))
check(any(!is.na(fr$locked_panel)), "locked four-gene panel appears as a reference row", "")

# the locked panel must match the tables it was read from (09a CV, 09b external), not a transcription
cv_l <- read.csv(file.path(tab_dir, "09a_panel_cv_summary.csv"), stringsAsFactors = FALSE)
ext_l <- read.csv(file.path(tab_dir, "09b_panel_external.csv"), stringsAsFactors = FALSE)
lock <- fr[!is.na(fr$locked_panel), ]
check(abs(lock$cv_auc_mean - cv_l$mean_auc[cv_l$method == "lasso"]) < 1e-6,
      "locked panel CV AUC reproduces 09a exactly", "")
check(abs(lock$ercb_auc - ext_l$auc[ext_l$cohort == "ERCB external"]) < 1e-6,
      "locked panel ERCB AUC reproduces 09b exactly", "")

# the two claims the frontier is used to support
check(tk$cv_auc_mean[tk$K == 12] > tk$cv_auc_mean[tk$K == 50],
      "internal CV peaks before K = 50 (larger panels overfit)",
      paste(round(tk$cv_auc_mean[tk$K == 12], 3), ">", round(tk$cv_auc_mean[tk$K == 50], 3)))
check(tk$ercb_auc[tk$K == 50] < tk$ercb_auc[tk$K == 4],
      "external ERCB degrades from K = 4 to K = 50",
      paste(round(tk$ercb_auc[tk$K == 4], 3), "->", round(tk$ercb_auc[tk$K == 50], 3)))
check(abs(lock$cv_auc_mean - tk$cv_auc_mean[tk$K == 4]) < 0.05,
      "three-method locked panel is within 0.05 CV AUC of the best top-4",
      sprintf("%.4f vs %.4f", lock$cv_auc_mean, tk$cv_auc_mean[tk$K == 4]))

# --- figure + 0-byte guard ---------------------------------------------------
for (stem in c("10c_attrition", "10d_panel_size_frontier")) {
  for (ext in c("pdf", "jpeg")) {
    f <- file.path(fig_dir, paste0(stem, ".", ext))
    check(file.exists(f) && file.size(f) > 0,
          sprintf("%s exists as %s", stem, toupper(ext)),
          if (file.exists(f)) paste(round(file.size(f) / 1024), "kB") else "missing")
  }
}

cat(sprintf("\n=========== %d passed, %d failed ===========\n", n_pass, n_fail))
if (n_fail > 0) {
  cat("TRACK B ASSERTIONS FAILED\n")
  quit(status = 1)
}

cat("ALL TRACK B ASSERTIONS PASSED\n")
