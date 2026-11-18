#!/usr/bin/env Rscript
# Paper 1 - Phase A step 3d: freeze the QC outlier policy decided by the user (2026-10-18).
# Every sample flagged by the median-correlation screen (median - 3 * MAD) is RETAINED; the flags
# are carried forward as annotation so that the methods section can report them, and so that a
# later sensitivity analysis can exclude them without recomputing anything.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

tab_dir <- file.path(OUT, "tables")

qc <- read.csv(file.path(tab_dir, "qc_summary.csv"), stringsAsFactors = FALSE)

policy <- data.frame(set = qc$set,
                     samples = qc$samples,
                     flagged_outliers = qc$flagged_outliers,
                     flagged_pct = round(100 * qc$flagged_outliers / qc$samples, 1),
                     min_median_corr = qc$min_median_corr,
                     decision = "retain",
                     rationale = "User decision 2026-10-18: keep every sample; document flags only",
                     stringsAsFactors = FALSE)

write.csv(policy, file.path(tab_dir, "qc_outlier_policy.csv"), row.names = FALSE)

# sample-level flags, one row per flagged sample, for the methods appendix
rows <- list()

for (s in qc$set) {

  f <- file.path(OUT, "qc", paste0("03c_", s, "_sample_qc.csv"))

  if (!file.exists(f)) next

  d <- read.csv(f, stringsAsFactors = FALSE)

  d <- d[d$flagged, , drop = FALSE]

  if (nrow(d) == 0) next

  rows[[s]] <- data.frame(set = s, gsm = d$gsm, disease = d$disease,
                          median_corr = d$median_corr, decision = "retain",
                          stringsAsFactors = FALSE)
}

flags <- if (length(rows)) do.call(rbind, rows) else
  data.frame(set = character(), gsm = character(), disease = character(),
             median_corr = numeric(), decision = character())

write.csv(flags, file.path(tab_dir, "qc_flagged_samples_retained.csv"), row.names = FALSE)

cat("\n================ QC OUTLIER POLICY (FROZEN) ================\n")

print(policy[, c("set", "samples", "flagged_outliers", "flagged_pct", "decision")],
      row.names = FALSE)

cat("\ntotal flagged but retained:", nrow(flags), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 03d_qc_policy.R"),
             capture.output(print(policy[, c("set", "samples", "flagged_outliers",
                                             "flagged_pct", "decision")], row.names = FALSE)),
             paste("total flagged but retained:", nrow(flags))),
           file.path(PROJ, "logs", "qc_policy.log"))

cat("\nDone.\n")
