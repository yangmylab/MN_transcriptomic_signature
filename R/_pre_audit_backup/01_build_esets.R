#!/usr/bin/env Rscript
# Paper 1 - Phase A step 2: read every GEO series matrix into an ExpressionSet, persist as RDS,
# and export per-sample phenotype tables. Reports dimension + whether expression values exist.
local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })
suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })

proc_dir <- file.path(DATA, "processed"); dir.create(proc_dir, showWarnings = FALSE)
meta_dir <- file.path(PROJ, "metadata");   dir.create(meta_dir, showWarnings = FALSE)

sm <- list.files(DATA, pattern = "series_matrix\\.txt\\.gz$", recursive = TRUE, full.names = TRUE)
cat("Found", length(sm), "series matrices\n")
summary_rows <- list()

for (f in sm) {
  tag <- sub("_series_matrix\\.txt\\.gz$", "", basename(f))
  cat("\n=====", tag, "=====\n")
  eset <- tryCatch(suppressWarnings(getGEO(filename = f)), error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(eset)) next
  if (is(eset, "ExpressionSet")) {
    ex <- exprs(eset); pd <- pData(eset)
    cat("  platform:", annotation(eset), " features:", nrow(ex), " samples:", ncol(ex), "\n")
    saveRDS(eset, file.path(proc_dir, paste0(tag, ".rds")))
    write.csv(pd, file.path(meta_dir, paste0("pdata_", tag, ".csv")), row.names = TRUE)
    summary_rows[[tag]] <- data.frame(tag = tag, platform = annotation(eset),
                                      features = nrow(ex), samples = ncol(ex),
                                      has_expr = nrow(ex) > 0 && any(is.finite(ex)),
                                      class = class(eset)[1])
  } else {
    cat("  not an ExpressionSet:", class(eset)[1], "\n")
  }
}

if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows)
  write.csv(s, file.path(meta_dir, "series_matrix_summary.csv"), row.names = FALSE)
  cat("\n================ SUMMARY ================\n"); print(s, row.names = FALSE)
}
cat("\nDone.\n")
