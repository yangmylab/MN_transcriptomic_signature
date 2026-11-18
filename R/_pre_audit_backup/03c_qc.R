#!/usr/bin/env Rscript
# Paper 1 - Phase A step 3c: per-cohort QC on the assembled sets - PCA, sample clustering and an
# outlier screen - with every plot saved as PDF (vector) and JPEG, and every sample decision
# written to a table for the methods section.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(ggplot2); library(pheatmap) })

set_dir <- file.path(PROJ, "work", "sets")
qc_dir <- file.path(OUT, "qc")
tab_dir <- file.path(OUT, "tables")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

sets <- sub("_expr\\.rds$", "", list.files(set_dir, pattern = "_expr\\.rds$"))

cat("QC on", length(sets), "cohort sets:", paste(sets, collapse = ", "), "\n")

qc_rows <- list()

save_plot <- function(p, stem) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = 7, height = 5.5)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = 7, height = 5.5)
}

for (s in sets) {
  cat("\n=====", s, "=====\n")
  ms <- readRDS(file.path(set_dir, paste0(s, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(s, "_samples.csv")), stringsAsFactors = FALSE)

  # genes absent from one of the platforms of a merged release carry NA values, so keep only genes
  # measured on every sample before any variance / correlation computation
  keep_g <- rowSums(is.na(ms)) == 0
  cat("  complete genes:", sum(keep_g), "of", nrow(ms), "\n")
  ms <- ms[keep_g, , drop = FALSE]

  # top-variable genes drive the PCA
  v <- apply(ms, 1, sd, na.rm = TRUE)
  top <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
  pc <- prcomp(t(ms[top, , drop = FALSE]), scale. = TRUE)

  ve <- round(100 * summary(pc)$importance[2, 1:2], 1)
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], disease = ann$disease, gsm = ann$gsm)

  p <- ggplot(df, aes(PC1, PC2, colour = disease)) +
    geom_point(size = 2.2, alpha = 0.85) +
    labs(title = paste0(s, " - PCA (top 2000 variable genes)"),
         x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_blank())

  save_plot(p, file.path(qc_dir, paste0("03c_", s, "_pca")))

  # sample-sample correlation screen: flag samples whose median correlation to the cohort is low
  cc <- cor(ms[top, , drop = FALSE], use = "pairwise.complete.obs")
  med <- apply(cc, 2, median, na.rm = TRUE)
  out <- names(med)[med < (median(med) - 3 * mad(med))]

  cat("  samples:", ncol(ms), "| median-corr range:",
      paste(round(range(med), 3), collapse = " .. "),
      "| flagged:", length(out), if (length(out)) paste(out, collapse = ",") else "", "\n")

  qc_rows[[s]] <- data.frame(set = s, samples = ncol(ms), genes = nrow(ms),
                             pc1_var = ve[1], pc2_var = ve[2],
                             min_median_corr = round(min(med), 3),
                             flagged_outliers = length(out),
                             flagged_samples = paste(out, collapse = ","),
                             stringsAsFactors = FALSE)

  write.csv(data.frame(gsm = colnames(ms), disease = ann$disease,
                       median_corr = round(med, 4), flagged = colnames(ms) %in% out),
            file.path(qc_dir, paste0("03c_", s, "_sample_qc.csv")), row.names = FALSE)
}

qc <- do.call(rbind, qc_rows)

write.csv(qc, file.path(tab_dir, "qc_summary.csv"), row.names = FALSE)

cat("\n================ QC SUMMARY ================\n")

print(qc, row.names = FALSE)

cat("\ntotal flagged samples:", sum(qc$flagged_outliers), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 03c_qc.R"),
             capture.output(print(qc, row.names = FALSE))),
           file.path(PROJ, "logs", "qc.log"))

cat("\nDone.\n")
