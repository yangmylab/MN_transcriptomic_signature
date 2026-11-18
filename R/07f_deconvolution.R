#!/usr/bin/env Rscript
# Paper 1 - scRNA arm, step 07f: NNLS deconvolution of every bulk cohort with the atlas-derived
# signature matrix. This is the CIBERSORTx-free route that unblocks Phase C.
#
# Model: mixture_s = sum_k w_k * signature_k, solved per sample with non-negative least squares
# (nnls). Signature columns and mixture are each scaled to unit sum, so weights are relative
# proportions; they are NOT absolute cell numbers (no total-cell normalisation exists here).
#
# Inputs : output/tables/signature_matrix.csv (07e), work/sets/*_expr_symbol.rds
# Outputs: output/tables/deconvolution_fractions.csv
#          output/tables/deconvolution_tests.csv
#          output/scrna/07f_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(ggplot2); library(limma) })

if (!requireNamespace("nnls", quietly = TRUE)) stop("package 'nnls' missing - run R/07z_install_scrna_pkgs.R")

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "scrna")

sig <- as.matrix(read.csv(file.path(tab_dir, "signature_matrix.csv"), row.names = 1, check.names = FALSE))

cat("signature:", nrow(sig), "genes x", ncol(sig), "cell types\n")

# --- primary signature: cell-type-specific genes only ------------------------
#
# On the full shared gene set every column correlates 0.99 with the grand mean, so NNLS has almost
# no signal to separate cell types and dumps 66-88 % of every sample on whichever column carries the
# most abundant genes ("Neutrophil": top S100A9/LYZ in the atlas). Restricting the signature to genes
# with cell-type specificity tau > 0.5 removes that shared component and restores the expected
# biology (glomerular samples become Podocyte-dominant, matching Phase B). The unfiltered fit is kept
# as an explicit sensitivity analysis instead of being silently dropped.
#
# tau MUST be computed on the CPM scale, BEFORE the unit-sum scaling below: afterwards the values are
# ~1e-5 and log2(x + 1) collapses them to ~0, which falsely marks almost every gene as specific
# (6 342 instead of 2 450 genes retained).
tau <- function(x) {
  n <- length(x)
  if (n < 2 || all(x == 0)) return(0)
  sum(1 - x / max(x)) / (n - 1)
}

gene_tau <- apply(log2(sig + 1), 1, tau)

sig <- sweep(sig, 2, colSums(sig), "/")

sig_specific <- sig[gene_tau > 0.5, , drop = FALSE]

cat("specificity filter tau > 0.5:", nrow(sig_specific), "of", nrow(sig), "genes retained\n")

deconvolve <- function(set, sig_use) {

  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)

  genes <- intersect(rownames(sig_use), rownames(ms))
  genes <- genes[rowSums(is.na(ms[genes, , drop = FALSE])) == 0]

  if (length(genes) < 500) return(NULL)

  mix <- 2^ms[genes, , drop = FALSE]

  sig_g <- sig_use[genes, , drop = FALSE]

  # The bulk matrices are on their own platform scale (RMA log2 intensity -> 2^, or voom log2 CPM
  # -> 2^) while the signature is CPM-derived, so raw NNLS is dominated by whichever cell-type
  # column has the largest magnitude at the brightest genes: every discovery array sample collapsed
  # to a single cell type at weight 1.000 (see P-28 in handoff.md).
  # Fix (the step CIBERSORT applies for the same reason): quantile-normalise the mixtures TOGETHER
  # with the signature columns, so every column of both sides shares one value distribution, then
  # unit-sum both sides and solve. Normalising the mixture to the signature's MEAN profile instead
  # makes it a permutation of that mean, which a single column then fits almost exactly (RMS 2e-4).
  # Quantile-normalise the two sides SEPARATELY. Normalising the signature columns together with the
  # mixture columns (the previous version) replaces each signature value by a cross-column quantile
  # mean dominated by the far more numerous mixture samples, which flattens the cell-type profiles
  # that the tau filter exists to preserve.
  sig_q <- pmax(normalizeQuantiles(sig_g), 0)

  mix_q <- pmax(normalizeQuantiles(mix), 0)

  sig_q <- sweep(sig_q, 2, colSums(sig_q), "/")

  mix_q <- sweep(mix_q, 2, colSums(mix_q), "/")

  fits <- lapply(seq_len(ncol(mix_q)), function(j) {
    fit <- nnls::nnls(sig_q, mix_q[, j])
    x <- fit$x
    s <- sum(x)
    # a mixture whose NNLS weights are all zero cannot be turned into proportions; flag it instead
    # of emitting NaN
    list(w = if (s > 0) x / s else rep(NA_real_, length(x)),
         zero = s <= 0,
         rms = sqrt(fit$deviance / nrow(sig_q)))
  })

  w <- t(vapply(fits, function(z) z$w, numeric(ncol(sig_g))))

  colnames(w) <- colnames(sig_g)
  df <- data.frame(sample = colnames(mix), set = set,
                   disease = ann$disease[match(colnames(mix), ann$gsm)],
                   w, check.names = FALSE, stringsAsFactors = FALSE)
  df$rms_residual <- round(vapply(fits, function(z) z$rms, numeric(1)), 5)
  df$max_weight <- round(apply(w, 1, max), 4)
  df$zero_sum <- vapply(fits, function(z) z$zero, logical(1))

  df
}

sets <- c("disc_glom_array", "disc_ti_array", "val_glom_array", "val_ti_array",
          "cross_glom_rnaseq", "cross_ti_rnaseq", "pla2r_lcm")

frac <- do.call(rbind, lapply(sets, function(s) {
  cat("  deconvolving", s, "(specific signature)\n")
  deconvolve(s, sig_specific)
}))

write.csv(frac, file.path(tab_dir, "deconvolution_fractions.csv"), row.names = FALSE)

ct_cols <- setdiff(names(frac), c("sample", "set", "disease", "rms_residual", "max_weight",
                                  "zero_sum"))

cat("\nfractions:", nrow(frac), "samples x", length(ct_cols), "cell types\n")

cat("degenerate single-cell-type fits (max weight > 0.99):",
    sum(frac$max_weight > 0.99), "of", nrow(frac), "\n")

print(tapply(frac$max_weight, frac$set, median))

# --- sensitivity: the unfiltered (all shared genes) signature ----------------

frac_all <- do.call(rbind, lapply(sets, function(s) {
  cat("  deconvolving", s, "(all-gene sensitivity)\n")
  deconvolve(s, sig)
}))

write.csv(frac_all, file.path(tab_dir, "deconvolution_fractions_allgenes.csv"), row.names = FALSE)

top_of <- function(d) {
  d$top_celltype <- ct_cols[max.col(as.matrix(d[, ct_cols]), ties.method = "first")]
  d
}

cmp <- merge(top_of(frac)[, c("sample", "set", "top_celltype", "max_weight")],
             top_of(frac_all)[, c("sample", "set", "top_celltype", "max_weight")],
             by = c("sample", "set"), suffixes = c("_specific", "_allgenes"))

write.csv(cmp, file.path(tab_dir, "deconvolution_signature_sensitivity.csv"), row.names = FALSE)

cat("\nsignature sensitivity - samples whose dominant cell type changes:",
    sum(cmp$top_celltype_specific != cmp$top_celltype_allgenes), "of", nrow(cmp), "\n")

print(table(specific = cmp$top_celltype_specific, allgenes = cmp$top_celltype_allgenes))

# --- tests: MN vs each comparator, per set and cell type ---------------------

frac$group <- ifelse(frac$disease == "UNKNOWN", NA_character_,
                     ifelse(frac$disease %in% c("Control-LD", "Control-TN"), "Control",
                            ifelse(frac$disease %in% c("MN", "FSGS", "MCD"), frac$disease,
                                   "Other")))

test_rows <- list()

for (s in unique(frac$set)) {

  d <- frac[frac$set == s & !is.na(frac$group), ]
  if (!any(d$group == "MN")) next

  for (grp in intersect(c("Control", "FSGS", "MCD", "Other"), unique(d$group))) {

    if (sum(d$group == "MN") < 3 || sum(d$group == grp) < 3) next

    for (ct in ct_cols) {

      x <- d[[ct]][d$group == "MN"]
      y <- d[[ct]][d$group == grp]

      # A cell type that is absent from a compartment deconvolves to exactly 0 in every sample, so
      # wilcox.test has zero variance and returns NaN. Those comparisons are not evidence of "no
      # difference" - they are untestable - so they are flagged instead of being silently averaged
      # into the BH correction.
      degenerate <- all(x == 0) && all(y == 0)

      test_rows[[paste(s, grp, ct)]] <- data.frame(
        set = s, celltype = ct, contrast = paste0("MN_vs_", grp),
        median_MN = round(median(x), 4), median_ref = round(median(y), 4),
        diff = round(median(x) - median(y), 4),
        p_value = if (degenerate) NA_real_ else wilcox.test(x, y)$p.value,
        n_MN = length(x), n_ref = length(y),
        testable = !degenerate,
        note = if (degenerate) "fraction 0 in both groups (cell type absent)" else "",
        stringsAsFactors = FALSE)
    }
  }
}

tests <- do.call(rbind, test_rows)
tests$p_adj <- p.adjust(tests$p_value, method = "BH")

write.csv(tests, file.path(tab_dir, "deconvolution_tests.csv"), row.names = FALSE)

cat("\nsignificant (BH < 0.05):", sum(tests$p_adj < 0.05, na.rm = TRUE), "of",
    sum(tests$testable), "testable comparisons (",
    sum(!tests$testable), "all-zero tests excluded)\n")

print(tests[which(tests$p_adj < 0.05), c("set", "contrast", "celltype", "median_MN", "median_ref", "p_adj")],
      row.names = FALSE)

# --- figures -----------------------------------------------------------------

frac_plot <- frac
frac_plot$group[is.na(frac_plot$group)] <- "Unknown"

long <- reshape(frac_plot[, c("sample", "set", "group", ct_cols)], direction = "long",
                varying = ct_cols, v.names = "fraction", timevar = "celltype",
                times = ct_cols, idvar = c("sample", "set", "group"))

p <- ggplot(long[long$set == "disc_glom_array", ], aes(sample, fraction, fill = celltype)) +
  geom_col(width = 0.9) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(title = "Deconvolved glomerular cell-type proportions (discovery array set)",
       x = NULL, y = "relative proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_plot <- function(p, stem, width = 9, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

save_plot(p, file.path(fig_dir, "07f_deconvolution_stacked"), width = 11, height = 5.5)

p <- ggplot(long, aes(group, fraction, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_point(position = position_jitter(width = 0.12), size = 1, alpha = 0.6) +
  facet_grid(celltype ~ set, scales = "free_y") +
  scale_fill_manual(values = c(Control = "grey60", FSGS = "#7FB3D5", MCD = "#F5B041",
                               MN = "#C0392B", Other = "#AAB7B8"), guide = "none") +
  labs(title = "Deconvolved proportions by disease group", x = NULL, y = "relative proportion") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p, file.path(fig_dir, "07f_deconvolution_boxplots"), width = 14, height = 10)

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 07f_deconvolution.R"),
             paste("signature: cell-type-specific genes (tau > 0.5):", nrow(sig_specific), "of", nrow(sig)),
             paste("samples deconvolved:", nrow(frac)),
             paste("significant tests:", sum(tests$p_adj < 0.05, na.rm = TRUE), "of",
                   sum(tests$testable), "testable;", sum(!tests$testable), "all-zero tests excluded")),
           file.path(PROJ, "logs", "scrna_deconvolution.log"))

cat("\nDone.\n")
