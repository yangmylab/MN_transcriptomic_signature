#!/usr/bin/env Rscript
# Paper 1 - deferred provenance test for MYOZ2 (handoff P1b), NOT a Paper 1 pipeline method.
#
# Question. MYOZ2 is the strongest single panel gene in bulk (ERCB AUC 0.903) yet is detected in
# almost no kidney single-cell/single-nucleus profile. Does its bulk signal track a perivascular /
# smooth-muscle / mesangial programme (ACTA2, MYH11, TAGLN, PDGFRB), a podocyte programme
# (NPHS1/NPHS2), a fibroblast programme (COL1A1), or none of them? The answer decides whether
# "cellular source unresolved" can be replaced by a tested statement.
#
# Design (pre-specified in this header, before the run)
#   Bulk, four cohorts (disc array, ERCB array, RNA-seq, LCM) x two sample sets (all non-UNKNOWN
#   samples; MN only). Spearman rho(MYOZ2, every gene) gives an empirical null: the percentile of
#   each marker is its rank among all genes, and BH q is computed over all genes. A partial
#   Spearman controlling for each sample's median expression (rank-residual method) is reported as
#   a sensitivity analysis for global expression/composition shifts. Composite scores (mean z of
#   the member genes within the analysis set) avoid depending on one marker.
#   Verdict rule, fixed in advance:
#     "tracks"      = >= 3 of the 4 SMC/perivascular markers with rho >= 0.5 and q < 0.05
#                     in >= 3 of the 4 cohorts (same rule separately for each sample set);
#     "does not track" = fewer than 2 such markers in >= 3 of the 4 cohorts;
#     anything else = "cohort-dependent".
#   snRNA subset (GSE271188) and scRNA subset (GSE171458 + GSE241302), read from the per-sample
#   QC objects because the atlas HVG assay (2,019 genes) does not carry MYOZ2:
#     - MYOZ2 detection fraction per cell type, and enrichment of MYOZ2+ cells over the cell-type
#       background (exact binomial test);
#     - co-detection with each marker (Fisher exact; Haldane-corrected OR when a cell is zero);
#     - donor-level pseudobulk CPM correlation across snRNA donors.
#   With ~0.1-0.8 % MYOZ2+ cells per sample the cell-level tests are explicitly low-powered; the
#   script reports counts, not just p-values, so "no evidence" cannot be read as "evidence of no".
#
# Inputs (read-only)
#   work/sets/{disc_glom_array,val_glom_array,cross_glom_rnaseq,pla2r_lcm}_expr_symbol.rds
#   work/sets/*_samples.csv
#   work/scrna/qc/*.rds                 per-sample Seurat objects (full gene universe)
#   work/scrna/atlas_cells.csv          atlas cell -> celltype/donor/accession key
#
# Outputs
#   output/tables/09g_moyz2_bulk_correlations.csv   every target/reference gene, all cohorts/sets
#   output/tables/09g_moyz2_bulk_summary.csv        per-cohort composite correlations + verdict
#   output/tables/09g_moyz2_panel_reference.csv     each panel gene vs each composite (context)
#   output/tables/09g_moyz2_cell_detection.csv      MYOZ2 by platform x cell type + enrichment
#   output/tables/09g_moyz2_cell_codetection.csv    MYOZ2 vs each marker (Fisher)
#   output/tables/09g_moyz2_donor_correlation.csv   snRNA donor-level pseudobulk correlations
#   output/models_panel/09g_moyz2_{bulk_correlations,bulk_scatter,cell_detection}.{pdf,jpeg}

local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")])
  here <- if (length(f)) dirname(normalizePath(f[1])) else getwd()
  source(file.path(dirname(here), "R", "_lib.R"))
})

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_panel")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, stem, width = 8, height = 5.2) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  # explicit grDevices device: the string "jpeg" routes through ragg 1.5.2, which is built
  # against an older graphics API and fails on this R 4.6.1 ("Graphics API version mismatch")
  ggsave(paste0(stem, ".jpeg"), plot = p, device = grDevices::jpeg, quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

TARGET <- "MYOZ2"
SMC <- c("ACTA2", "MYH11", "TAGLN", "PDGFRB")           # vascular smooth muscle / pericyte / mesangial
POD <- c("NPHS1", "NPHS2")                              # podocyte
FIB <- "COL1A1"                                         # fibroblast
REF <- c(POD, FIB, "RGS5", "PLVAP")                     # pericyte, endothelium
PANEL_GENES <- c("SLAMF8", "TRPC6", "MYOZ2", "PPP1R3C")

role_of <- function(g) {
  ifelse(g %in% SMC, "SMC/perivascular",
  ifelse(g %in% POD, "podocyte",
  ifelse(g == FIB, "fibroblast",
  ifelse(g == "RGS5", "pericyte",
  ifelse(g == "PLVAP", "endothelium",
  ifelse(g %in% PANEL_GENES, "panel gene", "other"))))))
}

sets <- c("disc_glom_array", "val_glom_array", "cross_glom_rnaseq", "pla2r_lcm")
set_label <- c(disc_glom_array = "Discovery array", val_glom_array = "ERCB array",
               cross_glom_rnaseq = "RNA-seq", pla2r_lcm = "LCM")

load_set <- function(s) {
  ms <- readRDS(file.path(set_dir, paste0(s, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(s, "_samples.csv")), stringsAsFactors = FALSE)
  stopifnot(ncol(ms) == nrow(ann))
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  keep <- !is.na(grp) & grp != "UNKNOWN"
  list(ms = ms[, keep, drop = FALSE], grp = grp[keep])
}

# --- helpers ----------------------------------------------------------------

# Spearman rho and t-approximation p for one gene vector vs MYOZ2.
spearman_one <- function(v, x) {
  ok <- is.finite(v) & is.finite(x)
  n <- sum(ok)
  if (n < 5 || sd(v[ok]) == 0 || sd(x[ok]) == 0) return(c(n = n, rho = NA, p = NA))
  r <- cor(v[ok], x[ok], method = "spearman")
  tstat <- r * sqrt((n - 2) / max(1e-12, 1 - r^2))
  c(n = n, rho = r, p = 2 * pt(-abs(tstat), df = n - 2))
}

fisher_ci <- function(r, n) {
  if (!is.finite(r) || n < 4 || abs(r) >= 1) return(c(NA, NA))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(c(z - 1.96 * se, z + 1.96 * se))
}

# Partial Spearman: residualise ranks of gene and MYOZ2 on ranks of the covariate.
partial_spearman <- function(v, x, z) {
  ok <- is.finite(v) & is.finite(x) & is.finite(z)
  if (sum(ok) < 10) return(NA_real_)
  rv <- rank(v[ok]); rx <- rank(x[ok]); rz <- rank(z[ok])
  rxy <- cor(rv, rx); rxz <- cor(rv, rz); ryz <- cor(rx, rz)
  den <- sqrt(max(1e-12, (1 - rxz^2) * (1 - ryz^2)))
  (rxy - rxz * ryz) / den
}

composite <- function(ms, genes) {
  g <- intersect(genes, rownames(ms))
  if (length(g) < 2) return(rep(NA_real_, ncol(ms)))
  z <- t(scale(t(ms[g, , drop = FALSE])))
  # an NA/constant gene must not be turned into a zero expression value: keep it missing
  z[!is.finite(z)] <- NA_real_
  v <- colMeans(z, na.rm = TRUE)
  v[is.nan(v)] <- NA_real_
  v
}

# --- bulk correlations ------------------------------------------------------

bulk_rows <- list()
bulk_summary <- list()
panel_ref <- list()

for (s in sets) {

  o <- load_set(s)
  stopifnot(TARGET %in% rownames(o$ms))
  stopifnot(all(SMC %in% rownames(o$ms)))

  for (ss in c("all", "MN")) {

    sel <- if (ss == "all") rep(TRUE, length(o$grp)) else o$grp == "MN"
    ms <- o$ms[, sel, drop = FALSE]
    x <- ms[TARGET, ]
    if (sum(is.finite(x)) < 5 || sd(x, na.rm = TRUE) == 0) {
      warning(sprintf("%s / %s: MYOZ2 has no usable variance (%d finite)", s, ss, sum(is.finite(x))))
      next
    }

    # empirical null over every gene
    rho_all <- suppressWarnings(as.numeric(cor(t(ms), x, method = "spearman",
                                               use = "pairwise.complete.obs")))
    names(rho_all) <- rownames(ms)
    n_all <- apply(ms, 1, function(v) sum(is.finite(v) & is.finite(x)))
    t_all <- rho_all * sqrt((n_all - 2) / pmax(1e-12, 1 - rho_all^2))
    p_all <- 2 * pt(-abs(t_all), df = n_all - 2)
    q_all <- p.adjust(p_all, method = "BH")
    pct_all <- 100 * rank(abs(rho_all), na.last = "keep", ties.method = "average") /
      sum(is.finite(rho_all))
    med_expr <- apply(ms, 2, median, na.rm = TRUE)

    targets <- intersect(c(TARGET, SMC, REF, PANEL_GENES), rownames(ms))
    for (g in targets) {
      i <- match(g, rownames(ms))
      ci <- fisher_ci(rho_all[i], n_all[i])
      bulk_rows[[paste(s, ss, g)]] <- data.frame(
        cohort = s, cohort_label = set_label[[s]], sample_set = ss,
        n_samples = sum(is.finite(x)), gene = g, role = role_of(g),
        n = n_all[i], rho = round(rho_all[i], 4),
        ci_low = round(ci[1], 4), ci_high = round(ci[2], 4),
        p = signif(p_all[i], 4), q_BH_all_genes = signif(q_all[i], 4),
        percentile_abs_rho = round(pct_all[i], 2),
        partial_rho_adj_sample_median =
          if (g == TARGET) NA_real_ else round(partial_spearman(ms[g, ], x, med_expr), 4),
        stringsAsFactors = FALSE)
    }

    # composites (MYOZ2 and, for context, every panel gene)
    for (pg in intersect(PANEL_GENES, rownames(ms))) {
      y <- ms[pg, ]
      for (comp_name in c("SMC/perivascular", "podocyte", "fibroblast")) {
        genes <- switch(comp_name, "SMC/perivascular" = SMC, "podocyte" = POD, "fibroblast" = FIB)
        comp <- composite(ms, genes)
        if (all(!is.finite(comp))) next
        st <- spearman_one(y, comp)
        ci <- fisher_ci(st[["rho"]], st[["n"]])
        panel_ref[[paste(s, ss, pg, comp_name)]] <- data.frame(
          cohort = s, cohort_label = set_label[[s]], sample_set = ss, gene = pg,
          composite = comp_name, n = st[["n"]], rho = round(st[["rho"]], 4),
          ci_low = round(ci[1], 4), ci_high = round(ci[2], 4),
          p = signif(st[["p"]], 4),
          partial_rho_adj_sample_median = round(partial_spearman(y, comp, med_expr), 4),
          stringsAsFactors = FALSE)
        if (pg == TARGET) {
          bulk_summary[[paste(s, ss, comp_name)]] <- data.frame(
            cohort = s, cohort_label = set_label[[s]], sample_set = ss,
            n_samples = st[["n"]], composite = comp_name,
            rho = round(st[["rho"]], 4), p = signif(st[["p"]], 4),
            partial_rho = round(partial_spearman(y, comp, med_expr), 4),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
}

bulk <- do.call(rbind, bulk_rows)
ref <- do.call(rbind, panel_ref)
summ <- do.call(rbind, bulk_summary)

# --- verdict, from the counts (never hard-coded) -----------------------------

verdict_for <- function(d) {
  smc <- d[d$role == "SMC/perivascular", ]
  strong <- sum(smc$rho >= 0.5 & smc$q_BH_all_genes < 0.05, na.rm = TRUE)
  if (strong >= 3) "tracks SMC/perivascular programme"
  else if (strong < 2) "does not track SMC/perivascular programme"
  else "cohort-dependent"
}

for (ss in c("all", "MN")) {
  d <- bulk[bulk$sample_set == ss, ]
  per <- do.call(rbind, lapply(sets, function(s) {
    dd <- d[d$cohort == s, ]
    data.frame(cohort = s, sample_set = ss,
               n_SMC_markers_rho_ge_0.5_q_lt_0.05 =
                 sum(dd$rho[dd$role == "SMC/perivascular"] >= 0.5 &
                       dd$q_BH_all_genes[dd$role == "SMC/perivascular"] < 0.05, na.rm = TRUE),
               verdict = verdict_for(dd), stringsAsFactors = FALSE)
  }))
  cat(sprintf("\n--- verdict, sample set = %s ---\n", ss)); print(per, row.names = FALSE)
  # global rule: the cohort verdict must be consistent in >= 3 of 4 cohorts
  v <- table(per$verdict)
  global <- names(v)[which.max(v)]
  cat(sprintf("  global: %s in %d/4 cohorts\n", global, max(v)))
}

summ <- merge(summ, {
  v <- do.call(rbind, lapply(c("all", "MN"), function(ss) {
    d <- bulk[bulk$sample_set == ss, ]
    do.call(rbind, lapply(sets, function(s) data.frame(
      cohort = s, sample_set = ss, verdict = verdict_for(d[d$cohort == s, ]),
      stringsAsFactors = FALSE)))
  }))
  v
}, by = c("cohort", "sample_set"), all.x = TRUE, sort = FALSE)

write.csv(bulk, file.path(tab_dir, "09g_moyz2_bulk_correlations.csv"), row.names = FALSE)
write.csv(summ, file.path(tab_dir, "09g_moyz2_bulk_summary.csv"), row.names = FALSE)
write.csv(ref, file.path(tab_dir, "09g_moyz2_panel_reference.csv"), row.names = FALSE)

cat("\n--- bulk:", nrow(bulk), "gene x cohort x set rows ---\n")
print(bulk[bulk$gene == TARGET, c("cohort", "sample_set", "n", "gene", "rho", "q_BH_all_genes",
                                  "percentile_abs_rho", "partial_rho_adj_sample_median")],
      row.names = FALSE)

# --- cell-level: snRNA subset + scRNA comparison ----------------------------

cells <- read.csv(file.path(PROJ, "work", "scrna", "atlas_cells.csv"), stringsAsFactors = FALSE)
cells$platform <- ifelse(cells$accession == "GSE271188", "snRNA", "scRNA")

qc_files <- list.files(file.path(PROJ, "work", "scrna", "qc"), pattern = "\\.rds$", full.names = TRUE)
qc_files <- qc_files[!grepl("_qc\\.rds$", qc_files)]
stopifnot(length(qc_files) > 0)

genes_use <- c(TARGET, SMC, REF)
cache <- file.path(PROJ, "work", "scrna", "09g_moyz2_cell_matrix.rds")
if (file.exists(cache) && !"--refresh" %in% commandArgs(FALSE)) {
  cl <- readRDS(cache)
  cat("  cell matrix from resume cache:", nrow(cl), "cells\n")
} else {
  cell_rows <- list()
  for (f in qc_files) {
    obj <- readRDS(f)
    uid <- paste0(sub("\\.rds$", "", basename(f)), "__", colnames(obj))
    idx <- match(uid, cells$cell)
    keep <- !is.na(idx)
    if (!any(keep)) { rm(obj); invisible(gc()); next }
    meta <- cells[idx[keep], ]
    g <- intersect(genes_use, rownames(obj))
    if (!TARGET %in% g) { rm(obj); invisible(gc()); next }
    cm <- LayerData(obj, assay = "RNA", layer = "counts")
    cm <- as.matrix(cm[g, keep, drop = FALSE])
    df <- data.frame(
      platform = meta$platform, celltype = meta$celltype, donor = meta$donor,
      disease_group = meta$disease_group,
      total_counts = obj$nCount_RNA[keep],
      stringsAsFactors = FALSE)
    if (is.null(df$total_counts) || all(is.na(df$total_counts))) {
      df$total_counts <- Matrix::colSums(LayerData(obj, assay = "RNA", layer = "counts")[, keep, drop = FALSE])
    }
    for (gg in g) df[[gg]] <- as.numeric(cm[gg, ])
    # keep one column set across files: a gene absent from this platform is NA, not a missing column
    for (gg in genes_use) if (!gg %in% names(df)) df[[gg]] <- NA_real_
    cell_rows[[basename(f)]] <- df
    rm(obj, cm, df); invisible(gc())
    cat("  cells from", basename(f), "\n")
  }
  # a marker absent from a platform is reported as NA, not silently dropped from the table
  cl <- do.call(rbind, cell_rows)
  saveRDS(cl, cache)
}
for (gg in genes_use) if (!gg %in% names(cl)) cl[[gg]] <- NA_real_
stopifnot(nrow(cl) > 1e5, sum(cl[[TARGET]] > 0, na.rm = TRUE) > 0)

ct_levels <- sort(unique(cl$celltype))

det <- do.call(rbind, lapply(split(cl, list(cl$platform, cl$celltype), drop = TRUE), function(d) {
  pos <- d[d[[TARGET]] > 0, ]
  data.frame(
    platform = d$platform[1], celltype = d$celltype[1],
    n_cells = nrow(d), n_myoz2_pos = nrow(pos),
    detection_fraction = round(mean(d[[TARGET]] > 0), 5),
    mean_counts = round(mean(d[[TARGET]]), 4),
    mean_cpm = round(mean(1e6 * d[[TARGET]] / d$total_counts), 3),
    frac_of_all_myoz2_pos_cells = NA_real_,
    enrichment_vs_all_cells = NA_real_,
    binom_p = NA_real_,
    stringsAsFactors = FALSE)
}))

for (pl in unique(det$platform)) {
  idx <- which(det$platform == pl)
  dd <- det[idx, ]
  n_pos_total <- sum(dd$n_myoz2_pos)
  base <- sum(dd$n_cells)
  for (i in seq_along(idx)) {
    if (!n_pos_total || !dd$n_cells[i]) next
    a <- dd$n_myoz2_pos[i]
    expected <- dd$n_cells[i] / base
    det$frac_of_all_myoz2_pos_cells[idx[i]] <- round(a / n_pos_total, 4)
    det$enrichment_vs_all_cells[idx[i]] <- round((a / n_pos_total) / expected, 3)
    det$binom_p[idx[i]] <- signif(binom.test(a, n_pos_total, p = expected)$p.value, 4)
  }
}
det$enrichment_q <- NA_real_
for (pl in unique(det$platform)) {
  i <- which(det$platform == pl)
  det$enrichment_q[i] <- p.adjust(det$binom_p[i], method = "BH")
}
write.csv(det, file.path(tab_dir, "09g_moyz2_cell_detection.csv"), row.names = FALSE)

cat("\n--- MYOZ2 detection by platform x cell type ---\n")
print(det[, c("platform", "celltype", "n_cells", "n_myoz2_pos", "detection_fraction",
              "enrichment_vs_all_cells", "binom_p", "enrichment_q")], row.names = FALSE)

codet_rows <- list()
codet_one <- function(d, scope, marker) {
  if (!marker %in% names(d) || all(is.na(d[[marker]]))) return(NULL)
  # a gene can be missing from part of the platform's objects; restrict the 2x2 to cells where
  # both genes were measured so NA never enters fisher.test (and report the restricted n)
  ok <- is.finite(d[[TARGET]]) & is.finite(d[[marker]])
  d <- d[ok, , drop = FALSE]
  if (!nrow(d) || sum(d[[TARGET]] > 0) == 0) return(NULL)
  mz <- d[[TARGET]] > 0; mk <- d[[marker]] > 0
  a <- sum(mz & mk); b <- sum(mz & !mk); c <- sum(!mz & mk); dd <- sum(!mz & !mk)
  ft <- fisher.test(matrix(c(a, b, c, dd), nrow = 2, byrow = TRUE))
  or <- if (any(c(a, b, c, dd) == 0)) {
    a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- dd + 0.5
    (a2 * d2) / (b2 * c2)
  } else (a * dd) / (b * c)
  data.frame(scope = scope, marker = marker, n_cells = nrow(d),
             myoz2_pos = sum(mz), myoz2_pos_marker_pos = a,
             myoz2_neg_marker_pos = c, odds_ratio = round(or, 3),
             fisher_p = signif(ft$p.value, 4), haldane_corrected = any(c(a, b, c, dd) == 0),
             stringsAsFactors = FALSE)
}

for (pl in c("snRNA", "scRNA")) {
  d <- cl[cl$platform == pl, ]
  for (g in c(SMC, REF)) codet_rows[[paste(pl, "all", g)]] <- codet_one(d, paste0(pl, ": all cells"), g)
  # within cell types with >= 20 MYOZ2+ cells (low-power by construction; counts reported)
  for (ct in unique(d$celltype)) {
    dc <- d[d$celltype == ct, ]
    if (sum(dc[[TARGET]] > 0, na.rm = TRUE) < 20) next
    for (g in SMC) codet_rows[[paste(pl, ct, g)]] <-
      codet_one(dc, paste0(pl, ": ", ct), g)
  }
}
codet <- do.call(rbind, codet_rows)
write.csv(codet, file.path(tab_dir, "09g_moyz2_cell_codetection.csv"), row.names = FALSE)
cat("\n--- co-detection (MYOZ2+ vs marker+) ---\n")
print(codet[codet$scope %in% c("snRNA: all cells", "scRNA: all cells"), ], row.names = FALSE)

# --- depth- and composition-adjusted cell-level models ----------------------
# Raw co-detection is confounded by sequencing depth: a MYOZ2+ cell tends to have more counts,
# so every marker is more likely to be called positive. Adjust the cell-type enrichment and the
# co-detection OR for log10(total counts) and cell type, in the well-powered snRNA platform.
sn0 <- cl[cl$platform == "snRNA", ]
sn0$myoz2_pos <- as.integer(sn0[[TARGET]] > 0)
sn0$log_counts <- log10(pmax(sn0$total_counts, 1))
sn0$celltype <- relevel(factor(sn0$celltype), ref = "PT")

depth <- do.call(rbind, lapply(split(sn0, sn0$myoz2_pos), function(d) data.frame(
  platform = "snRNA", status = if (d$myoz2_pos[1] == 1) "MYOZ2+" else "MYOZ2-",
  n_cells = nrow(d), median_total_counts = median(d$total_counts),
  stringsAsFactors = FALSE)))
write.csv(depth, file.path(tab_dir, "09g_moyz2_depth_check.csv"), row.names = FALSE)

adj_rows <- list()
fit_ct <- glm(myoz2_pos ~ log_counts + celltype, data = sn0, family = binomial)
cf <- summary(fit_ct)$coefficients
row_of <- function(term, model_name, n, n_pos = NA_integer_) {
  if (!term %in% rownames(cf)) return(NULL)
  est <- cf[term, "Estimate"]; se <- cf[term, "Std. Error"]
  data.frame(model = model_name, term = term, n = n, n_myoz2_pos = n_pos,
             odds_ratio = round(exp(est), 3),
             ci_low = round(exp(est - 1.96 * se), 3),
             ci_high = round(exp(est + 1.96 * se), 3),
             p = signif(cf[term, "Pr(>|z|)"], 4), stringsAsFactors = FALSE)
}
for (ct in levels(sn0$celltype)[-1]) {
  n_pos_ct <- sum(sn0$myoz2_pos[sn0$celltype == ct])
  # a cell type with < 10 positive cells is separation-prone: keep it in the model, do not
  # report its coefficient as if it were estimable
  if (n_pos_ct < 10) next
  r <- row_of(paste0("celltype", ct), "MYOZ2 detection (ref: PT)", nrow(sn0), n_pos_ct)
  if (!is.null(r)) adj_rows[[paste("ct", ct)]] <- r
}
adj_rows[["depth"]] <- row_of("log_counts", "MYOZ2 detection (ref: PT)", nrow(sn0))

for (g in c(SMC, REF)) {
  d <- sn0[is.finite(sn0[[g]]), , drop = FALSE]
  if (!nrow(d) || length(unique(d[[g]] > 0)) < 2) next
  d$marker_pos <- as.integer(d[[g]] > 0)
  fit <- glm(marker_pos ~ myoz2_pos + log_counts + celltype, data = d, family = binomial)
  cf2 <- summary(fit)$coefficients
  if (!"myoz2_pos" %in% rownames(cf2)) next
  est <- cf2["myoz2_pos", "Estimate"]; se <- cf2["myoz2_pos", "Std. Error"]
  adj_rows[[paste("codet", g)]] <- data.frame(
    model = paste0(g, " detection (adjusted)"), term = "MYOZ2+", n = nrow(d),
    n_myoz2_pos = sum(d$myoz2_pos),
    odds_ratio = round(exp(est), 3),
    ci_low = round(exp(est - 1.96 * se), 3),
    ci_high = round(exp(est + 1.96 * se), 3),
    p = signif(cf2["myoz2_pos", "Pr(>|z|)"], 4), stringsAsFactors = FALSE)
}
adj <- do.call(rbind, adj_rows)
write.csv(adj, file.path(tab_dir, "09g_moyz2_cell_adjusted.csv"), row.names = FALSE)
cat("\n--- depth-adjusted snRNA models (log10 counts + cell type) ---\n")
print(adj, row.names = FALSE)

# --- donor-level pseudobulk (snRNA subset) ----------------------------------

sn <- cl[cl$platform == "snRNA", ]
agg_one <- function(d) {
  out <- data.frame(donor = d$donor[1], disease_group = d$disease_group[1],
                    n_cells = nrow(d), total_counts = sum(d$total_counts),
                    stringsAsFactors = FALSE)
  for (g in c(TARGET, SMC, REF)) out[[g]] <- if (g %in% names(d)) sum(d[[g]], na.rm = TRUE) else NA_real_
  out
}
don <- do.call(rbind, lapply(split(sn, sn$donor), agg_one))
cpm_cols <- c(TARGET, SMC, REF)
for (g in cpm_cols) don[[paste0(g, "_cpm")]] <- 1e6 * don[[g]] / don$total_counts

don_rows <- list()
for (scope in c("all snRNA donors", "MN snRNA donors")) {
  d <- if (scope == "all snRNA donors") don else don[don$disease_group == "MN", ]
  for (g in c(SMC, REF)) {
    if (!all(is.finite(d[[paste0(g, "_cpm")]])) || sd(d[[paste0(g, "_cpm")]]) == 0) next
    st <- spearman_one(d[[paste0(TARGET, "_cpm")]], d[[paste0(g, "_cpm")]])
    don_rows[[paste(scope, g)]] <- data.frame(
      scope = scope, gene = g, n_donors = st[["n"]], rho = round(st[["rho"]], 4),
      p = signif(st[["p"]], 4), stringsAsFactors = FALSE)
  }
}
don_cor <- do.call(rbind, don_rows)
write.csv(don_cor, file.path(tab_dir, "09g_moyz2_donor_correlation.csv"), row.names = FALSE)
cat("\n--- snRNA donor-level pseudobulk correlations ---\n"); print(don_cor, row.names = FALSE)

# --- figures ----------------------------------------------------------------

# 1. bulk correlation heatmap (MYOZ2 vs SMC markers + composites)
hm <- bulk[bulk$gene %in% SMC, ]
hm <- rbind(hm, data.frame(
  cohort = summ$cohort, cohort_label = summ$cohort_label, sample_set = summ$sample_set,
  n_samples = summ$n_samples, gene = paste0("composite: ", summ$composite),
  role = "composite", n = NA, rho = summ$rho, ci_low = NA, ci_high = NA, p = summ$p,
  q_BH_all_genes = NA, percentile_abs_rho = NA, partial_rho_adj_sample_median = summ$partial_rho,
  stringsAsFactors = FALSE))
hm$gene <- factor(hm$gene, levels = rev(c(SMC, "composite: SMC/perivascular",
                                           "composite: podocyte", "composite: fibroblast")))
hm$set_label <- factor(paste0(set_label[hm$cohort], "\n(", hm$sample_set, ")"),
                       levels = unlist(lapply(set_label[sets], function(x)
                         paste0(x, "\n(", c("all", "MN"), ")"))))
save_plot(ggplot(hm, aes(x = set_label, y = gene, fill = rho)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = ifelse(is.na(rho), "NA", sprintf("%.2f", rho))), size = 2.6) +
            scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                 midpoint = 0, limits = c(-1, 1)) +
            labs(title = "MYOZ2 vs perivascular / smooth-muscle markers in bulk kidney tissue",
                 subtitle = "Spearman rho; composites are mean z of the member genes. SMC = ACTA2, MYH11, TAGLN, PDGFRB",
                 x = NULL, y = NULL, fill = "rho") +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09g_moyz2_bulk_correlations"), width = 8.2, height = 4.6)

# 2. scatter of MYOZ2 vs SMC composite, four cohorts
sc <- do.call(rbind, lapply(sets, function(s) {
  o <- load_set(s)
  comp <- composite(o$ms, SMC)
  data.frame(cohort = set_label[[s]], myoz2 = o$ms[TARGET, ], smc = comp,
             group = ifelse(o$grp == "MN", "MN", ifelse(o$grp == "Control", "Control", "Other")),
             stringsAsFactors = FALSE)
}))
sc <- sc[is.finite(sc$myoz2) & is.finite(sc$smc), ]
sc$cohort <- factor(sc$cohort, levels = set_label[sets])
lab <- do.call(rbind, lapply(split(sc, sc$cohort), function(d) {
  data.frame(cohort = d$cohort[1], rho = spearman_one(d$myoz2, d$smc)[["rho"]])
}))
lab$label <- sprintf("rho = %.2f", lab$rho)
save_plot(ggplot(sc, aes(x = smc, y = myoz2, colour = group)) +
            geom_point(size = 0.9, alpha = 0.65) +
            geom_smooth(data = sc, aes(group = 1), method = "lm", se = FALSE,
                        colour = "grey20", linewidth = 0.5) +
            geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
                      hjust = -0.1, vjust = 1.6, size = 2.8) +
            facet_wrap(~ cohort, scales = "free") +
            scale_colour_brewer(palette = "Dark2") +
            labs(title = "Is bulk MYOZ2 a readout of perivascular tissue content?",
                 x = "SMC/perivascular composite (mean z)", y = "MYOZ2 expression",
                 colour = NULL) +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), legend.position = "bottom",
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09g_moyz2_bulk_scatter"), width = 8.4, height = 4.6)

# 3. cell-level detection
dp <- det[det$n_myoz2_pos > 0 | det$n_cells >= 20, ]
dp$celltype <- factor(dp$celltype, levels = rev(ct_levels))
save_plot(ggplot(dp, aes(x = celltype, y = 100 * detection_fraction, fill = platform)) +
            geom_col(position = position_dodge(0.75), width = 0.7) +
            geom_text(aes(label = ifelse(n_myoz2_pos > 0, n_myoz2_pos, "")),
                      position = position_dodge(0.75), hjust = -0.15, size = 2.2) +
            coord_flip() +
            scale_fill_brewer(palette = "Set2") +
            labs(title = "MYOZ2 detection across the kidney single-cell / single-nucleus atlas",
                 subtitle = "labels = positive cells; snRNA = GSE271188, scRNA = GSE171458 + GSE241302",
                 x = NULL, y = "% cells with MYOZ2 > 0", fill = NULL) +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), legend.position = "bottom",
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09g_moyz2_cell_detection"), width = 8, height = 5.2)

# --- reproducibility guards (P-28, P-68) ------------------------------------

outs <- c("09g_moyz2_bulk_correlations.csv", "09g_moyz2_bulk_summary.csv",
          "09g_moyz2_panel_reference.csv", "09g_moyz2_cell_detection.csv",
          "09g_moyz2_cell_codetection.csv", "09g_moyz2_cell_adjusted.csv",
          "09g_moyz2_depth_check.csv", "09g_moyz2_donor_correlation.csv")
for (f in outs) {
  p <- file.path(tab_dir, f)
  stopifnot(file.exists(p), file.size(p) > 0)
}
for (f in c("09g_moyz2_bulk_correlations", "09g_moyz2_bulk_scatter", "09g_moyz2_cell_detection")) {
  for (ext in c("pdf", "jpeg")) {
    p <- file.path(fig_dir, paste0(f, ".", ext))
    stopifnot(file.exists(p), file.size(p) > 0)
  }
}
# re-read and reproduce one quoted value from the artefact
chk <- read.csv(file.path(tab_dir, "09g_moyz2_bulk_correlations.csv"), stringsAsFactors = FALSE)
cat(sprintf("\nre-read check: %d rows; MYOZ2 / discovery / SMC markers: %s\n", nrow(chk),
            paste(sprintf("%s=%.3f", chk$gene[chk$cohort == "disc_glom_array" &
                                               chk$sample_set == "all" &
                                               chk$gene %in% SMC],
                          chk$rho[chk$cohort == "disc_glom_array" &
                                    chk$sample_set == "all" & chk$gene %in% SMC]),
                  collapse = ", ")))

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] tools/moyz2_provenance.R"),
             paste0("bulk rows: ", nrow(bulk), "; cell rows: ", nrow(cl)),
             paste0("MYOZ2+ cells: ", sum(cl[[TARGET]] > 0, na.rm = TRUE))),
           file.path(PROJ, "logs", "09g_moyz2_provenance.log"))

cat("\n09g (tools/moyz2_provenance.R) done\n")
