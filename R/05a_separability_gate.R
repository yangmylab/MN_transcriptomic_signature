#!/usr/bin/env Rscript
# Paper 1 - downstream gate: can MN be separated from OTHER glomerulopathies at all?
#
# Phase B showed that at |log2FC| > 1 only ~18 (FSGS), 60 (MCD) and 18 (non-MN) glomerular genes
# separate MN from its mimics, against 584 for MN vs healthy donors. Before investing in WGCNA and
# a LASSO/SVM-RFE/Boruta panel, this script measures the thing the panel is supposed to deliver:
# discrimination of MN from other proteinuric glomerulopathies, trained on discovery and tested
# blind in ERCB (0 % patient overlap).
#
# Design
#   training   : disc_glom_array, non-control samples only (MN 51, FSGS 31, MCD 19, Other 46)
#   external   : val_glom_array, non-control samples (MN 21, FSGS 18, MCD 14, Other 117)
#   independent : cross_glom_rnaseq. Patient-level linkage to the array cohorts is NOT possible
#                from the downloaded metadata (GSE197307/GSE182380 titles are local sample indices and
#                the series matrices carry only disease + tissue), so this cohort is treated as
#                independent - NOT as the same patients measured twice.
#   features   : top-K genes by limma |t| for MN vs the other classes, reselected INSIDE every
#                outer CV fold (nested), K in 50/100/200/500/1000/all
#   models     : elastic-net logistic (alpha = 0.5) for binary MN vs rest, and elastic-net
#                multinomial MN/FSGS/MCD/Other; lambda by inner 5-fold cv.glmnet
#   scaling    : genes z-scored within each cohort (cross-platform transfer)
#
# Outputs
#   output/tables/separability_nestedcv.csv
#   output/tables/separability_external.csv
#   output/tables/separability_confusion.csv
#   output/tables/separability_selected_genes.csv
#   output/models_gate/05a_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(limma); library(glmnet); library(pROC); library(ggplot2) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_gate")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20261018)

K_GRID <- c(50, 100, 200, 500, 1000, NA)
ALPHA <- 0.5
N_OUTER <- 5
N_INNER <- 5

save_plot <- function(p, stem, width = 7.5, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

# --- data assembly -----------------------------------------------------------

mg <- function(set, drop_control = TRUE) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  keep <- ann$disease != "UNKNOWN"
  if (drop_control) keep <- keep & grp != "Control"
  ms <- ms[, keep, drop = FALSE]
  grp <- grp[keep]
  lv <- if (drop_control) c("FSGS", "MCD", "Other", "MN") else c("Control", "FSGS", "MCD", "Other", "MN")
  list(ms = ms, grp = factor(grp, levels = lv), samples = ann$gsm[keep])
}

tr <- mg("disc_glom_array")
ex <- mg("val_glom_array")
xp <- mg("cross_glom_rnaseq")

cat("training (disc_glom_array):", ncol(tr$ms), "samples\n")
print(table(tr$grp))
cat("\nexternal (val_glom_array):", ncol(ex$ms), "samples\n")
print(table(ex$grp))
cat("\ncross-platform (cross_glom_rnaseq):", ncol(xp$ms), "samples\n")
print(table(xp$grp))

genes <- Reduce(intersect, list(rownames(tr$ms), rownames(ex$ms), rownames(xp$ms)))
cat("\ncommon genes:", length(genes), "\n")

# Two steps, in this order, and both use the DISCOVERY cohort only:
#
#   1. cross-platform harmonisation - each validation cohort is mapped, gene by gene, onto the
#      discovery reference distribution (rank-preserving quantile mapping). This is required because
#      the discovery cohort is HuGene-2.1-ST and ERCB is HG-U133: applying discovery mu/sd without it
#      leaves the typical gene at |z| ~ 4 instead of ~0.6 (audit 2026-11-18), the linear predictor
#      saturates and the transfer measures platform shift rather than biology.
#   2. z-scoring with the discovery mu/sd. Cohort-internal z-scoring (the pre-audit version) is
#      transductive and must not come back - see P-58.
#
# The mapping uses no labels: it transfers the marginal distribution of each gene, so it cannot
# encode case/control status.
tr_mu <- rowMeans(tr$ms[genes, , drop = FALSE])
tr_sd <- apply(tr$ms[genes, , drop = FALSE], 1, sd)
tr_sd[!is.finite(tr_sd) | tr_sd == 0] <- 1

ref_sorted <- t(apply(tr$ms[genes, , drop = FALSE], 1, sort))

# rank -> reference quantile, per gene. NA values stay NA (06b/03b leaves NA where a release does not
# measure a gene) and are zeroed by std() below, exactly as before.
quantile_map <- function(Mte, ref) {

  out <- Mte

  for (i in seq_len(nrow(Mte))) {

    v <- Mte[i, ]
    ok <- is.finite(v)

    if (sum(ok) < 3) { out[i, ] <- NA_real_; next }

    r <- rank(v[ok], ties.method = "average")
    idx <- pmin(ncol(ref), pmax(1L, round(r / (sum(ok) + 1) * ncol(ref))))

    out[i, ] <- NA_real_
    out[i, ok] <- ref[i, idx]

  }

  out

}

std <- function(M) {
  Z <- sweep(sweep(M, 1, tr_mu), 1, tr_sd, "/")
  # genes that are all-NA in a release survive the gene filter; glmnet/pROC cannot take them
  Z[!is.finite(Z)] <- 0
  Z
}

X_tr <- std(tr$ms[genes, , drop = FALSE])
X_ex <- std(quantile_map(ex$ms[genes, , drop = FALSE], ref_sorted))
X_xp <- std(quantile_map(xp$ms[genes, , drop = FALSE], ref_sorted))

# --- feature ranking (reselected inside every fold) --------------------------

rank_genes <- function(X, y, k) {
  d <- model.matrix(~ 0 + y)
  colnames(d) <- levels(y)
  fit <- eBayes(lmFit(X, d), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = which(colnames(d) == "MN"), number = Inf, sort.by = "t")
  tt <- tt[tt$logFC != 0, ]
  g <- rownames(tt)
  if (is.na(k)) g else head(g, min(k, length(g)))
}

fit_binary <- function(Xtr, ytr) {
  ybin <- factor(ifelse(ytr == "MN", "MN", "rest"), levels = c("rest", "MN"))
  cv <- cv.glmnet(t(Xtr), ybin, family = "binomial", alpha = ALPHA,
                  nfolds = min(N_INNER, min(table(ybin))), type.measure = "auc")
  list(model = cv$glmnet.fit, lambda = cv$lambda.min, levels = levels(ybin))
}

pred_binary <- function(fit, Xte) {
  as.numeric(predict(fit$model, newx = t(Xte), s = fit$lambda, type = "response"))
}

fit_multi <- function(Xtr, ytr) {
  cv <- cv.glmnet(t(Xtr), ytr, family = "multinomial", alpha = ALPHA,
                  nfolds = min(N_INNER, min(table(ytr))), type.measure = "deviance")
  list(model = cv$glmnet.fit, lambda = cv$lambda.min)
}

pred_multi <- function(fit, Xte) {
  p <- predict(fit$model, newx = t(Xte), s = fit$lambda, type = "response")
  p <- drop(p)
  data.frame(p)
}

# --- nested CV on the discovery cohort ---------------------------------------

strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (lv in levels(y)) {
    idx <- which(y == lv)
    f[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  f
}

nested_rows <- list()

# 10 repeats x 5 folds, with one fold partition drawn per repeat. A single 5-fold draw leaves the
# AUC standard deviation resting on five values; repeated CV also makes the fold assignment explicit.
N_REPEAT_GATE <- 10L
set.seed(20261018)
folds_by_rep <- lapply(seq_len(N_REPEAT_GATE), function(r) strat_folds(tr$grp, N_OUTER))

for (K in K_GRID) {

  aucs <- numeric(0)

  for (r in seq_len(N_REPEAT_GATE)) {

    f <- folds_by_rep[[r]]

    for (fold in seq_len(N_OUTER)) {

      te <- which(f == fold)
      trn <- which(f != fold)

      g <- rank_genes(X_tr[, trn, drop = FALSE], tr$grp[trn], K)

      fit <- fit_binary(X_tr[g, trn, drop = FALSE], tr$grp[trn])
      pr <- pred_binary(fit, X_tr[g, te, drop = FALSE])

      ybin_te <- factor(ifelse(tr$grp[te] == "MN", "MN", "rest"), levels = c("rest", "MN"))
      r <- roc(response = ybin_te, predictor = pr, direction = "<", quiet = TRUE)

      aucs <- c(aucs, as.numeric(auc(r)))
    }
  }

  nested_rows[[as.character(K)]] <- data.frame(
    model = "binary MN vs rest", K = ifelse(is.na(K), "all", K),
    nested_cv_auc_mean = round(mean(aucs), 4),
    nested_cv_auc_sd = round(sd(aucs), 4),
    nested_cv_auc_min = round(min(aucs), 4),
    nested_cv_auc_max = round(max(aucs), 4),
    n_evaluations = length(aucs), repeats = N_REPEAT_GATE,
    stringsAsFactors = FALSE)
}

nested <- do.call(rbind, nested_rows)

# --- final models on the whole discovery cohort, tested externally -----------

ext_rows <- list()
roc_list <- list()
conf_list <- list()
gene_rows <- list()

# best K by nested CV is reported, but every K is transferred for transparency
for (K in K_GRID) {

  g <- rank_genes(X_tr, tr$grp, K)

  gene_rows[[as.character(K)]] <- data.frame(
    K = ifelse(is.na(K), "all", K), gene = g, stringsAsFactors = FALSE)

  fit <- fit_binary(X_tr[g, , drop = FALSE], tr$grp)

  for (cohort in c("ERCB (external)", "RNA-seq (independent cohort)")) {

    Xte <- if (cohort == "ERCB (external)") X_ex[g, , drop = FALSE] else X_xp[g, , drop = FALSE]
    yte <- if (cohort == "ERCB (external)") ex$grp else xp$grp

    pr <- pred_binary(fit, Xte)
    ybin_te <- factor(ifelse(yte == "MN", "MN", "rest"), levels = c("rest", "MN"))
    r <- roc(response = ybin_te, predictor = pr, direction = "<", quiet = TRUE)

    ext_rows[[paste(K, cohort)]] <- data.frame(
      model = "binary MN vs rest", K = ifelse(is.na(K), "all", K), cohort = cohort,
      auc = round(as.numeric(auc(r)), 4),
      ci_low = round(ci.auc(r, method = "delong")[1], 4),
      ci_high = round(ci.auc(r, method = "delong")[3], 4),
      n_MN = sum(yte == "MN"), n_other = sum(yte != "MN"), stringsAsFactors = FALSE)

    roc_list[[paste(K, cohort)]] <- data.frame(
      K = ifelse(is.na(K), "all", K), cohort = cohort,
      specificity = rev(r$specificities), sensitivity = rev(r$sensitivities))
  }

  # multinomial: does the model know WHICH glomerulopathy it is looking at?
  fm <- fit_multi(X_tr[g, , drop = FALSE], tr$grp)
  pm <- pred_multi(fm, X_ex[g, , drop = FALSE])
  pred_cls <- colnames(pm)[apply(pm, 1, which.max)]

  conf_list[[as.character(K)]] <- data.frame(
    K = ifelse(is.na(K), "all", K), truth = ex$grp, predicted = pred_cls,
    stringsAsFactors = FALSE)
}

external <- do.call(rbind, ext_rows)
roc_df <- do.call(rbind, roc_list)
conf_df <- do.call(rbind, conf_list)
genes_df <- do.call(rbind, gene_rows)

write.csv(nested, file.path(tab_dir, "separability_nestedcv.csv"), row.names = FALSE)
write.csv(external, file.path(tab_dir, "separability_external.csv"), row.names = FALSE)
write.csv(conf_df, file.path(tab_dir, "separability_confusion.csv"), row.names = FALSE)
write.csv(genes_df, file.path(tab_dir, "separability_selected_genes.csv"), row.names = FALSE)

# --- secondary endpoint: MN vs everything including healthy donors -----------
# (this is the axis Phase B already showed to be strong; run it so the contrast is explicit)

tr_all <- mg("disc_glom_array", drop_control = FALSE)
ex_all <- mg("val_glom_array", drop_control = FALSE)

genes_all <- intersect(rownames(tr_all$ms), rownames(ex_all$ms))
a_mu <- rowMeans(tr_all$ms[genes_all, , drop = FALSE])
a_sd <- apply(tr_all$ms[genes_all, , drop = FALSE], 1, sd)
a_sd[!is.finite(a_sd) | a_sd == 0] <- 1
a_ref <- t(apply(tr_all$ms[genes_all, , drop = FALSE], 1, sort))
astd <- function(M) {
  Z <- sweep(sweep(M, 1, a_mu), 1, a_sd, "/")
  Z[!is.finite(Z)] <- 0
  Z
}
Xa_tr <- astd(tr_all$ms[genes_all, , drop = FALSE])
Xa_ex <- astd(quantile_map(ex_all$ms[genes_all, , drop = FALSE], a_ref))

g_all <- rank_genes(Xa_tr, tr_all$grp, 500)
fit_all <- fit_binary(Xa_tr[g_all, , drop = FALSE], tr_all$grp)
pr_all <- pred_binary(fit_all, Xa_ex[g_all, , drop = FALSE])
y_all <- factor(ifelse(ex_all$grp == "MN", "MN", "rest"), levels = c("rest", "MN"))
r_all <- roc(response = y_all, predictor = pr_all, direction = "<", quiet = TRUE)

secondary <- data.frame(
  endpoint = "MN vs all non-MN including healthy donors", K = 500,
  cohort = "ERCB (external)", auc = round(as.numeric(auc(r_all)), 4),
  n_MN = sum(ex_all$grp == "MN"), n_other = sum(ex_all$grp != "MN"),
  stringsAsFactors = FALSE)

write.csv(secondary, file.path(tab_dir, "separability_secondary.csv"), row.names = FALSE)

# --- figures -----------------------------------------------------------------

p <- ggplot(roc_df, aes(specificity, sensitivity, colour = K)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  facet_wrap(~ cohort) +
  scale_x_reverse(limits = c(1, 0)) +
  coord_equal() +
  labs(title = "MN vs other glomerulopathies: ROC of the elastic-net gate",
       subtitle = "trained on NEPTUNE glomerular array; K = genes reselected per model",
       x = "specificity", y = "sensitivity", colour = "genes (K)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "05a_gate_roc"), width = 8.5, height = 5)

kk <- external[external$cohort == "ERCB (external)", ]
kk$K <- factor(kk$K, levels = c("50", "100", "200", "500", "1000", "all"))

p <- ggplot(kk, aes(K, auc, group = 1)) +
  geom_point(size = 2.6, colour = "#C0392B") +
  geom_line(colour = "#C0392B") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, colour = "grey40") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  ylim(0.4, 1) +
  labs(title = "External discrimination in ERCB by number of features",
       subtitle = "error bars: DeLong 95 % CI; dashed line = chance",
       x = "genes selected in discovery (K)", y = "AUC (MN vs FSGS/MCD/Other)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "05a_gate_auc_vs_k"), width = 7, height = 4.6)

# K is chosen by nested CV on the discovery cohort, never by the external AUC: selecting the
# feature count that looks best in ERCB and then drawing the ERCB confusion matrix at that K is
# selection on the validation data.
best_k <- nested$K[which.max(nested$nested_cv_auc_mean)]

cf <- conf_df[conf_df$K == best_k, ]
cf$predicted <- factor(cf$predicted, levels = c("MN", "FSGS", "MCD", "Other"))
cf$truth <- factor(cf$truth, levels = c("MN", "FSGS", "MCD", "Other"))

cf_counts <- as.data.frame(table(truth = cf$truth, predicted = cf$predicted),
                          stringsAsFactors = FALSE)

p <- ggplot(cf_counts, aes(truth, predicted, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_gradient(low = "#F4F6F7", high = "#C0392B") +
  labs(title = paste0("Multiclass gate on ERCB (K = ", best_k, ")"),
       x = "true diagnosis", y = "predicted class", fill = "samples") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "05a_gate_confusion"), width = 6, height = 5)

# --- report ------------------------------------------------------------------

cat("\n================ NESTED CV (discovery) ================\n")
print(nested, row.names = FALSE)

cat("\n================ EXTERNAL ================\n")
print(external, row.names = FALSE)

cat("\n================ SECONDARY ENDPOINT ================\n")
print(secondary, row.names = FALSE)

cat("\n================ MULTICLASS CONFUSION (best K =", best_k, ") ================\n")
print(table(truth = cf$truth, predicted = cf$predicted))

cat("\naccuracy:", round(mean(cf$truth == cf$predicted), 3),
    "| balanced accuracy:",
    round(mean(sapply(levels(cf$truth), function(l) mean(cf$predicted[cf$truth == l] == l))), 3), "\n")

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 05a_separability_gate.R"),
             capture.output(print(nested, row.names = FALSE)),
             capture.output(print(external, row.names = FALSE)),
             capture.output(print(secondary, row.names = FALSE))),
           file.path(PROJ, "logs", "separability_gate.log"))

cat("\nDone.\n")
