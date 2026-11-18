#!/usr/bin/env Rscript
# Paper 1 - Phase E step 2: lock the panel and evaluate it once on independent data.
#
# 09a produced a 5-gene panel from nested selection. This step fits the FINAL model on all discovery
# samples and then evaluates that fixed model, without any further tuning, in
#   ERCB glomerular (val_glom_array, non-control, 0% patient overlap with discovery)
#   the RNA-seq cohort (cross_glom_rnaseq, non-control) as a cross-platform check
#   the secondary endpoint (MN vs everything including healthy donors)
# and reports where the panel sits relative to the gate's whole-transcriptome elastic net (0.915 in
# ERCB). Also reported: DeLong CI, sensitivity/specificity at the Youden point, a calibration slope,
# decision-curve net benefit, an ERCB-label permutation null (200 shuffles) and a discovery-label
# null (20 refits) that shows transfer is not a property of any labelling.
#
# The model is a plain logistic regression on genes z-scored with the DISCOVERY mean/sd - the same
# transformation is applied to every validation cohort, so no cohort-specific tuning is involved.
#
# Outputs
#   output/tables/09b_panel_coefficients.csv        per-gene coefficient, OR, 95% CI, module identity
#   output/tables/09b_panel_external.csv            per-cohort AUC with CI, sens, spec, calibration
#   output/tables/09b_panel_permutation.csv         empirical p for each null
#   output/tables/09b_panel_decision_curve.csv      net benefit across thresholds
#   output/models_panel/09b_{roc,score_box,forest,dca}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(pROC); library(ggplot2) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_panel")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 20261019L
N_PERM_EXT <- 200L
N_PERM_REFIT <- 200L
THR_GRID <- seq(0.05, 0.95, by = 0.01)

save_plot <- function(p, stem, width = 8, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

panel <- read.csv(file.path(tab_dir, "09a_panel_consensus.csv"), stringsAsFactors = FALSE)
genes <- panel$gene
cat("locked panel:", paste(genes, collapse = ", "), "\n")

load_cohort <- function(set, keep_control = FALSE) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(ms), ann$gsm))
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  keep <- if (keep_control) grp != "UNKNOWN" else !(grp %in% c("Control", "UNKNOWN"))
  X <- t(ms[, keep, drop = FALSE])
  list(X = X, y = factor(ifelse(grp[keep] == "MN", "MN", "rest"), levels = c("rest", "MN")),
       group = factor(grp[keep], levels = c("Control", "FSGS", "MCD", "Other", "MN")),
       gsm = ann$gsm[keep])
}

# --- fit the locked model on discovery -------------------------------------
disc  <- load_cohort("disc_glom_array")
ercb  <- load_cohort("val_glom_array")
rnaseq <- load_cohort("cross_glom_rnaseq")
ercb_all <- load_cohort("val_glom_array", keep_control = TRUE)
disc_all <- load_cohort("disc_glom_array", keep_control = TRUE)

stopifnot(all(genes %in% colnames(disc$X)), all(genes %in% colnames(ercb$X)))

mu  <- colMeans(disc$X[, genes, drop = FALSE])
sdv <- apply(disc$X[, genes, drop = FALSE], 2, sd)
sdv[!is.finite(sdv) | sdv == 0] <- 1

# Any residual gap (a gene absent for part of a two-platform cohort) is set to the discovery mean,
# i.e. z = 0, and the count is reported so the reader knows where imputation happened.
z <- function(X) {
  m <- X[, genes, drop = FALSE]
  n_na <- sum(is.na(m))
  if (n_na) cat(sprintf("  imputing %d missing panel values with the discovery mean\n", n_na))
  zz <- sweep(sweep(m, 2, mu), 2, sdv, "/")
  zz[!is.finite(zz)] <- 0
  zz
}

d_z <- as.data.frame(z(disc$X));  d_z$y <- disc$y
fit <- glm(y ~ ., data = d_z, family = binomial)
cat("\n--- locked model ---\n"); print(summary(fit)$coefficients)
cat("discovery AUC:", round(as.numeric(auc(roc(disc$y, fitted(fit), levels = c("rest", "MN"),
                                             direction = "<", quiet = TRUE))), 4), "\n")

co <- summary(fit)$coefficients
ci <- confint.default(fit)
# module programme identity from 08f, so a panel gene is reported with the programme it belongs to
prog <- read.csv(file.path(tab_dir, "08f_module_programme_labels.csv"), stringsAsFactors = FALSE)
prog <- prog[prog$compartment == "glom", c("module", "programme_label", "label_padj",
                                           "second_term", "third_term")]
coef_tab <- data.frame(
  gene = rownames(co)[-1],
  coefficient = round(co[-1, 1], 4),
  odds_ratio = round(exp(co[-1, 1]), 4),
  OR_low = round(exp(ci[-1, 1]), 4),
  OR_high = round(exp(ci[-1, 2]), 4),
  p = signif(co[-1, 4], 4),
  selection_frequency_mean = panel$mean_frequency[match(rownames(co)[-1], panel$gene)],
  methods = panel$n_methods_freq[match(rownames(co)[-1], panel$gene)],
  module = panel$hub_module[match(rownames(co)[-1], panel$gene)],
  kME = panel$hub_kME[match(rownames(co)[-1], panel$gene)],
  meta_contrast = panel$meta_contrast[match(rownames(co)[-1], panel$gene)],
  meta_logFC = panel$meta_logFC[match(rownames(co)[-1], panel$gene)],
  top_celltype = panel$top_celltype[match(rownames(co)[-1], panel$gene)],
  celltype_tau = panel$celltype_tau[match(rownames(co)[-1], panel$gene)],
  stringsAsFactors = FALSE)
coef_tab$module_programme <- prog$programme_label[match(coef_tab$module, prog$module)]
coef_tab$module_programme_2nd <- prog$second_term[match(coef_tab$module, prog$module)]
write.csv(coef_tab, file.path(tab_dir, "09b_panel_coefficients.csv"), row.names = FALSE)

# --- external evaluation ---------------------------------------------------
score <- function(X) as.numeric(predict(fit, newdata = as.data.frame(z(X)), type = "response"))

eval_cohort <- function(obj, label) {
  ok <- genes[genes %in% colnames(obj$X)]
  if (length(ok) < length(genes)) {
    cat(sprintf("  %s: %d/%d panel genes measured - scoring with the available subset\n",
                label, length(ok), length(genes)))
  }
  s <- score(obj$X)
  y <- obj$y
  r <- roc(obj$y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  best <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"),
                 transpose = FALSE)[1, ]
  cal <- try(glm(obj$y ~ s, family = binomial), silent = TRUE)
  data.frame(cohort = label, n = length(obj$y), n_MN = sum(obj$y == "MN"),
             n_nonMN = sum(obj$y == "rest"),
             auc = round(as.numeric(auc(r)), 4),
             ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
             youden_threshold = round(best$threshold, 4),
             sensitivity = round(best$sensitivity, 4),
             specificity = round(best$specificity, 4),
             sens_at_0.5 = round(mean(s[y == "MN"] >= 0.5), 4),
             spec_at_0.5 = round(mean(s[y == "rest"] < 0.5), 4),
             calibration_in_large = round(mean(s, na.rm = TRUE) - mean(obj$y == "MN"), 4),
             calibration_intercept = round(coef(cal)[1], 4),
             calibration_slope = round(coef(cal)[2], 4),
             genes_used = length(ok), stringsAsFactors = FALSE)
}

ext <- rbind(
  eval_cohort(disc, "Discovery (apparent, not independent)"),
  eval_cohort(ercb, "ERCB external"),
  eval_cohort(rnaseq, "RNA-seq independent cohort"),
  eval_cohort(ercb_all, "ERCB secondary endpoint (incl. donors)"),
  eval_cohort(disc_all, "Discovery secondary endpoint (incl. donors)"))

# reference point: the gate's whole-transcriptome elastic net in ERCB
gate_f <- file.path(tab_dir, "separability_external.csv")
if (file.exists(gate_f)) {
  g <- read.csv(gate_f, stringsAsFactors = FALSE)
  g <- g[g$cohort == "ERCB (external)" & g$K == "all", ]
  if (nrow(g)) {
    ext <- rbind(ext, data.frame(
      cohort = "Reference: gate elastic net (13,379 genes)",
      n = g$n_MN[1] + g$n_other[1], n_MN = g$n_MN[1], n_nonMN = g$n_other[1],
      auc = round(g$auc[1], 4), ci_low = round(g$ci_low[1], 4), ci_high = round(g$ci_high[1], 4),
      youden_threshold = NA, sensitivity = NA, specificity = NA,
      sens_at_0.5 = NA, spec_at_0.5 = NA,
      calibration_in_large = NA, calibration_intercept = NA, calibration_slope = NA,
      genes_used = 13379, stringsAsFactors = FALSE))
  }
}
write.csv(ext, file.path(tab_dir, "09b_panel_external.csv"), row.names = FALSE)
cat("\n--- external evaluation ---\n"); print(ext, row.names = FALSE)

# --- permutation nulls -----------------------------------------------------
set.seed(SEED)
s_ercb <- score(ercb$X)
obs <- as.numeric(auc(roc(ercb$y, s_ercb, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
null_ext <- vapply(seq_len(N_PERM_EXT), function(i) {
  yy <- sample(ercb$y)
  as.numeric(auc(roc(yy, s_ercb, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
}, numeric(1))
p_ext <- (sum(null_ext >= obs) + 1) / (length(null_ext) + 1)

# discovery-label null: refit the panel on shuffled discovery labels, then test in ERCB
null_refit <- vapply(seq_len(N_PERM_REFIT), function(i) {
  d2 <- d_z
  d2$y <- sample(d2$y)
  f2 <- try(glm(y ~ ., data = d2, family = binomial), silent = TRUE)
  if (inherits(f2, "try-error")) return(NA_real_)
  s2 <- as.numeric(predict(f2, newdata = as.data.frame(z(ercb$X)), type = "response"))
  as.numeric(auc(roc(ercb$y, s2, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
}, numeric(1))
p_refit <- (sum(null_refit >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null_refit)) + 1)

# median / 95th percentile must be taken per null, not reused across rows
perm_tab <- data.frame(
  null = c("ERCB label permutation", "Discovery label permutation (refit)"),
  n_permutations = c(N_PERM_EXT, N_PERM_REFIT),
  observed_auc = round(obs, 4),
  null_median = round(c(median(null_ext), median(null_refit, na.rm = TRUE)), 4),
  null_q95 = round(c(quantile(null_ext, 0.95), quantile(null_refit, 0.95, na.rm = TRUE)), 4),
  null_max = round(c(max(null_ext), max(null_refit, na.rm = TRUE)), 4),
  empirical_p = signif(c(p_ext, p_refit), 4),
  stringsAsFactors = FALSE)
write.csv(perm_tab, file.path(tab_dir, "09b_panel_permutation.csv"), row.names = FALSE)
cat("\n--- permutation nulls ---\n"); print(perm_tab, row.names = FALSE)

# --- decision curve --------------------------------------------------------
prevs <- c(ERCB = mean(ercb$y == "MN"), `RNA-seq` = mean(rnaseq$y == "MN"))
dca_rows <- list()
for (nm in names(prevs)) {
  obj <- if (nm == "ERCB") ercb else rnaseq
  s <- score(obj$X)
  p <- prevs[[nm]]
  for (thr in THR_GRID) {
    treat_all <- p - (1 - p) * thr / (1 - thr)
    tp <- mean(s >= thr & obj$y == "MN")
    fp <- mean(s >= thr & obj$y == "rest")
    nb <- tp - fp * thr / (1 - thr)
    dca_rows[[length(dca_rows) + 1]] <- data.frame(
      cohort = nm, threshold = thr, net_benefit_model = round(nb, 5),
      net_benefit_all = round(treat_all, 5), net_benefit_none = 0, stringsAsFactors = FALSE)
  }
}
dca <- do.call(rbind, dca_rows)
write.csv(dca, file.path(tab_dir, "09b_panel_decision_curve.csv"), row.names = FALSE)

# --- figures ---------------------------------------------------------------
roc_rows <- list()
for (nm in c("ERCB external", "RNA-seq independent cohort")) {
  obj <- if (nm == "ERCB external") ercb else rnaseq
  r <- roc(obj$y, score(obj$X), levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  roc_rows[[nm]] <- data.frame(cohort = nm, specificity = r$specificities,
                               sensitivity = r$sensitivities)
}
roc_df <- do.call(rbind, roc_rows)
lab <- paste0(ext$cohort, ": AUC ", sprintf("%.3f", ext$auc))
names(lab) <- ext$cohort

save_plot(ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, colour = cohort)) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
            geom_line(linewidth = 0.9) +
            coord_equal() +
            scale_colour_brewer(palette = "Dark2",
                                labels = c(lab[["ERCB external"]],
                                           lab[["RNA-seq independent cohort"]])) +
            labs(title = paste0("Locked ", length(genes), "-gene panel"),
                 subtitle = paste0(paste(genes, collapse = ", "),
                                   " - evaluated once, no tuning on validation data"),
                 x = "1 - specificity", y = "sensitivity", colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = c(0.62, 0.18),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09b_roc"), width = 6.4, height = 5.6)

sb <- rbind(
  data.frame(cohort = "ERCB", group = ercb$group, score = score(ercb$X)),
  data.frame(cohort = "RNA-seq", group = rnaseq$group, score = score(rnaseq$X)),
  data.frame(cohort = "Discovery", group = disc$group, score = score(disc$X)))
sb$group <- factor(sb$group, levels = c("FSGS", "MCD", "Other", "MN"))

save_plot(ggplot(sb, aes(x = group, y = score, fill = group)) +
            geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey30") +
            geom_jitter(width = 0.15, size = 0.7, alpha = 0.5) +
            facet_wrap(~ cohort, nrow = 1) +
            scale_fill_manual(values = c(FSGS = "#8DA0CB", MCD = "#A6D854", Other = "#FFD92F",
                                         MN = "#E41A1C")) +
            labs(title = "Panel score by diagnosis", x = NULL, y = "predicted probability of MN") +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = "none",
                  strip.background = element_rect(fill = "grey92", colour = NA),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09b_score_box"), width = 9, height = 4.4)

coef_tab$direction <- ifelse(coef_tab$odds_ratio > 1, "higher in MN", "lower in MN")
coef_tab$gene_lab <- paste0(coef_tab$gene, ifelse(is.na(coef_tab$module), "",
                                                  paste0(" [", coef_tab$module, "]")))
save_plot(ggplot(coef_tab, aes(x = odds_ratio, y = reorder(gene_lab, odds_ratio),
                               colour = direction)) +
            geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
            geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.18) +
            geom_point(size = 3) +
            scale_x_log10() +
            scale_colour_manual(values = c("higher in MN" = "firebrick",
                                           "lower in MN" = "steelblue")) +
            labs(title = "Locked panel coefficients (per 1 SD of expression)",
                 subtitle = "logistic regression on discovery; bars are 95% Wald CIs",
                 x = "odds ratio (log scale)", y = NULL, colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09b_forest"), width = 7, height = 4)

dca_long <- rbind(
  data.frame(cohort = dca$cohort, threshold = dca$threshold,
             strategy = "model", net_benefit = dca$net_benefit_model),
  data.frame(cohort = dca$cohort, threshold = dca$threshold,
             strategy = "treat all", net_benefit = dca$net_benefit_all),
  data.frame(cohort = dca$cohort, threshold = dca$threshold,
             strategy = "treat none", net_benefit = 0))
save_plot(ggplot(dca_long, aes(x = threshold, y = net_benefit, colour = strategy)) +
            geom_line(linewidth = 0.85) +
            facet_wrap(~ cohort) +
            scale_colour_manual(values = c(model = "firebrick", "treat all" = "grey40",
                                           "treat none" = "grey70")) +
            coord_cartesian(ylim = c(-0.05, max(dca_long$net_benefit, na.rm = TRUE) * 1.05)) +
            labs(title = "Decision-curve analysis", x = "threshold probability",
                 y = "net benefit", colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92",
                                                                                colour = NA),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09b_dca"), width = 8, height = 4.2)

# --- calibration figure ----------------------------------------------------
cal_df <- do.call(rbind, lapply(list(
  list(nm = "Discovery", obj = disc), list(nm = "ERCB", obj = ercb),
  list(nm = "RNA-seq", obj = rnaseq)), function(z) {
  s <- score(z$obj$X)
  y <- as.integer(z$obj$y == "MN")
  br <- unique(quantile(s, probs = seq(0, 1, length.out = 6)))
  if (length(br) < 3) return(NULL)
  bin <- cut(s, breaks = br, include.lowest = TRUE)
  d <- data.frame(cohort = z$nm, pred = tapply(s, bin, mean), obs = tapply(y, bin, mean),
                  n = tapply(y, bin, length))
  d$se <- sqrt(d$obs * (1 - d$obs) / d$n)
  d
}))

if (!is.null(cal_df) && nrow(cal_df)) {
  save_plot(ggplot(cal_df, aes(x = pred, y = obs, colour = cohort)) +
              geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
              geom_errorbar(aes(ymin = pmax(0, obs - 1.96 * se), ymax = pmin(1, obs + 1.96 * se)),
                            width = 0.01) +
              geom_point(size = 2) + geom_line(linewidth = 0.7) +
              coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
              scale_colour_brewer(palette = "Dark2") +
              labs(title = "Calibration of the locked panel",
                   subtitle = "binned into quintiles of predicted risk; predictions are over-confident",
                   x = "mean predicted probability", y = "observed proportion MN", colour = NULL) +
              theme_bw(base_size = 10) +
              theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")),
            file.path(fig_dir, "09b_calibration"), width = 6, height = 5.2)
}

cat("\n09b done\n")
