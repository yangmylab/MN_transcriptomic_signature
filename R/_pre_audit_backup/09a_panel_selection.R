#!/usr/bin/env Rscript
# Paper 1 - Phase E step 1: three-method biomarker panel selection on the Phase D candidate space.
#
# The separability gate (05a/05b) showed MN is learnable from other glomerulopathies (external AUC
# 0.87-0.92) with a 13,379-gene elastic net. Phase D then produced a defensible shortlist:
# output/tables/phaseE_candidate_core.csv, genes supported by >= 2 independent evidence lines, every
# hub gene belonging to an ERCB-preserved module.
#
# This step shrinks it to a 4-10 gene panel with three orthogonal selectors, all run INSIDE the
# cross-validation folds so the reported AUC is not selection-optimistic:
#   LASSO    glmnet alpha = 1, lambda = lambda.1se (5-fold inner CV)
#   SVM-RFE  linear-SVM recursive feature elimination on |w|, panel size chosen by 3-fold inner CV
#   Boruta   all-relevant random-forest selection, maxRuns = 100, Confirmed features only
# Consensus: a gene enters the panel when it is selected by >= 2 methods at selection frequency
# >= FREQ_CUT, measured over every fold of every repeat, so one lucky fold cannot carry a gene in.
#
# Scaling is fit on the TRAINING fold only (mean/sd) and applied to the test fold - no leakage.
# Endpoint: MN vs other glomerulopathy (controls excluded), matching the gate.
#
# Outputs
#   output/tables/09a_panel_nestedcv.csv            per fold/repeat AUC, panel size, chosen genes
#   output/tables/09a_panel_selection_frequency.csv gene x method selection frequency
#   output/tables/09a_panel_consensus.csv           the panel with its Phase D evidence columns
#   output/tables/09a_panel_cv_summary.csv          per-method nested-CV AUC summary
#   output/tables/09a_panel_runinfo.csv             design, seed, counts, runtime
#   output/models_panel/09a_{frequency,consensus_counts,cv_auc}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({
  library(glmnet); library(Boruta); library(randomForest); library(e1071)
  library(pROC); library(ggplot2); library(future); library(furrr)
})

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_panel")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

SEED        <- 20261019L
N_REPEAT    <- 10L
N_FOLD      <- 5L
FREQ_CUT    <- 0.60
MAX_PANEL   <- 10L
MIN_PANEL   <- 4L
SVM_SIZES   <- c(4, 6, 8, 10, 15, 20, 30, 40)
BORUTA_RUNS <- 100L
N_WORKERS   <- max(1L, min(10L, parallel::detectCores() - 2L))

save_plot <- function(p, stem, width = 8, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

# --- data ------------------------------------------------------------------
ms  <- readRDS(file.path(set_dir, "disc_glom_array_expr_symbol.rds"))
ann <- read.csv(file.path(set_dir, "disc_glom_array_samples.csv"), stringsAsFactors = FALSE)
stopifnot(identical(colnames(ms), ann$gsm))

grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
keep <- grp != "Control"

Xraw <- t(ms[, keep, drop = FALSE])
y_all <- factor(ifelse(grp[keep] == "MN", "MN", "rest"), levels = c("rest", "MN"))

cand  <- read.csv(file.path(tab_dir, "phaseE_candidate_core.csv"), stringsAsFactors = FALSE)
cand  <- cand[cand$compartment == "glom", ]

# A panel gene must be measurable in the cohort used to validate it: the ERCB glomerular release is a
# two-platform union, so a gene can be present in discovery and absent for part of ERCB (GPC6 is NA
# in 125/196 samples). Availability is therefore a SELECTION CONSTRAINT, not a post-hoc drop.
val_ms <- readRDS(file.path(set_dir, "val_glom_array_expr_symbol.rds"))
COVERAGE_CUT <- 0.95

genes <- intersect(cand$gene, colnames(Xraw))
cov_val <- vapply(genes, function(g) {
  if (!g %in% rownames(val_ms)) return(0)
  mean(!is.na(val_ms[g, ]))
}, numeric(1))
dropped <- genes[cov_val < COVERAGE_CUT]
genes <- genes[cov_val >= COVERAGE_CUT]
cat("candidate genes:", length(genes) + length(dropped), "| measurable in >=",
    round(100 * COVERAGE_CUT), "% of ERCB:", length(genes), "\n")
if (length(dropped)) cat("  dropped for ERCB coverage:", paste(dropped, collapse = ", "), "\n")
stopifnot(length(genes) >= 50)
Xraw <- Xraw[, genes, drop = FALSE]
Xraw[!is.finite(Xraw)] <- 0

cat("discovery:", nrow(Xraw), "samples x", ncol(Xraw), "candidate genes\n")
print(table(y_all))

# --- selectors (training data only) ---------------------------------------
select_lasso <- function(Xtr, ytr) {
  cv <- try(cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1,
                      nfolds = min(5L, min(table(ytr))), type.measure = "auc"), silent = TRUE)
  if (inherits(cv, "try-error")) return(character(0))
  b <- as.matrix(coef(cv$glmnet.fit, s = cv$lambda.1se))
  b <- b[rownames(b) != "(Intercept)", , drop = FALSE]
  rownames(b)[b[, 1] != 0]
}

# linear-SVM RFE: drop the weakest |w| chunk repeatedly, keep the elimination ranking
svm_rfe_rank <- function(Xtr, ytr, step_frac = 0.10) {
  keep <- colnames(Xtr)
  elim <- character(0)
  while (length(keep) > 3) {
    fit <- try(e1071::svm(x = Xtr[, keep, drop = FALSE], y = ytr, kernel = "linear",
                          cost = 1, scale = FALSE), silent = TRUE)
    if (inherits(fit, "try-error")) break
    w <- as.numeric(t(fit$coefs) %*% fit$SV)
    names(w) <- keep
    n_drop <- max(1L, floor(length(keep) * step_frac))
    drop <- names(sort(abs(w)))[seq_len(n_drop)]
    elim <- c(elim, drop)                 # weakest eliminated first
    keep <- setdiff(keep, drop)
  }
  if (length(keep) > 1) {
    fit <- try(e1071::svm(x = Xtr[, keep, drop = FALSE], y = ytr, kernel = "linear",
                          cost = 1, scale = FALSE), silent = TRUE)
    w <- if (inherits(fit, "try-error")) rep(0, length(keep)) else
      as.numeric(t(fit$coefs) %*% fit$SV)
    names(w) <- keep
    keep <- names(sort(abs(w)))           # weakest first
  }
  c(elim, keep)                           # ordered weakest -> strongest
}

select_svm <- function(Xtr, ytr) {
  rank_order <- svm_rfe_rank(Xtr, ytr)
  sizes <- SVM_SIZES[SVM_SIZES <= length(rank_order)]
  if (!length(sizes)) return(character(0))
  # choose the panel size by 3-fold inner CV on the training fold only
  set.seed(SEED)
  inner <- caret_folds(ytr, 3L)
  score <- vapply(sizes, function(k) {
    sel <- tail(rank_order, k)
    aucs <- vapply(inner, function(te) {
      tr <- setdiff(seq_len(nrow(Xtr)), te)
      fit <- try(glm(ytr[tr] ~ ., family = binomial,
                     data = data.frame(Xtr[tr, sel, drop = FALSE], ytr = ytr[tr])), silent = TRUE)
      if (inherits(fit, "try-error")) return(NA_real_)
      pr <- try(predict(fit, newdata = data.frame(Xtr[te, sel, drop = FALSE]), type = "response"),
                silent = TRUE)
      if (inherits(pr, "try-error") || length(unique(pr)) < 2) return(NA_real_)
      as.numeric(pROC::auc(pROC::roc(ytr[te], pr, levels = c("rest", "MN"),
                                     direction = "<", quiet = TRUE)))
    }, numeric(1))
    mean(aucs, na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(score))) return(tail(rank_order, min(sizes)))
  tail(rank_order, sizes[which.max(score)])
}

select_boruta <- function(Xtr, ytr) {
  b <- try(Boruta(x = as.data.frame(Xtr), y = ytr, maxRuns = BORUTA_RUNS,
                  doTrace = 0, getImp = getImpRfZ), silent = TRUE)
  if (inherits(b, "try-error")) return(character(0))
  names(b$finalDecision)[as.character(b$finalDecision) == "Confirmed"]
}

# stratified folds without caret
caret_folds <- function(y, k) {
  idx <- seq_along(y)
  out <- vector("list", k)
  for (lv in levels(y)) {
    ids <- which(y == lv)
    ids <- sample(ids)
    parts <- split(ids, rep(seq_len(k), length.out = length(ids)))
    for (i in seq_len(k)) out[[i]] <- c(out[[i]], parts[[i]])
  }
  lapply(out, sort)
}

# --- nested CV -------------------------------------------------------------
one_fold <- function(r, k) {
  set.seed(SEED + 1000L * r + k)
  folds <- caret_folds(y_all, N_FOLD)
  te <- folds[[k]]
  tr <- setdiff(seq_len(nrow(Xraw)), te)

  mu  <- colMeans(Xraw[tr, , drop = FALSE])
  sdv <- apply(Xraw[tr, , drop = FALSE], 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  Xtr <- sweep(sweep(Xraw[tr, , drop = FALSE], 2, mu), 2, sdv, "/")
  Xte <- sweep(sweep(Xraw[te, , drop = FALSE], 2, mu), 2, sdv, "/")
  ytr <- y_all[tr]
  yte <- y_all[te]

  out <- vector("list", 3)
  for (i in seq_along(c("lasso", "svm", "boruta"))) {
    m <- c("lasso", "svm", "boruta")[i]
    sel <- switch(m,
                  lasso  = select_lasso(Xtr, ytr),
                  svm    = select_svm(Xtr, ytr),
                  boruta = select_boruta(Xtr, ytr))
    auc <- NA_real_
    if (length(sel) >= 1) {
      fit <- try(glm(ytr ~ ., family = binomial,
                     data = data.frame(Xtr[, sel, drop = FALSE], ytr = ytr)), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        pr <- try(predict(fit, newdata = data.frame(Xte[, sel, drop = FALSE]), type = "response"),
                  silent = TRUE)
        if (!inherits(pr, "try-error") && length(unique(pr)) > 1) {
          auc <- as.numeric(pROC::auc(pROC::roc(yte, pr, levels = c("rest", "MN"),
                                                direction = "<", quiet = TRUE)))
        }
      }
    }
    out[[i]] <- data.frame(repeat_id = r, fold = k, method = m, n_selected = length(sel),
                           auc = round(auc, 4), genes = paste(sort(sel), collapse = ";"),
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

cat("\n=== nested CV:", N_REPEAT, "repeats x", N_FOLD, "folds x 3 selectors ===\n")
t0 <- Sys.time()

use_multicore <- .Platform$OS.type == "unix"
if (use_multicore) plan(multicore, workers = N_WORKERS) else plan(sequential)

tasks <- expand.grid(r = seq_len(N_REPEAT), k = seq_len(N_FOLD))
res <- future_pmap_dfr(list(tasks$r, tasks$k), function(r, k) one_fold(r, k),
                       .options = furrr_options(seed = TRUE))
plan(sequential)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
cat("nested CV finished in", elapsed, "min\n")
write.csv(res, file.path(tab_dir, "09a_panel_nestedcv.csv"), row.names = FALSE)

# --- selection frequency ---------------------------------------------------
freq_long <- do.call(rbind, lapply(split(res, res$method), function(d) {
  g <- unlist(strsplit(d$genes, ";"))
  g <- g[nzchar(g)]
  if (!length(g)) return(NULL)
  tab <- table(g)
  data.frame(method = d$method[1], gene = names(tab), n_times = as.integer(tab),
             frequency = round(as.integer(tab) / nrow(d), 4), stringsAsFactors = FALSE)
}))
rownames(freq_long) <- NULL
freq_long <- freq_long[order(freq_long$method, -freq_long$frequency), ]
write.csv(freq_long, file.path(tab_dir, "09a_panel_selection_frequency.csv"), row.names = FALSE)

wide <- reshape(freq_long[, c("gene", "method", "frequency")], idvar = "gene",
                timevar = "method", direction = "wide")
names(wide) <- sub("^frequency\\.", "", names(wide))
for (m in c("lasso", "svm", "boruta")) if (!m %in% names(wide)) wide[[m]] <- 0
wide[is.na(wide)] <- 0
wide$n_methods_freq <- (wide$lasso >= FREQ_CUT) + (wide$svm >= FREQ_CUT) + (wide$boruta >= FREQ_CUT)
wide$mean_frequency <- round(rowMeans(wide[, c("lasso", "svm", "boruta")]), 4)

sel <- wide[wide$n_methods_freq >= 2, ]
sel <- sel[order(-sel$n_methods_freq, -sel$mean_frequency), ]
if (nrow(sel) > MAX_PANEL) sel <- head(sel, MAX_PANEL)
if (nrow(sel) < MIN_PANEL) {
  cat("consensus below", MIN_PANEL, "genes; relaxing to >= 1 method at frequency >= 0.8\n")
  alt <- wide[(wide$lasso >= 0.8) + (wide$svm >= 0.8) + (wide$boruta >= 0.8) >= 1, ]
  alt <- alt[order(-alt$n_methods_freq, -alt$mean_frequency), ]
  sel <- head(alt, MAX_PANEL)
}

core_anno <- cand[, c("gene", "is_hub", "hub_module", "hub_kME", "module_Z", "is_metaDEG",
                      "meta_contrast", "meta_logFC", "meta_padj", "gate_min_K",
                      "top_celltype", "celltype_tau", "n_evidence")]
sel <- merge(sel, core_anno, by = "gene", all.x = TRUE, sort = FALSE)
sel <- sel[order(-sel$n_methods_freq, -sel$mean_frequency), ]
write.csv(sel, file.path(tab_dir, "09a_panel_consensus.csv"), row.names = FALSE)

cat("\npanel (", nrow(sel), "genes ):\n")
print(sel[, c("gene", "lasso", "svm", "boruta", "n_methods_freq", "mean_frequency",
              "hub_module", "meta_contrast", "n_evidence")], row.names = FALSE)

cv_sum <- do.call(rbind, lapply(split(res, res$method), function(d) {
  data.frame(method = d$method[1], n_folds = nrow(d),
             mean_auc = round(mean(d$auc, na.rm = TRUE), 4),
             sd_auc = round(sd(d$auc, na.rm = TRUE), 4),
             n_failed = sum(is.na(d$auc)),
             median_panel = median(d$n_selected), stringsAsFactors = FALSE)
}))
write.csv(cv_sum, file.path(tab_dir, "09a_panel_cv_summary.csv"), row.names = FALSE)

runinfo <- data.frame(
  item = c("seed", "repeats", "folds", "candidate_genes", "candidate_measured",
           "genes_dropped_ercb_coverage", "endpoint", "samples", "n_MN", "n_nonMN",
           "freq_cut", "max_panel", "boruta_runs", "panel_size", "elapsed_min",
           "parallel_workers"),
  value = as.character(c(SEED, N_REPEAT, N_FOLD, nrow(cand), length(genes),
                         if (length(dropped)) paste(dropped, collapse = ";") else "none",
                         "MN vs other GN (controls excluded)", nrow(Xraw),
                         sum(y_all == "MN"), sum(y_all == "rest"), FREQ_CUT, MAX_PANEL,
                         BORUTA_RUNS, nrow(sel), elapsed, N_WORKERS)),
  stringsAsFactors = FALSE)
write.csv(runinfo, file.path(tab_dir, "09a_panel_runinfo.csv"), row.names = FALSE)

cat("\nper-method nested CV:\n"); print(cv_sum, row.names = FALSE)

# --- figures ---------------------------------------------------------------
top_genes <- unique(head(freq_long$gene[order(-freq_long$frequency)], 25))
pd <- freq_long[freq_long$gene %in% top_genes, ]
pd$method <- factor(pd$method, levels = c("lasso", "svm", "boruta"))
ord <- names(sort(tapply(pd$frequency, pd$gene, mean), decreasing = TRUE))
pd$gene <- factor(pd$gene, levels = rev(ord))

save_plot(ggplot(pd, aes(x = method, y = gene, fill = frequency)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = ifelse(frequency > 0, sprintf("%.2f", frequency), "")),
                      size = 2.4, colour = "grey15") +
            scale_fill_gradient(low = "grey95", high = "firebrick") +
            labs(title = "Selection frequency across nested CV folds",
                 subtitle = paste0(N_REPEAT, " repeats x ", N_FOLD,
                                   " folds; top 25 genes; red = always selected"),
                 x = NULL, y = NULL, fill = "frequency") +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09a_frequency"), width = 6.5, height = 7)

share <- data.frame(
  pattern = factor(c("1 method", "2 methods", "3 methods"),
                   levels = c("1 method", "2 methods", "3 methods")),
  n = c(sum(wide$n_methods_freq == 1), sum(wide$n_methods_freq == 2),
        sum(wide$n_methods_freq == 3)))

save_plot(ggplot(share, aes(x = pattern, y = n, fill = pattern)) +
            geom_col(width = 0.65, show.legend = FALSE) +
            geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
            scale_fill_brewer(palette = "YlOrRd") +
            labs(title = "Consensus across the three selectors",
                 subtitle = paste0("genes at selection frequency >= ", FREQ_CUT,
                                   "; panel size ", nrow(sel)),
                 x = NULL, y = "genes") +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09a_consensus_counts"), width = 6, height = 4.5)

cvp <- res[!is.na(res$auc), ]
save_plot(ggplot(cvp, aes(x = method, y = auc, fill = method)) +
            geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.75) +
            geom_jitter(width = 0.12, size = 0.8, alpha = 0.6) +
            geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
            coord_cartesian(ylim = c(0.4, 1)) +
            labs(title = "Nested-CV AUC by selector",
                 subtitle = "features reselected inside every fold; controls excluded",
                 x = NULL, y = "AUC") +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = "none",
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09a_cv_auc"), width = 6, height = 4.5)

cat("\n09a done\n")
