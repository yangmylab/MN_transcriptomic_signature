#!/usr/bin/env Rscript
# Track B - the panel-size frontier: is a small panel a choice or a ceiling? (2026-11-18)
#
# The manuscript locks four genes. A reviewer's first question is "why four when your whole-transcriptome
# model scores higher?". This script answers it quantitatively: for K = 2..50, rank genes by moderated t
# INSIDE each training fold (no leakage), fit a plain logistic model, and report both the internal CV
# AUC and the external ERCB AUC of the locked-at-full-discovery model. The locked four-gene panel is
# plotted as a reference point, so the frontier shows what parsimony costs.
#
# Inputs   work/bench/10b_harmonised.rds, output/tables/09a_panel_consensus.csv
# Outputs  output/tables/10d_panel_size_frontier.csv
#          output/models_trackb/10d_panel_size_frontier.{pdf,jpeg}
# parallel: serial - 50 folds x one limma fit each; fork overhead would dominate

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(limma); library(pROC); library(ggplot2) })

tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_trackb")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20261118)

B <- readRDS(file.path(PROJ, "work", "bench", "10b_harmonised.rds"))
locked <- read.csv(file.path(tab_dir, "09a_panel_consensus.csv"), stringsAsFactors = FALSE)$gene

K_GRID <- c(2L, 3L, 4L, 6L, 8L, 12L, 20L, 50L)
N_REPEAT <- 10L
N_FOLD <- 5L

y <- B$y$disc
Xi <- B$z_harmonised$disc
Xe <- B$z_harmonised$ercb
Xr <- B$z_harmonised$rna

# one sample per patient (within a compartment these are already unique; asserted below)
keep <- !duplicated(B$patient$disc)
stopifnot(identical(keep, rep(TRUE, length(keep))))   # documented: 147 samples = 147 patients

auc_of <- function(yy, s) {

  if (length(unique(yy)) < 2) return(NA_real_)
  as.numeric(auc(roc(yy, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))

}

# rank genes by moderated t on a training block only
rank_on <- function(Xtr, ytr) {

  d <- model.matrix(~ 0 + factor(ytr))
  colnames(d) <- levels(factor(ytr))
  fit <- eBayes(lmFit(Xtr, d), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = which(colnames(d) == "MN"), number = Inf, sort.by = "t")
  rownames(tt)[tt$logFC != 0]

}

score <- function(Xtr, ytr, Xte, g) {

  d <- data.frame(y = as.integer(ytr == "MN"), t(Xtr[g, , drop = FALSE]))
  fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial))
  as.numeric(predict(fit, newdata = data.frame(t(Xte[g, , drop = FALSE])), type = "response"))

}

strat_folds <- function(yy, k) {

  f <- integer(length(yy))
  for (lv in levels(yy)) {
    idx <- which(yy == lv)
    f[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }

  f

}

cv_rows <- list()

for (rep_i in seq_len(N_REPEAT)) {

  f <- strat_folds(y, N_FOLD)

  for (k in seq_len(N_FOLD)) {

    tr <- f != k
    g_all <- rank_on(Xi[, tr, drop = FALSE], y[tr])

    for (K in K_GRID) {

      g <- head(g_all, K)
      pr <- score(Xi[, tr, drop = FALSE], y[tr], Xi[, !tr, drop = FALSE], g)
      cv_rows[[length(cv_rows) + 1]] <- data.frame(
        K = K, repeat_id = rep_i, fold = k, auc = auc_of(y[!tr], pr),
        stringsAsFactors = FALSE)

    }

  }

}

cv <- do.call(rbind, cv_rows)

cv_sum <- do.call(rbind, lapply(split(cv, cv$K), function(x) {
  data.frame(K = x$K[1], n_eval = nrow(x),
             cv_auc_mean = round(mean(x$auc), 4),
             cv_auc_sd = round(sd(x$auc), 4),
             cv_auc_lo = round(quantile(x$auc, 0.025), 4),
             cv_auc_hi = round(quantile(x$auc, 0.975), 4),
             stringsAsFactors = FALSE)
}))

# external: lock each K on the full discovery cohort, then test the validation cohorts
ext_rows <- list()
g_all <- rank_on(Xi, y)

for (K in K_GRID) {

  g <- head(g_all, K)
  pr_e <- score(Xi, y, Xe, g)
  pr_r <- score(Xi, y, Xr, g)

  ext_rows[[length(ext_rows) + 1]] <- data.frame(
    K = K,
    ercb_auc = round(auc_of(B$y$ercb, pr_e), 4),
    rna_auc = round(auc_of(B$y$rna, pr_r), 4),
    genes = paste(g, collapse = ";"), stringsAsFactors = FALSE)

}

ext <- do.call(rbind, ext_rows)

front <- merge(cv_sum, ext, by = "K", sort = FALSE)
front$locked_panel <- NA

# The locked four-gene panel as a reference point. Its numbers are READ from the tables that produced
# them (09a CV summary, 09b external) - never transcribed from memory or from a handoff note.
cv_locked <- read.csv(file.path(tab_dir, "09a_panel_cv_summary.csv"), stringsAsFactors = FALSE)
ext_locked <- read.csv(file.path(tab_dir, "09b_panel_external.csv"), stringsAsFactors = FALSE)

cv_l <- cv_locked[cv_locked$method == "lasso", ]
stopifnot(nrow(cv_l) == 1)

ercb_l <- ext_locked$auc[ext_locked$cohort == "ERCB external"]
rna_l  <- ext_locked$auc[ext_locked$cohort == "RNA-seq independent cohort"]
stopifnot(length(ercb_l) == 1, length(rna_l) == 1)

front <- rbind(
  front,
  data.frame(K = length(locked), n_eval = cv_l$n_folds,
             cv_auc_mean = cv_l$mean_auc, cv_auc_sd = cv_l$sd_auc,
             cv_auc_lo = NA, cv_auc_hi = NA,
             ercb_auc = ercb_l, rna_auc = rna_l,
             genes = paste(locked, collapse = ";"), locked_panel = "locked four-gene panel",
             stringsAsFactors = FALSE))

cat(sprintf("\nlocked panel read from 09a/09b: CV %.4f (sd %.4f, %d folds), ERCB %.4f, RNA-seq %.4f\n",
            cv_l$mean_auc, cv_l$sd_auc, cv_l$n_folds, ercb_l, rna_l))

write.csv(front, file.path(tab_dir, "10d_panel_size_frontier.csv"), row.names = FALSE)

cat("\n================ PANEL-SIZE FRONTIER ================\n")
print(front[, c("K", "cv_auc_mean", "cv_auc_sd", "ercb_auc", "rna_auc", "locked_panel")], row.names = FALSE)

# the top-K curve only: K == 4 also matches the locked-panel reference row, so filter it out
tk <- front[is.na(front$locked_panel), ]
tk <- tk[order(tk$K), ]
stopifnot(!any(duplicated(tk$K)))

cat(sprintf("\nTop-K curve - internal CV: K=2 %.3f | K=4 %.3f | K=12 %.3f | K=50 %.3f\n",
            tk$cv_auc_mean[tk$K == 2], tk$cv_auc_mean[tk$K == 4],
            tk$cv_auc_mean[tk$K == 12], tk$cv_auc_mean[tk$K == 50]))
cat(sprintf("Top-K curve - external ERCB: best K=%d (%.3f) | K=4 %.3f | K=50 %.3f (delta %.3f)\n",
            tk$K[which.max(tk$ercb_auc)], max(tk$ercb_auc),
            tk$ercb_auc[tk$K == 4], tk$ercb_auc[tk$K == 50],
            tk$ercb_auc[tk$K == 50] - tk$ercb_auc[tk$K == 4]))
cat(sprintf("Locked panel vs top-4 at the same K: CV %.4f vs %.4f (delta %+.4f); ERCB %.4f vs %.4f (%+.4f)\n",
            cv_l$mean_auc, tk$cv_auc_mean[tk$K == 4], cv_l$mean_auc - tk$cv_auc_mean[tk$K == 4],
            ercb_l, tk$ercb_auc[tk$K == 4], ercb_l - tk$ercb_auc[tk$K == 4]))

# --- figure ----------------------------------------------------------------

long <- rbind(
  data.frame(K = front$K, auc = front$cv_auc_mean, series = "internal CV (nested)"),
  data.frame(K = front$K, auc = front$ercb_auc, series = "external ERCB"),
  data.frame(K = front$K, auc = front$rna_auc, series = "external RNA-seq"))

p <- ggplot(long, aes(x = K, y = auc, colour = series)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_vline(xintercept = length(locked), linetype = "dashed", colour = "grey45") +
  annotate("text", x = length(locked), y = min(long$auc, na.rm = TRUE),
           label = "locked panel (4)", hjust = -0.08, vjust = 0, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = K_GRID) +
  labs(title = "Panel-size frontier: what does parsimony cost?",
       subtitle = "Genes ranked inside every fold (no leakage); external bars lock each K on full discovery",
       x = "genes in the panel (K)", y = "AUC (MN vs other)", colour = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")

print(p)

pdf(file.path(fig_dir, "10d_panel_size_frontier.pdf"), width = 6.8, height = 4.6)
print(p)
dev.off()

jpeg(file.path(fig_dir, "10d_panel_size_frontier.jpeg"), width = 6.8, height = 4.6,
     units = "in", res = 300, quality = 95, bg = "white")
print(p)
dev.off()

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 10d_panel_size_frontier.R"),
             capture.output(print(front[, c("K", "cv_auc_mean", "ercb_auc", "rna_auc")], row.names = FALSE))),
           file.path(PROJ, "logs", "10d_panel_size_frontier.log"))

cat("\nDone.\n")
