#!/usr/bin/env Rscript
# Track B - the attrition ladder: how much of a published MN panel's AUC survives rigour? (2026-11-18)
#
# For every panel in metadata/trackB_panel_register.csv, the SAME gene list is carried down six rungs
# of increasing rigour. Nothing is re-selected: the published genes are fixed, so each drop is
# attributable to the validation design, not to a different feature set.
#
#   R1 apparent      fit on all discovery samples, own-cohort scaling, score the same samples
#   R2 + dedup       one sample per patient (patient-level de-duplication)
#   R3 + nested CV   repeated stratified CV inside the deduplicated discovery cohort
#   R4 + transductive external   train -> ERCB, but ERCB z-scored against ITS OWN mean/sd
#                                (the convention most papers leave unstated)
#   R5 + harmonised external     train -> ERCB after quantile mapping onto the discovery reference
#                                (the valid cross-platform transfer, F62)
#   R6 + permutation null        refit null for the R5 ERCB AUC (200 label shuffles)
#
# R0 is the AUC as reported by the paper; it is transcribed, never recomputed.
#
# Inputs   metadata/trackB_panel_register.csv, work/bench/10b_harmonised.rds
# Outputs  output/tables/10c_attrition.csv, output/tables/10c_panel_summary.csv
#          output/models_trackb/10c_attrition.{pdf,jpeg}
# parallel: serial - fits are 4-7-predictor logistic models; fork overhead would dominate

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(pROC); library(ggplot2) })

tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_trackb")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20261118)

B <- readRDS(file.path(PROJ, "work", "bench", "10b_harmonised.rds"))
reg <- read.csv(file.path(PROJ, "metadata", "trackB_panel_register.csv"), stringsAsFactors = FALSE)

N_PERM <- 200L
N_REPEAT <- 10L
N_FOLD <- 5L

auc_of <- function(y, s) {

  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(auc(roc(y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))

}

ci_of <- function(y, s) {

  r <- roc(y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  c(ci[1], ci[3])

}

# fit a fixed-gene logistic model on the training matrix and score the test matrix
fit_score <- function(Xtr, ytr, Xte) {

  d <- data.frame(y = as.integer(ytr == "MN"), t(Xtr))
  fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial))
  as.numeric(predict(fit, newdata = data.frame(t(Xte)), type = "response"))

}

# stratified folds, fixed seed so every panel sees the same partitioning
strat_folds <- function(y, k) {

  f <- integer(length(y))
  for (lv in levels(y)) {
    idx <- which(y == lv)
    f[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }

  f

}

rows <- list()

for (p in seq_len(nrow(reg))) {

  pid <- reg$panel_id[p]
  genes_pub <- trimws(strsplit(reg$genes[p], ";")[[1]])
  present <- intersect(genes_pub, B$genes)
  absent <- setdiff(genes_pub, B$genes)

  cat(sprintf("\n%s %s | reported: %s\n", pid, reg$first_author_year[p], reg$reported_auc_text[p]))
  cat(sprintf("   genes: %d published, %d measurable, %d absent (%s)\n",
              length(genes_pub), length(present), length(absent),
              if (length(absent)) paste(absent, collapse = ",") else "none"))

  add <- function(rung, label, auc, lo, hi, n, n_mn, note = "") {

    rows[[length(rows) + 1]] <<- data.frame(
      panel_id = pid, year = reg$first_author_year[p], pmid = reg$pmid[p],
      rung = rung, rung_label = label,
      n_genes_published = length(genes_pub), n_genes_used = length(present),
      n = n, n_MN = n_mn,
      auc = round(auc, 4), ci_low = round(lo, 4), ci_high = round(hi, 4),
      note = note, stringsAsFactors = FALSE)

  }

  yd <- B$y$disc
  Xd_trans <- B$z_transductive$disc[present, , drop = FALSE]
  Xd_harm  <- B$z_harmonised$disc[present, , drop = FALSE]
  Xe_trans <- B$z_transductive$ercb[present, , drop = FALSE]
  Xe_harm  <- B$z_harmonised$ercb[present, , drop = FALSE]
  Xr_harm  <- B$z_harmonised$rna[present, , drop = FALSE]

  # --- R1 apparent ---------------------------------------------------------
  s <- fit_score(Xd_trans, yd, Xd_trans)
  ci <- ci_of(yd, s)
  add("R1", "apparent (own cohort, own scaling)", auc_of(yd, s), ci[1], ci[2],
      length(yd), sum(yd == "MN"))

  # --- R2 + patient de-duplication -----------------------------------------
  keep <- !duplicated(B$patient$disc)
  Xd2 <- Xd_trans[, keep, drop = FALSE]
  yd2 <- yd[keep]
  s2 <- fit_score(Xd2, yd2, Xd2)
  ci2 <- ci_of(yd2, s2)
  add("R2", "+ patient de-duplication", auc_of(yd2, s2), ci2[1], ci2[2],
      length(yd2), sum(yd2 == "MN"),
      sprintf("%d samples collapsed to %d patients", length(yd), length(yd2)))

  # --- R3 + repeated stratified CV -----------------------------------------
  cv_auc <- numeric(0)
  for (rep_i in seq_len(N_REPEAT)) {
    f <- strat_folds(yd2, N_FOLD)
    for (k in seq_len(N_FOLD)) {
      tr <- f != k
      pr <- fit_score(Xd2[, tr, drop = FALSE], yd2[tr], Xd2[, !tr, drop = FALSE])
      a <- auc_of(yd2[!tr], pr)
      if (!is.na(a)) cv_auc <- c(cv_auc, a)
    }
  }
  add("R3", "+ repeated 5-fold CV (selection fixed)", mean(cv_auc), quantile(cv_auc, 0.025),
      quantile(cv_auc, 0.975), length(yd2), sum(yd2 == "MN"),
      sprintf("%d evaluations, sd = %.3f", length(cv_auc), sd(cv_auc)))

  # --- R4 transductive external -------------------------------------------
  s4 <- fit_score(Xd2, yd2, Xe_trans)
  ci4 <- ci_of(B$y$ercb, s4)
  add("R4", "+ external ERCB, cohort-internal scaling", auc_of(B$y$ercb, s4), ci4[1], ci4[2],
      length(B$y$ercb), sum(B$y$ercb == "MN"), "the convention most papers leave unstated")

  # --- R5 harmonised external ---------------------------------------------
  s5 <- fit_score(Xd2, yd2, Xe_harm)
  ci5 <- ci_of(B$y$ercb, s5)
  obs <- auc_of(B$y$ercb, s5)
  add("R5", "+ external ERCB, harmonised transfer", obs, ci5[1], ci5[2],
      length(B$y$ercb), sum(B$y$ercb == "MN"), "valid cross-platform transfer (F62)")

  # RNA-seq replication at the same rung
  s5r <- fit_score(Xd2, yd2, Xr_harm)
  ci5r <- ci_of(B$y$rna, s5r)
  add("R5r", "+ independent RNA-seq cohort", auc_of(B$y$rna, s5r), ci5r[1], ci5r[2],
      length(B$y$rna), sum(B$y$rna == "MN"), "same model, third sample source")

  # --- R6 refit permutation null ------------------------------------------
  null <- vapply(seq_len(N_PERM), function(i) {
    ysh <- sample(yd2)
    pr <- fit_score(Xd2, ysh, Xe_harm)
    auc_of(B$y$ercb, pr)
  }, numeric(1))

  p_emp <- (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
  add("R6", "+ refit permutation null (ERCB)", p_emp, NA, NA,
      length(B$y$ercb), sum(B$y$ercb == "MN"),
      sprintf("empirical p, null median %.3f, q95 %.3f",
              median(null, na.rm = TRUE), quantile(null, 0.95, na.rm = TRUE)))

  cat(sprintf("   R1 %.3f -> R3 %.3f -> R5 %.3f (RNA-seq %.3f), R6 p = %.3f\n",
              auc_of(yd, s), mean(cv_auc), obs, auc_of(B$y$rna, s5r), p_emp))

}

att <- do.call(rbind, rows)

# R6 stores a p-value rather than an AUC; keep the column honest
att$metric <- ifelse(att$rung == "R6", "empirical_p", "auc")

write.csv(att, file.path(tab_dir, "10c_attrition.csv"), row.names = FALSE)

# --- per-panel summary: where the drop happens ------------------------------

wide <- reshape(att[att$metric == "auc", c("panel_id", "rung", "auc")],
                idvar = "panel_id", timevar = "rung", direction = "wide")
names(wide) <- sub("^auc\\.", "", names(wide))
pv <- att[att$rung == "R6", c("panel_id", "auc")]
names(pv) <- c("panel_id", "perm_p")

summ <- merge(wide, pv, by = "panel_id", sort = FALSE)
summ$drop_R1_R5 <- round(summ$R1 - summ$R5, 4)
summ$drop_R3_R4 <- round(summ$R3 - summ$R4, 4)
summ$drop_R4_R5 <- round(summ$R4 - summ$R5, 4)

write.csv(summ, file.path(tab_dir, "10c_panel_summary.csv"), row.names = FALSE)

cat("\n================ PER-PANEL ATTRITION ================\n")
print(summ, row.names = FALSE)
cat(sprintf("\nmedian R1 %.3f | R3 %.3f | R5 %.3f | median drop R1->R5 %.3f | panels with R6 p < 0.05: %d/%d\n",
            median(summ$R1), median(summ$R3), median(summ$R5),
            median(summ$drop_R1_R5), sum(summ$perm_p < 0.05), nrow(summ)))

# --- figure ----------------------------------------------------------------

plot_df <- att[att$rung %in% c("R1", "R2", "R3", "R4", "R5", "R5r"), ]
plot_df$rung <- factor(plot_df$rung, levels = c("R1", "R2", "R3", "R4", "R5", "R5r"))

p <- ggplot(plot_df, aes(x = rung, y = auc, group = panel_id, colour = panel_id)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  labs(title = "Attrition of published MN biomarker panels under increasing validation rigour",
       subtitle = "Same published gene list at every rung; nothing re-selected",
       x = NULL, y = "AUC (MN vs other)", colour = "panel") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

print(p)

pdf(file.path(fig_dir, "10c_attrition.pdf"), width = 7.2, height = 4.6)
print(p)
dev.off()

jpeg(file.path(fig_dir, "10c_attrition.jpeg"), width = 7.2, height = 4.6,
     units = "in", res = 300, quality = 95, bg = "white")
print(p)
dev.off()

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 10c_panel_reproduction.R"),
             capture.output(print(summ, row.names = FALSE))),
           file.path(PROJ, "logs", "10c_panel_reproduction.log"))

cat("\nDone.\n")
