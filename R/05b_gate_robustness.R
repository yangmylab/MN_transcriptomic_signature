#!/usr/bin/env Rscript
# Paper 1 - robustness checks for the separability gate (05a).
#
# 05a showed that MN-vs-other-GN is learnable, but its external AUC is only trustworthy because the
# validation cohort is standardised with TRAINING statistics (F07, audit 2026-10-20). This script
# applies the same rule: the scaling is fitted on the training matrix and applied unchanged to the
# test cohort, so these robustness numbers describe the same model that 05a reports. Cohort-internal
# z-scoring (the previous version) is a transductive step and must never come back - see P-58.
#
# Three objections have to be answered before the gate number is used to justify Phase D/E:
#
#   1. Is it better than chance?           -> permutation null (labels shuffled in the training set)
#   2. Is it a QC artefact?                -> recompute after removing the flagged ERCB samples
#                                             (71 of 170 non-control samples were flagged in Phase A;
#                                             the frozen policy keeps them, this is a sensitivity test)
#   3. Is it symmetric?                    -> train on ERCB, test on the NEPTUNE discovery cohort
#
# plus a per-comparator breakdown (which glomerulopathy is the hard one) using disease_raw.
#
# Outputs
#   output/tables/gate_robustness.csv
#   output/tables/gate_auc_by_comparator.csv
#   output/models_gate/05b_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(limma); library(glmnet); library(pROC); library(ggplot2) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_gate")

set.seed(20261018)

K <- NA          # all common genes, the best external configuration in 05a
ALPHA <- 0.5
B_PERM <- 20

save_plot <- function(p, stem, width = 7.5, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

mg <- function(set) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  keep <- ann$disease != "UNKNOWN" & grp != "Control"
  list(ms = ms[, keep, drop = FALSE], grp = grp[keep], raw = ann$disease_raw[keep],
       gsm = ann$gsm[keep])
}

# Training statistics only: the scaling is estimated on the training matrix and applied unchanged to
# the test cohort, which therefore never contributes to its own standardisation (F07, P-58).
std_fit <- function(tr_ms) {

  mu <- rowMeans(tr_ms)
  sdv <- apply(tr_ms, 1, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  list(mu = mu, sd = sdv)

}

std <- function(ms, sc) {

  z <- sweep(sweep(ms, 1, sc$mu), 1, sc$sd, "/")
  z[!is.finite(z)] <- 0

  z

}

# Cross-platform harmonisation, required before any scaling: discovery is HuGene-2.1-ST, ERCB is
# HG-U133. Without this the typical ERCB gene sits at |z| ~ 4 against discovery sd, the linear
# predictor saturates and the "transfer" measures platform shift (audit 2026-11-18, F62). Rank-based
# mapping onto the discovery reference distribution uses no labels.
quantile_map <- function(ms, ref) {

  out <- ms

  for (i in seq_len(nrow(ms))) {

    v <- ms[i, ]
    ok <- is.finite(v)

    if (sum(ok) < 3) { out[i, ] <- NA_real_; next }

    r <- rank(v[ok], ties.method = "average")
    idx <- pmin(ncol(ref), pmax(1L, round(r / (sum(ok) + 1) * ncol(ref))))

    out[i, ] <- NA_real_
    out[i, ok] <- ref[i, idx]

  }

  out

}

rank_genes <- function(X, y, k) {
  d <- model.matrix(~ 0 + factor(y))
  colnames(d) <- levels(factor(y))
  fit <- eBayes(lmFit(X, d), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = which(colnames(d) == "MN"), number = Inf, sort.by = "t")
  tt <- tt[tt$logFC != 0, ]
  g <- rownames(tt)
  if (is.na(k)) g else head(g, min(k, length(g)))
}

run_fit <- function(Mtr, ytr, Mte) {

  sc <- std_fit(Mtr)
  Xtr <- std(Mtr, sc)
  Xte <- std(Mte, sc)

  ybin <- factor(ifelse(ytr == "MN", "MN", "rest"), levels = c("rest", "MN"))
  g <- rank_genes(Xtr, ytr, K)
  cv <- cv.glmnet(t(Xtr[g, , drop = FALSE]), ybin, family = "binomial", alpha = ALPHA,
                  nfolds = min(5, min(table(ybin))), type.measure = "auc")
  pr <- as.numeric(predict(cv$glmnet.fit, newx = t(Xte[g, , drop = FALSE]),
                           s = cv$lambda.min, type = "response"))

  list(pr = pr, genes = g, lambda = cv$lambda.min)

}

auc_of <- function(pr, y) {
  ybin <- factor(ifelse(y == "MN", "MN", "rest"), levels = c("rest", "MN"))
  as.numeric(auc(roc(ybin, pr, direction = "<", quiet = TRUE)))
}

# --- data --------------------------------------------------------------------

tr <- mg("disc_glom_array")
ex <- mg("val_glom_array")

genes <- intersect(rownames(tr$ms), rownames(ex$ms))

M_tr <- tr$ms[genes, , drop = FALSE]
M_ex <- quantile_map(ex$ms[genes, , drop = FALSE], t(apply(M_tr, 1, sort)))

# --- reference model and its AUC ---------------------------------------------

ref <- run_fit(M_tr, tr$grp, M_ex)
auc_ref <- auc_of(ref$pr, ex$grp)

cat("reference model: AUC", round(auc_ref, 4), "| genes used", length(ref$genes), "\n")

# reconcile with 05a: the K = "all" ERCB transfer is the same model, so the two AUCs must agree.
ext_file <- file.path(tab_dir, "separability_external.csv")

if (file.exists(ext_file)) {

  ext <- read.csv(ext_file, stringsAsFactors = FALSE)
  ext_ref <- ext$auc[ext$K == "all" & grepl("ERCB", ext$cohort)]

  if (length(ext_ref) == 1) {

    cat("reconciliation vs 05a (K = all, ERCB):", round(ext_ref, 4),
        "| delta", round(auc_ref - ext_ref, 4), "\n")
    if (abs(auc_ref - ext_ref) > 0.05) warning("05b and 05a disagree by more than 0.05 AUC")

  }

}


# --- 1. permutation null -----------------------------------------------------

null <- numeric(B_PERM)

for (b in seq_len(B_PERM)) {
  yperm <- sample(tr$grp)
  fit_b <- run_fit(M_tr, yperm, M_ex)
  null[b] <- auc_of(fit_b$pr, ex$grp)
  cat("  permutation", b, "AUC", round(null[b], 4), "\n")
}

p_emp <- (1 + sum(null >= auc_ref)) / (1 + B_PERM)

# --- 2. QC sensitivity -------------------------------------------------------

qc_flags <- read.csv(file.path(OUT, "qc", "03c_val_glom_array_sample_qc.csv"),
                     stringsAsFactors = FALSE)

clean_ex <- !(ex$gsm %in% qc_flags$gsm[qc_flags$flagged])

auc_clean <- auc_of(ref$pr[clean_ex], ex$grp[clean_ex])

tr_flags <- read.csv(file.path(OUT, "qc", "03c_disc_glom_array_sample_qc.csv"),
                     stringsAsFactors = FALSE)

clean_tr <- !(tr$gsm %in% tr_flags$gsm[tr_flags$flagged])

fit_clean <- run_fit(M_tr[, clean_tr, drop = FALSE], tr$grp[clean_tr], M_ex[, clean_ex, drop = FALSE])
auc_both <- auc_of(fit_clean$pr, ex$grp[clean_ex])

# --- 3. reverse direction ----------------------------------------------------

rev <- run_fit(M_ex, ex$grp, M_tr)
auc_rev <- auc_of(rev$pr, tr$grp)

# --- per-comparator AUCs -----------------------------------------------------

# ERCB labels MN as "MGN" in disease_raw, so MN status comes from the harmonised group and the
# comparator breakdown keeps every ERCB diagnosis separate (IgA, SLE, RPGN, DN, HT, TMD, FSGS_MCD).
is_mn <- ex$grp == "MN"

raw <- ex$raw
raw[is.na(raw) | raw == ""] <- "Other"
comparators <- sort(setdiff(unique(raw[!is_mn]), c("MGN", "MN")))

by_comp <- do.call(rbind, lapply(comparators, function(cp) {
  sel <- is_mn | raw == cp
  data.frame(comparator = cp, n = sum(raw[sel] == cp),
             auc = auc_of(ref$pr[sel], ex$grp[sel]), stringsAsFactors = FALSE)
}))

by_comp$comparator <- factor(by_comp$comparator, levels = by_comp$comparator[order(by_comp$auc)])

# --- report ------------------------------------------------------------------

rob <- data.frame(
  check = c("reference (discovery -> ERCB)", "permutation null: median AUC",
            "permutation null: 95th percentile", "empirical p",
            "ERCB QC-flagged samples removed", "QC-flagged removed in both cohorts",
            "reverse direction (ERCB -> discovery)"),
  value = c(round(auc_ref, 4), round(median(null), 4), round(quantile(null, 0.95), 4),
            signif(p_emp, 3), round(auc_clean, 4), round(auc_both, 4), round(auc_rev, 4)),
  n_test = c(sum(ex$grp != ""), NA, NA, NA, sum(clean_ex), sum(clean_ex), length(tr$grp)),
  stringsAsFactors = FALSE)

write.csv(rob, file.path(tab_dir, "gate_robustness.csv"), row.names = FALSE)
write.csv(by_comp, file.path(tab_dir, "gate_auc_by_comparator.csv"), row.names = FALSE)

p <- ggplot(data.frame(auc = c(null, auc_ref), type = c(rep("label-shuffled", B_PERM),
                                                         "observed")),
            aes(auc, fill = type)) +
  geom_histogram(bins = 12, position = "identity", alpha = 0.75, colour = "white") +
  geom_vline(xintercept = auc_ref, colour = "#C0392B", linewidth = 0.9) +
  labs(title = "Permutation null for the MN-vs-other-GN gate",
       subtitle = paste0("observed AUC ", round(auc_ref, 3), ", empirical p = ", signif(p_emp, 2),
                         " (", B_PERM, " shuffles)"),
       x = "AUC in ERCB", y = "shuffles", fill = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "05b_gate_permutation_null"), width = 7, height = 4.8)

p <- ggplot(by_comp, aes(auc, comparator)) +
  geom_col(fill = "#2471A3", width = 0.65) +
  geom_text(aes(label = paste0("AUC ", round(auc, 3), "  (n=", n, ")")), hjust = -0.1, size = 3.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
  xlim(0, 1.12) +
  labs(title = "MN vs each comparator in ERCB",
       subtitle = "MN probability from the discovery-trained elastic net",
       x = "AUC (MN as case)", y = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_plot(p, file.path(fig_dir, "05b_gate_auc_by_comparator"), width = 7.5, height = 4.8)

cat("\n================ ROBUSTNESS ================\n")
print(rob, row.names = FALSE)

cat("\n================ AUC BY COMPARATOR (ERCB) ================\n")
print(by_comp, row.names = FALSE)

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 05b_gate_robustness.R"),
             capture.output(print(rob, row.names = FALSE)),
             capture.output(print(by_comp, row.names = FALSE))),
           file.path(PROJ, "logs", "gate_robustness.log"))

cat("\nDone.\n")
