#!/usr/bin/env Rscript
# Paper 1 - Phase E step 5: the four remaining sensitivity analyses.
#
#   A. PER-COMPARATOR AUC  - the gate reports per-comparator AUCs, the panel did not. Reviewers ask
#      how a marker behaves against FSGS specifically, not against a pooled "other GN" class.
#   B. PIPELINE INDEPENDENCE - the five genes are re-scored in the submitters' own processed matrix
#      (GSE200828_processed_data.txt.gz: Entrez-keyed, disease encoded in the column names). If the
#      per-gene MN-vs-other separation reappears there, the signal is not an artefact of our RMA /
#      Brainarray pipeline. MYOZ2 is the gene that most needs this check.
#   C. TECHNICAL CONFOUND  - does the panel score track expression PCs or ERCB batch rather than
#      disease? Reported as R2 of score ~ PCs and as AUC before/after adjusting for PCs + batch.
#   D. FORMAL "WITHOUT SLAMF8" MODEL - the leave-one-gene-out in 09d shows that dropping SLAMF8
#      (detected in only 9.5% of RNA-seq samples) gives a model with no gene that two platforms fail
#      to measure. Here that model is reported properly: external AUC with CI, a permutation null and
#      recalibration, so it can be footnoted honestly. The gene set is derived as
#      setdiff(panel, "SLAMF8"), so its size follows whatever 09a locks (it was 4 genes when the
#      panel was 5; after the 2026-10-20 audit re-run the panel is 4, so this model is 3 genes).
#
# Outputs
#   output/tables/09e_panel_by_comparator.csv
#   output/tables/09e_authors_matrix_check.csv
#   output/tables/09e_technical_confound.csv
#   output/tables/09e_no_slamf8_model.csv
#   output/models_panel/09e_{by_comparator,authors_check,no_slamf8_roc}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(pROC); library(ggplot2); library(org.Hs.eg.db) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_panel")
SEED <- 20261019L
set.seed(SEED)

save_plot <- function(p, stem, width = 8, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

panel <- read.csv(file.path(tab_dir, "09a_panel_consensus.csv"), stringsAsFactors = FALSE)$gene

load_set <- function(set, keep_control = FALSE) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(ms), ann$gsm))
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  keep <- if (keep_control) grp != "UNKNOWN" else !(grp %in% c("Control", "UNKNOWN"))
  # `raw` keeps the per-diagnosis labels (DN, SLE, IgAN, ...) for the comparator analysis; `group`
  # is the collapsed 5-class vocabulary and must NOT keep unmapped labels (they would become NA)
  raw <- grp[keep]
  collapsed <- ifelse(raw %in% c("MN", "FSGS", "MCD", "Control"), raw, "Other")
  list(ms = ms[, keep, drop = FALSE],
       y = factor(ifelse(raw == "MN", "MN", "rest"), levels = c("rest", "MN")),
       group = factor(collapsed, levels = c("Control", "FSGS", "MCD", "Other", "MN")),
       raw = raw)
}

disc  <- load_set("disc_glom_array")
ercb  <- load_set("val_glom_array")
rna   <- load_set("cross_glom_rnaseq")
lcm   <- load_set("pla2r_lcm")

mu  <- colMeans(t(disc$ms[panel, , drop = FALSE]))
sdv <- apply(t(disc$ms[panel, , drop = FALSE]), 2, sd)
sdv[!is.finite(sdv) | sdv == 0] <- 1

z <- function(o) {
  X <- t(o$ms)
  out <- matrix(0, nrow(X), length(panel), dimnames = list(NULL, panel))
  have <- intersect(panel, colnames(X))
  if (length(have)) {
    zz <- sweep(sweep(X[, have, drop = FALSE], 2, mu[have]), 2, sdv[have], "/")
    zz[!is.finite(zz)] <- 0
    out[, have] <- zz
  }
  out
}

fit5 <- glm(disc$y ~ ., data = data.frame(z(disc), y = disc$y), family = binomial)
sc5 <- function(o) as.numeric(predict(fit5, newdata = as.data.frame(z(o)), type = "response"))

# --- A. per-comparator AUC ------------------------------------------------
rows <- list()
for (nm in c("ERCB", "RNA-seq")) {
  o <- if (nm == "ERCB") ercb else rna
  s <- sc5(o)
  for (g in setdiff(unique(o$raw), c("MN", "Control-LD", "Control-TN"))) {
    idx <- o$raw %in% c("MN", g)
    if (sum(o$raw == g) < 3) next
    r <- roc(o$y[idx], s[idx], levels = c("rest", "MN"), direction = "<", quiet = TRUE)
    ci <- as.numeric(ci.auc(r, method = "delong"))
    rows[[paste(nm, g)]] <- data.frame(
      cohort = nm, comparator = as.character(g), n_MN = sum(o$raw[idx] == "MN"),
      n_comparator = sum(o$raw[idx] == g),
      auc = round(as.numeric(auc(r)), 4), ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
      stringsAsFactors = FALSE)
  }
}
bycomp <- do.call(rbind, rows)
bycomp <- bycomp[order(bycomp$cohort, bycomp$auc), ]
write.csv(bycomp, file.path(tab_dir, "09e_panel_by_comparator.csv"), row.names = FALSE)
cat("--- per-comparator AUC ---\n"); print(bycomp, row.names = FALSE)

bc <- rbind(
  data.frame(cohort = bycomp$cohort, label = paste0(bycomp$cohort, ": ", bycomp$comparator),
             auc = bycomp$auc, ci_low = bycomp$ci_low, ci_high = bycomp$ci_high))
bc$label <- factor(bc$label, levels = bc$label[order(bc$auc)])

# data-driven axis window: never crop a CI whisker (audit F37)
bc_rng <- range(c(bc$ci_low, bc$auc, bc$ci_high))
bc_ylim <- c(max(0, bc_rng[1] - 0.03), min(1.02, bc_rng[2] + 0.02))

save_plot(ggplot(bc, aes(x = label, y = auc, fill = cohort)) +
            geom_col(width = 0.68) +
            geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
            geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
            # Limits are computed from the data, not fixed at c(0.3, 1): with the re-derived panel the
            # FSGS_MCD lower CI is 0.194 and several upper bounds sit exactly on 1.0, either of which
            # a hard ylim cropped (audit F37 — whiskers were cut off in the composed plate).
            coord_flip(ylim = bc_ylim) +
            scale_fill_brewer(palette = "Set2") +
            labs(title = "Locked panel against each comparator diagnosis",
                 subtitle = "MN vs one comparator at a time; FSGS is the clinically hard case",
                 x = NULL, y = "AUC (95% CI)", fill = NULL) +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), legend.position = "bottom",
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09e_by_comparator"), width = 7.4, height = 5.2)

# --- B. pipeline independence (authors' processed matrix) -----------------
authors_f <- file.path(DATA, "GSE200828", "GSE200828_processed_data.txt.gz")
if (file.exists(authors_f)) {
  ids <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = panel,
                                                keytype = "SYMBOL", columns = "ENTREZID"))
  ids <- ids[!is.na(ids$ENTREZID), ]
  d <- data.table::fread(authors_f, sep = "\t", header = TRUE, showProgress = FALSE)
  d <- as.data.frame(d)
  sub <- d[d$ENTREZ_GENE_ID %in% ids$ENTREZID, ]
  cat("authors' matrix rows for panel:", nrow(sub), "of", length(panel), "genes\n")

  smp <- setdiff(colnames(d), c("ENTREZ_GENE_ID", "Expressed"))
  lab <- ifelse(grepl("_MN$", smp), "MN",
         ifelse(grepl("_LD|_TN|_Control", smp), "Control", "Other"))

  chk <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
    v <- as.numeric(sub[i, smp])
    sym <- ids$SYMBOL[match(sub$ENTREZ_GENE_ID[i], ids$ENTREZID)]
    keep <- lab %in% c("MN", "Other")
    yy <- factor(ifelse(lab[keep] == "MN", "MN", "rest"), levels = c("rest", "MN"))
    a_auth <- if (length(unique(yy)) == 2 && sd(v[keep]) > 0) {
      as.numeric(auc(roc(yy, v[keep], levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
    } else NA_real_
    # our own RMA for the same gene in the discovery cohort
    v_ours <- as.numeric(disc$ms[sym, ])
    a_ours <- as.numeric(auc(roc(disc$y, v_ours, levels = c("rest", "MN"), direction = "<",
                                 quiet = TRUE)))
    data.frame(gene = sym, authors_auc = round(a_auth, 4), our_rma_auc = round(a_ours, 4),
               difference = round(a_auth - a_ours, 4), stringsAsFactors = FALSE)
  }))
  chk <- chk[order(-chk$authors_auc), ]
  write.csv(chk, file.path(tab_dir, "09e_authors_matrix_check.csv"), row.names = FALSE)
  cat("\n--- authors' matrix vs our RMA (single-gene MN-vs-other AUC) ---\n"); print(chk, row.names = FALSE)

  chk_df <- rbind(
    data.frame(gene = chk$gene, source = "submitters' matrix", auc = chk$authors_auc),
    data.frame(gene = chk$gene, source = "our RMA + Brainarray", auc = chk$our_rma_auc))
  # data-driven window (same F37 defect class as the comparator panel): a fixed c(0.4, 1) cropped
  # any gene whose AUC fell below 0.4.
  chk_ylim <- c(max(0, min(chk_df$auc) - 0.05), 1)
  save_plot(ggplot(chk_df, aes(x = reorder(gene, -auc), y = auc, fill = source)) +
              geom_col(position = position_dodge(0.72), width = 0.66) +
              geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
              coord_cartesian(ylim = chk_ylim) +
              scale_fill_manual(values = c("submitters' matrix" = "grey45",
                                           "our RMA + Brainarray" = "firebrick")) +
              labs(title = "Panel genes are not an artefact of our preprocessing",
                   subtitle = "MN vs other-GN AUC computed independently in each matrix",
                   x = NULL, y = "AUC", fill = NULL) +
              theme_bw(base_size = 10) +
              theme(panel.grid = element_blank(), legend.position = "bottom",
                    plot.title = element_text(face = "bold")),
            file.path(fig_dir, "09e_authors_check"), width = 7, height = 4.4)
} else {
  cat("authors' matrix not found at", authors_f, "- pipeline-independence check skipped\n")
}

# --- C. technical confound -------------------------------------------------
tech_rows <- list()
for (nm in c("disc", "ercb")) {
  o <- if (nm == "disc") disc else ercb
  X <- t(o$ms)
  # ERCB carries union-filled NA gene columns from the two-platform release: drop them before PCA
  X <- X[, colSums(is.na(X)) == 0, drop = FALSE]
  v <- apply(X, 2, var)
  v[!is.finite(v)] <- 0
  top <- order(v, decreasing = TRUE)[1:min(2000, sum(v > 0))]
  X <- X[, top, drop = FALSE]
  X[!is.finite(X)] <- 0
  pc <- prcomp(X, scale. = TRUE, center = TRUE)
  pcs <- pc$x[, 1:5]
  s <- sc5(o)
  r2 <- summary(lm(scale(s) ~ pcs))$r.squared
  auc_raw <- as.numeric(auc(roc(o$y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
  adj <- glm(o$y ~ s + pcs, family = binomial)
  auc_adj <- as.numeric(auc(roc(o$y, fitted(adj), levels = c("rest", "MN"), direction = "<",
                               quiet = TRUE)))
  tech_rows[[nm]] <- data.frame(cohort = nm, n = length(s),
                                r2_score_on_PC1_5 = round(r2, 4),
                                auc_unadjusted = round(auc_raw, 4),
                                auc_adjusted_for_PC1_5 = round(auc_adj, 4),
                                MN_term_p = signif(summary(adj)$coefficients["s", 4], 4),
                                stringsAsFactors = FALSE)
}
tech <- do.call(rbind, tech_rows)
write.csv(tech, file.path(tab_dir, "09e_technical_confound.csv"), row.names = FALSE)
cat("\n--- technical confound ---\n"); print(tech, row.names = FALSE)

# --- D. formal "without SLAMF8" model (size follows the locked panel) ------
genes4 <- setdiff(panel, "SLAMF8")
z4 <- function(o) z(o)[, genes4, drop = FALSE]
fit4 <- glm(disc$y ~ ., data = data.frame(z4(disc), y = disc$y), family = binomial)
sc4 <- function(o) as.numeric(predict(fit4, newdata = as.data.frame(z4(o)), type = "response"))

four <- do.call(rbind, lapply(list(
  list(nm = "Discovery (apparent)", o = disc), list(nm = "ERCB external", o = ercb),
  list(nm = "RNA-seq independent", o = rna), list(nm = "LCM GSE216841", o = lcm)), function(z) {
  s <- sc4(z$o)
  r <- roc(z$o$y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  cal <- glm(z$o$y ~ s, family = binomial)
  data.frame(model = sprintf("%d-gene (%s)", length(genes4), paste(genes4, collapse = ", ")),
             cohort = z$nm,
             n = length(s), n_MN = sum(z$o$y == "MN"),
             auc = round(as.numeric(auc(r)), 4), ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
             calibration_in_large = round(mean(s) - mean(z$o$y == "MN"), 4),
             calibration_slope = round(coef(cal)[2], 4), stringsAsFactors = FALSE)
}))

s5 <- sc5(ercb); s4 <- sc4(ercb)
obs4 <- as.numeric(auc(roc(ercb$y, s4, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
null4 <- vapply(seq_len(200), function(i) {
  yy <- sample(ercb$y)
  as.numeric(auc(roc(yy, s4, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
}, numeric(1))
four$permutation_p_ERCB <- NA_real_
four$permutation_p_ERCB[four$cohort == "ERCB external"] <-
  signif((sum(null4 >= obs4) + 1) / (length(null4) + 1), 4)
write.csv(four, file.path(tab_dir, "09e_no_slamf8_model.csv"), row.names = FALSE)
cat(sprintf("\n--- formal %d-gene model without SLAMF8 (%s) ---\n",
            length(genes4), paste(genes4, collapse = ", "))); print(four, row.names = FALSE)

# --- E. cross-platform harmonisation sensitivity ---------------------------
# The gate (05a/05b) harmonises each validation cohort onto the discovery reference distribution
# before applying training statistics, because discovery and ERCB sit on different array platforms
# (F62, audit 2026-11-18). The panel has always used raw discovery mu/sd, which is why its calibration
# slope is extreme. The question here is whether that choice changes the panel ANSWER: it must not,
# otherwise the panel would depend on the standardisation convention.

ref_panel <- t(apply(disc$ms[panel, , drop = FALSE], 1, sort))

harmonise_panel <- function(o) {

  present <- intersect(panel, rownames(o$ms))
  M <- o$ms[present, , drop = FALSE]
  ref <- ref_panel[present, , drop = FALSE]
  out <- M

  for (i in seq_len(nrow(M))) {

    v <- M[i, ]
    ok <- is.finite(v)

    if (sum(ok) < 3) { out[i, ] <- NA_real_; next }

    r <- rank(v[ok], ties.method = "average")
    idx <- pmin(ncol(ref), pmax(1L, round(r / (sum(ok) + 1) * ncol(ref))))

    out[i, ] <- NA_real_
    out[i, ok] <- ref[i, idx]

  }

  out

}

zh <- function(o) {

  present <- intersect(panel, rownames(o$ms))
  X <- t(harmonise_panel(o))
  out <- matrix(0, nrow(X), length(panel), dimnames = list(NULL, panel))

  if (length(present)) {

    zz <- sweep(sweep(X[, present, drop = FALSE], 2, mu[present]), 2, sdv[present], "/")
    zz[!is.finite(zz)] <- 0
    out[, present] <- zz

  }

  out

}

fit_h <- glm(disc$y ~ ., data = data.frame(zh(disc), y = disc$y), family = binomial)
sc_h <- function(o) as.numeric(predict(fit_h, newdata = as.data.frame(zh(o)), type = "response"))

harm_rows <- lapply(list(list(nm = "ERCB external", o = ercb),
                         list(nm = "RNA-seq independent", o = rna)), function(z2) {

  s_raw <- sc5(z2$o)
  s_h <- sc_h(z2$o)
  r_raw <- roc(z2$o$y, s_raw, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  r_h <- roc(z2$o$y, s_h, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  obs <- as.integer(z2$o$y == "MN")

  data.frame(cohort = z2$nm,
             auc_raw = round(as.numeric(auc(r_raw)), 4),
             auc_harmonised = round(as.numeric(auc(r_h)), 4),
             cil_raw = round(mean(s_raw) - mean(obs), 4),
             cil_harmonised = round(mean(s_h) - mean(obs), 4),
             score_range_raw = sprintf("%.4f-%.4f", min(s_raw), max(s_raw)),
             score_range_harmonised = sprintf("%.4f-%.4f", min(s_h), max(s_h)),
             stringsAsFactors = FALSE)

})

harm <- do.call(rbind, harm_rows)

write.csv(harm, file.path(tab_dir, "09e_harmonisation_check.csv"), row.names = FALSE)

cat(sprintf("\n--- E. cross-platform harmonisation sensitivity of the locked %d-gene panel ---\n",
            length(panel)))
print(harm, row.names = FALSE)

lab_locked <- sprintf("%d-gene (locked)", length(panel))
lab_drop <- sprintf("%d-gene (without SLAMF8)", length(genes4))

roc_df <- rbind(
  data.frame(model = lab_locked, specificity = (roc(ercb$y, s5, levels = c("rest", "MN"),
                                                   direction = "<", quiet = TRUE))$specificities,
             sensitivity = (roc(ercb$y, s5, levels = c("rest", "MN"), direction = "<",
                                quiet = TRUE))$sensitivities),
  data.frame(model = lab_drop, specificity = (roc(ercb$y, s4, levels = c("rest", "MN"),
                                                   direction = "<", quiet = TRUE))$specificities,
             sensitivity = (roc(ercb$y, s4, levels = c("rest", "MN"), direction = "<",
                                quiet = TRUE))$sensitivities))
save_plot(ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, colour = model)) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
            geom_line(linewidth = 0.9) + coord_equal() +
            scale_colour_manual(values = setNames(c("firebrick", "seagreen"), c(lab_locked, lab_drop))) +
            labs(title = sprintf("Locked %d-gene vs %d-gene sensitivity model in ERCB",
                                 length(panel), length(genes4)),
                 subtitle = sprintf("AUC %.3f vs %.3f", as.numeric(auc(roc(ercb$y, s5,
                                    levels = c("rest", "MN"), direction = "<", quiet = TRUE))),
                                    obs4),
                 x = "1 - specificity", y = "sensitivity", colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = c(0.62, 0.16),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09e_no_slamf8_roc"), width = 6.2, height = 5.4)

cat("\n09e done\n")
