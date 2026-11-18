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
#   D. FORMAL 4-GENE MODEL - the leave-one-gene-out in 09d showed that dropping SLAMF8 (detected in
#      only 9.5% of RNA-seq samples) gives 0.900 / 0.895. Here that model is reported properly:
#      external AUC with CI, a permutation null and recalibration, so it can be footnoted honestly.
#
# Outputs
#   output/tables/09e_panel_by_comparator.csv
#   output/tables/09e_authors_matrix_check.csv
#   output/tables/09e_technical_confound.csv
#   output/tables/09e_four_gene_model.csv
#   output/models_panel/09e_{by_comparator,authors_check,four_gene_roc}.{pdf,jpeg}

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
save_plot(ggplot(bc, aes(x = label, y = auc, fill = cohort)) +
            geom_col(width = 0.68) +
            geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
            geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
            coord_flip(ylim = c(0.3, 1)) +
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
  save_plot(ggplot(chk_df, aes(x = reorder(gene, -auc), y = auc, fill = source)) +
              geom_col(position = position_dodge(0.72), width = 0.66) +
              geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
              coord_cartesian(ylim = c(0.4, 1)) +
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

# --- D. formal 4-gene model (drop SLAMF8) ---------------------------------
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
  data.frame(model = "4-gene (TRPC6, MYOZ2, ICAM1, PPP1R3C)", cohort = z$nm,
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
write.csv(four, file.path(tab_dir, "09e_four_gene_model.csv"), row.names = FALSE)
cat("\n--- formal 4-gene model ---\n"); print(four, row.names = FALSE)

roc_df <- rbind(
  data.frame(model = "5-gene (locked)", specificity = (roc(ercb$y, s5, levels = c("rest", "MN"),
                                                           direction = "<", quiet = TRUE))$specificities,
             sensitivity = (roc(ercb$y, s5, levels = c("rest", "MN"), direction = "<",
                                quiet = TRUE))$sensitivities),
  data.frame(model = "4-gene (without SLAMF8)", specificity = (roc(ercb$y, s4, levels = c("rest", "MN"),
                                                                   direction = "<", quiet = TRUE))$specificities,
             sensitivity = (roc(ercb$y, s4, levels = c("rest", "MN"), direction = "<",
                                quiet = TRUE))$sensitivities))
save_plot(ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, colour = model)) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
            geom_line(linewidth = 0.9) + coord_equal() +
            scale_colour_manual(values = c("5-gene (locked)" = "firebrick",
                                           "4-gene (without SLAMF8)" = "seagreen")) +
            labs(title = "Locked 5-gene vs 4-gene sensitivity model in ERCB",
                 subtitle = sprintf("AUC %.3f vs %.3f", as.numeric(auc(roc(ercb$y, s5,
                                    levels = c("rest", "MN"), direction = "<", quiet = TRUE))),
                                    obs4),
                 x = "1 - specificity", y = "sensitivity", colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = c(0.62, 0.16),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09e_four_gene_roc"), width = 6.2, height = 5.4)

cat("\n09e done\n")
