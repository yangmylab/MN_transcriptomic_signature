#!/usr/bin/env Rscript
# Paper 1 - Phase E step 4: where do the five panel genes live, and does a module-level model explain
# the panel's performance?
#
# Three questions the existing data can answer without new experiments:
#   1. DIAGNOSTICS - is every panel gene actually expressed in kidney tissue? A gene that is absent
#      from the single-cell atlas but carries a large bulk coefficient (MYOZ2 is the case here)
#      needs to be flagged rather than quietly shipped in a clinical panel.
#   2. LOCALISATION - which cell type expresses each panel gene in the kidney atlas, and does the gene
#      move with MN in that cell type's pseudobulk? This is the bridge to the IF validation and the
#      reason the panel is interpretable.
#   3. MODULE-LEVEL MODEL - do three module scores (pink, black, red) match the locked gene panel?
#      If yes, the panel is a compact readout of two/three named programmes; if no, it carries
#      information the modules do not.
#
# Outputs
#   output/tables/09d_panel_expression_diagnostics.csv   detection + atlas expression per gene/cohort
#   output/tables/09d_panel_pseudobulk.csv               MN-vs-control logFC per cell type
#   output/tables/09d_module_vs_panel.csv                AUC: module scores vs genes vs combined
#   output/models_panel/09d_{gene_celltype,pseudobulk,module_vs_panel}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(pROC); library(ggplot2) })

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
sm <- read.csv(file.path(tab_dir, "signature_matrix.csv"), row.names = 1, check.names = FALSE)
mods <- read.csv(file.path(tab_dir, "wgcna_modules_glom.csv"), stringsAsFactors = FALSE)

load_set <- function(set) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  keep <- !(grp %in% c("Control", "UNKNOWN"))
  list(ms = ms[, keep, drop = FALSE],
       y = factor(ifelse(grp[keep] == "MN", "MN", "rest"), levels = c("rest", "MN")))
}

sets <- c("disc_glom_array", "val_glom_array", "cross_glom_rnaseq", "pla2r_lcm")
obj <- lapply(sets, load_set); names(obj) <- sets

# --- 1. expression diagnostics --------------------------------------------
diag_rows <- list()
for (s in sets) {
  ms <- obj[[s]]$ms
  for (g in panel) {
    if (!g %in% rownames(ms)) {
      diag_rows[[paste(s, g)]] <- data.frame(cohort = s, gene = g, detected_fraction = NA,
                                             median_expr = NA, max_expr = NA,
                                             note = "gene absent from this platform",
                                             stringsAsFactors = FALSE)
      next
    }
    v <- ms[g, ]
    det <- mean(!is.na(v) & v > 5)
    diag_rows[[paste(s, g)]] <- data.frame(
      cohort = s, gene = g, detected_fraction = round(det, 4),
      median_expr = round(median(v, na.rm = TRUE), 3), max_expr = round(max(v, na.rm = TRUE), 3),
      note = if (det < 0.5) "low detection" else "", stringsAsFactors = FALSE)
  }
}
diag <- do.call(rbind, diag_rows)
# atlas expression: highest CPM and the cell type that carries it
atlas <- data.frame(gene = panel, stringsAsFactors = FALSE)
atlas$atlas_max_cpm <- round(apply(sm[panel, , drop = FALSE], 1, max), 2)
atlas$atlas_top_celltype <- colnames(sm)[apply(sm[panel, , drop = FALSE], 1, which.max)]
atlas$atlas_frac_zero <- round(rowMeans(sm[panel, , drop = FALSE] == 0), 3)
# name it atlas_note: a second column called "note" would be renamed note.x/note.y by merge()
atlas$atlas_note <- ifelse(atlas$atlas_max_cpm < 1, "essentially absent from the kidney atlas",
                    ifelse(atlas$atlas_max_cpm < 10, "weakly expressed", ""))
write.csv(merge(diag, atlas, by = "gene"), file.path(tab_dir, "09d_panel_expression_diagnostics.csv"),
          row.names = FALSE)
cat("--- expression diagnostics ---\n"); print(atlas, row.names = FALSE)
print(diag[diag$note != "" | is.na(diag$detected_fraction), ], row.names = FALSE)

# --- 2. pseudobulk localisation -------------------------------------------
pb <- read.csv(file.path(tab_dir, "scrna_pseudobulk_de.csv"), stringsAsFactors = FALSE)
pb <- pb[pb$gene %in% panel & pb$contrast == "MN_vs_control", ]
pb <- pb[, c("celltype", "gene", "n_case", "n_control", "logFC", "P.Value", "adj.P.Val")]
pb$logFC <- round(pb$logFC, 3)
write.csv(pb, file.path(tab_dir, "09d_panel_pseudobulk.csv"), row.names = FALSE)

if (nrow(pb)) {
  ord <- names(sort(tapply(abs(pb$logFC), pb$gene, max, na.rm = TRUE), decreasing = TRUE))
  pb$gene <- factor(pb$gene, levels = ord)
  save_plot(ggplot(pb, aes(x = celltype, y = gene, fill = logFC)) +
              geom_tile(colour = "white") +
              geom_text(aes(label = ifelse(is.na(logFC), "", sprintf("%.2f", logFC))), size = 2.3) +
              scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                   midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
              labs(title = "Panel genes in MN vs control pseudobulk, per cell type",
                   subtitle = "limma-voom on donor-level pseudobulk; **no** comparison survives BH (n = 15 vs 10 donors)",
                   x = NULL, y = NULL, fill = "log2FC") +
              theme_bw(base_size = 9) +
              theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_text(face = "bold")),
            file.path(fig_dir, "09d_pseudobulk"), width = 7.4, height = 4)
}

# --- 3. module scores vs the gene panel -----------------------------------
module_score <- function(ms, genes, train = FALSE, ref = NULL) {
  g <- intersect(genes, rownames(ms))
  x <- t(scale(t(ms[g, , drop = FALSE])))
  x[!is.finite(x)] <- 0
  pc <- prcomp(t(x), center = FALSE, scale. = FALSE)
  s <- pc$x[, 1]
  if (train) return(list(score = s, loadings = pc$rotation[, 1]))
  # orient the validation score with the training loadings so a flipped PC cannot fake an AUC
  as.numeric(t(x) %*% ref)
}

mod_sets <- list(
  pink = mods$gene[mods$module == "pink"],
  black = mods$gene[mods$module == "black"],
  red = mods$gene[mods$module == "red"])

feat_train <- sapply(names(mod_sets), function(m) {
  r <- module_score(obj$disc_glom_array$ms, mod_sets[[m]], train = TRUE)
  attr(r$score, "load") <- r$loadings
  r$score
})
loadings <- lapply(names(mod_sets), function(m)
  module_score(obj$disc_glom_array$ms, mod_sets[[m]], train = TRUE)$loadings)
names(loadings) <- names(mod_sets)

score_cohort <- function(o) {
  sapply(names(mod_sets), function(m) {
    g <- intersect(mod_sets[[m]], rownames(o$ms))
    x <- t(scale(t(o$ms[g, , drop = FALSE]))); x[!is.finite(x)] <- 0
    as.numeric(t(x) %*% loadings[[m]][g])
  })
}

res <- list()
for (s in c("disc_glom_array", "val_glom_array", "cross_glom_rnaseq")) {
  o <- obj[[s]]
  f <- score_cohort(o)
  d <- data.frame(y = o$y, f)
  fit <- glm(y ~ ., data = data.frame(score_cohort(obj$disc_glom_array), y = obj$disc_glom_array$y),
             family = binomial)
  s_mod <- as.numeric(predict(fit, newdata = as.data.frame(f), type = "response"))
  r <- roc(o$y, s_mod, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  res[[s]] <- data.frame(cohort = s, model = "3 module scores (pink+black+red)",
                         auc = round(as.numeric(auc(r)), 4),
                         ci_low = round(as.numeric(ci.auc(r)[1]), 4),
                         ci_high = round(as.numeric(ci.auc(r)[3]), 4), stringsAsFactors = FALSE)
}

ext <- read.csv(file.path(tab_dir, "09b_panel_external.csv"), stringsAsFactors = FALSE)
gene_rows <- data.frame(
  cohort = c("disc_glom_array", "val_glom_array", "cross_glom_rnaseq"),
  # label derived from the panel actually read from 09a, not hard-coded: the audit re-run
  # changed the consensus panel size (5 -> 4 genes), and a hard-coded "5-gene" label then
  # mis-describes a 4-gene model in both the table and the figure.
  model = sprintf("%d-gene panel (locked)", length(panel)),
  auc = ext$auc[match(c("Discovery (apparent, not independent)", "ERCB external",
                        "RNA-seq independent cohort"), ext$cohort)],
  ci_low = ext$ci_low[match(c("Discovery (apparent, not independent)", "ERCB external",
                              "RNA-seq independent cohort"), ext$cohort)],
  ci_high = ext$ci_high[match(c("Discovery (apparent, not independent)", "ERCB external",
                                "RNA-seq independent cohort"), ext$cohort)],
  stringsAsFactors = FALSE)

mvp <- rbind(gene_rows, do.call(rbind, res))
write.csv(mvp, file.path(tab_dir, "09d_module_vs_panel.csv"), row.names = FALSE)
cat("\n--- module scores vs gene panel ---\n"); print(mvp[, c("cohort", "model", "auc", "ci_low", "ci_high")],
                                                      row.names = FALSE)

mvp$label <- paste0(mvp$cohort, "\n",
                    ifelse(grepl("-gene panel", mvp$model), sprintf("%d genes", length(panel)), "3 modules"))
save_plot(ggplot(mvp, aes(x = label, y = auc, fill = model)) +
            geom_col(position = position_dodge(0.7), width = 0.62) +
            geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(0.7),
                          width = 0.18) +
            geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
            coord_cartesian(ylim = c(0.4, 1)) +
            scale_fill_brewer(palette = "Dark2") +
            labs(title = sprintf("Does a three-module score match the %d-gene panel?", length(panel)),
                 x = NULL, y = "AUC (95% CI)", fill = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = "bottom",
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09d_module_vs_panel"), width = 8, height = 4.6)

# --- 4. leave-one-gene-out sensitivity -----------------------------------
# obj$...$ms is GENES x SAMPLES, so gene-wise summaries are row-wise (colMeans/apply(...,2) would
# give one value per sample and the sweep below would silently produce all-NA columns)
mu  <- rowMeans(obj$disc_glom_array$ms[panel, , drop = FALSE])
sdv <- apply(obj$disc_glom_array$ms[panel, , drop = FALSE], 1, sd)
sdv[!is.finite(sdv) | sdv == 0] <- 1

zscore <- function(ms, genes) {
  out <- matrix(0, ncol(ms), length(genes), dimnames = list(NULL, genes))
  have <- intersect(genes, rownames(ms))
  if (length(have)) {
    zz <- sweep(sweep(t(ms[have, , drop = FALSE]), 2, mu[have]), 2, sdv[have], "/")
    zz[!is.finite(zz)] <- 0
    out[, have] <- zz
  }
  out
}

fit_and_auc <- function(train_obj, test_objs, genes) {
  fit <- glm(train_obj$y ~ ., data = data.frame(zscore(train_obj$ms, genes), y = train_obj$y),
             family = binomial)
  sapply(test_objs, function(o) {
    s <- as.numeric(predict(fit, newdata = as.data.frame(zscore(o$ms, genes)), type = "response"))
    as.numeric(auc(roc(o$y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
  })
}

test_objs <- list(ERCB = obj$val_glom_array, RNAseq = obj$cross_glom_rnaseq)
loo <- do.call(rbind, lapply(c("full", panel), function(g) {
  genes <- if (g == "full") panel else setdiff(panel, g)
  a <- fit_and_auc(obj$disc_glom_array, test_objs, genes)
  data.frame(model = if (g == "full") sprintf("all %d genes", length(panel)) else paste0("without ", g),
             n_genes = length(genes),
             ERCB_auc = round(a[["ERCB"]], 4), RNAseq_auc = round(a[["RNAseq"]], 4),
             stringsAsFactors = FALSE)
}))
single <- do.call(rbind, lapply(panel, function(g) {
  a <- fit_and_auc(obj$disc_glom_array, test_objs, g)
  data.frame(model = paste0("only ", g), n_genes = 1,
             ERCB_auc = round(a[["ERCB"]], 4), RNAseq_auc = round(a[["RNAseq"]], 4),
             stringsAsFactors = FALSE)
}))
loo <- rbind(loo, single)
write.csv(loo, file.path(tab_dir, "09d_leave_one_out.csv"), row.names = FALSE)
cat("\n--- leave-one-gene-out ---\n"); print(loo, row.names = FALSE)

loo_df <- rbind(
  data.frame(model = loo$model, cohort = "ERCB", auc = loo$ERCB_auc),
  data.frame(model = loo$model, cohort = "RNA-seq", auc = loo$RNAseq_auc))
loo_df$model <- factor(loo_df$model, levels = rev(loo$model))
save_plot(ggplot(loo_df, aes(x = model, y = auc, fill = cohort)) +
            geom_col(position = position_dodge(0.7), width = 0.65) +
            geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
            coord_flip(ylim = c(0.4, 1)) +
            scale_fill_brewer(palette = "Set2") +
            labs(title = "Leave-one-gene-out: does any single gene carry the panel?",
                 subtitle = "models refit on discovery; evaluated in ERCB and the RNA-seq cohort",
                 x = NULL, y = "AUC", fill = NULL) +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), legend.position = "bottom",
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09d_leave_one_out"), width = 7.6, height = 5.4)

# gene-by-cell-type atlas heatmap
hm <- sm[panel, , drop = FALSE]
hm_df <- data.frame(gene = rep(rownames(hm), times = ncol(hm)),
                    celltype = rep(colnames(hm), each = nrow(hm)),
                    cpm = as.vector(as.matrix(hm)))
hm_df$gene <- factor(hm_df$gene, levels = rev(panel))
save_plot(ggplot(hm_df, aes(x = celltype, y = gene, fill = log1p(cpm))) +
            geom_tile(colour = "white") +
            scale_fill_gradient(low = "white", high = "firebrick") +
            labs(title = "Panel genes across the kidney atlas (CPM, log1p)",
                 subtitle = "signature matrix from 07e; pDC and Neutrophil columns are unreliable (B8)",
                 x = NULL, y = NULL, fill = "log1p CPM") +
            theme_bw(base_size = 9) +
            theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09d_gene_celltype"), width = 7.6, height = 3.6)

cat("\n09d done\n")
