#!/usr/bin/env Rscript
# Paper 1 - Phase B step 8: random-effects meta-analysis of the per-cohort limma results.
#
# Design (fixed after the Phase 1 tiering):
#   PRIMARY meta  = discovery + cross-platform replication cohorts of the same compartment.
#                   Glom: GSE200828 (array) + GSE197307 (RNA-seq)
#                   TI  : GSE133288/200818 (array, ComBat) + GSE182380 (RNA-seq)
#   VALIDATION    = the ERCB cohort of the same compartment (0 % patient overlap with discovery)
#                   is held OUT of the primary meta and used only to test whether each meta-DEG
#                   replicates independently.
#   SENSITIVITY   = a meta that additionally includes the ERCB cohort, reported separately.
#
# Per-cohort effect sizes come from 04b (logFC and its limma SE). Genes absent from a cohort are
# simply not pooled for it (k is reported per gene). Random effects: metafor::rma, REML.
#
# Outputs
#   output/tables/meta/meta_<compartment>_<contrast>.csv
#   output/tables/meta/meta_<compartment>_<contrast>_metaDEG.csv
#   output/tables/meta_validation.csv
#   output/tables/meta_sensitivity_allcohorts.csv
#   output/deg/04c_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(metafor); library(ggplot2); library(ggrepel); library(pheatmap) })

deg_dir <- file.path(OUT, "tables", "deg")
meta_dir <- file.path(OUT, "tables", "meta")
fig_dir <- file.path(OUT, "deg")
tab_dir <- file.path(OUT, "tables")

dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

LFC <- 1
FDR <- 0.05
NCPU <- max(1, min(8, parallel::detectCores() - 1))

fdr_line <- function(mm) {
  p <- mm$pval[mm$p.adj < FDR]
  if (!length(p) || !any(is.finite(p))) return(1)
  max(p, na.rm = TRUE)
}

save_plot <- function(p, stem, width = 7.5, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

# base-graphics / grid figures go through explicit device blocks; fn is a closure so the objects it
# draws with stay in the calling environment
save_device_plot <- function(fn, stem, width = 7.5, height = 6) {
  while (!is.null(dev.list())) dev.off()
  dev.list()

  pdf(paste0(stem, ".pdf"), width = width, height = height)
  fn()
  dev.off()

  jpeg(paste0(stem, ".jpeg"), width = width, height = height, units = "in", res = 300,
       quality = 95, bg = "white")
  fn()
  dev.off()
}

deg_sum <- read.csv(file.path(tab_dir, "deg_summary.csv"), stringsAsFactors = FALSE)

# cohort tiers per compartment: primary meta vs held-out validation
tiers <- data.frame(
  compartment = c("Glom", "Glom", "Glom", "TI", "TI", "TI"),
  set = c("disc_glom_array", "cross_glom_rnaseq", "val_glom_array",
          "disc_ti_array", "cross_ti_rnaseq", "val_ti_array"),
  tier = c("primary", "primary", "validation", "primary", "primary", "validation"),
  stringsAsFactors = FALSE)

contrasts <- c("MN_vs_Control", "MN_vs_FSGS", "MN_vs_MCD", "MN_vs_Other", "MN_vs_nonMN")

# --- load every per-cohort table once ----------------------------------------
tabs <- list()

for (s in tiers$set) {
  for (ct in contrasts) {
    f <- file.path(deg_dir, paste0(s, "__", ct, ".csv"))
    if (file.exists(f)) {
      d <- read.csv(f, stringsAsFactors = FALSE)
      tabs[[paste(s, ct)]] <- d[, c("gene", "logFC", "SE", "P.Value", "adj.P.Val")]
    }
  }
}

cat("Loaded", length(tabs), "cohort x contrast tables\n")

# --- one random-effects meta per gene ----------------------------------------

meta_gene <- function(g, dat) {
  ok <- is.finite(dat$logFC) & is.finite(dat$SE) & dat$SE > 0
  if (!any(ok)) return(NULL)
  if (sum(ok) == 1) {
    d <- dat[ok, ]
    z <- d$logFC / d$SE
    return(data.frame(k = 1L, logFC = d$logFC, SE = d$SE, zval = z,
                      pval = 2 * pnorm(-abs(z)), ci.lb = d$logFC - 1.96 * d$SE,
                      ci.ub = d$logFC + 1.96 * d$SE,
                      QE = NA_real_, QEp = NA_real_, I2 = NA_real_, tau2 = NA_real_,
                      used_dl = FALSE))
  }
  # REML Fisher scoring does not converge for a small number of genes (near-identical SEs with
  # zero heterogeneity); fall back to DerSimonian-Laird rather than losing the gene from the pool
  # REML Fisher scoring can fail to converge. metafor 5.2.1 raises that as a WARNING (not an
  # error) and still returns a fit, so `tryCatch(error=)` alone never fired; the warning is now
  # captured, the gene is refitted with DerSimonian-Laird, and the switch is recorded per gene.
  used_dl <- FALSE
  warned <- FALSE
  fit <- tryCatch(
    withCallingHandlers(
      rma(yi = dat$logFC[ok], sei = dat$SE[ok], method = "REML", slab = rownames(dat)[ok]),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) warned <<- TRUE
        invokeRestart("muffleWarning")
      }),
    error = function(e) NULL)

  if (is.null(fit) || warned) {
    fit <- rma(yi = dat$logFC[ok], sei = dat$SE[ok], method = "DL", slab = rownames(dat)[ok])
    used_dl <- TRUE
  }

  data.frame(k = sum(ok), logFC = as.numeric(fit$beta), SE = fit$se, zval = fit$zval,
             pval = fit$pval, ci.lb = fit$ci.lb, ci.ub = fit$ci.ub,
             QE = fit$QE, QEp = fit$QEp, I2 = fit$I2, tau2 = fit$tau2,
             used_dl = used_dl)
}

run_meta <- function(comp, ct, sets) {
  sets <- sets[vapply(sets, function(s) !is.null(tabs[[paste(s, ct)]]), logical(1))]
  if (!length(sets)) return(NULL)
  mats <- lapply(sets, function(s) {
    d <- tabs[[paste(s, ct)]]
    rownames(d) <- d$gene
    d
  })
  names(mats) <- sets
  genes <- Reduce(intersect, lapply(mats, rownames))
  dat_all <- do.call(cbind, lapply(mats, function(d) d[genes, c("logFC", "SE")]))
  res <- parallel::mclapply(seq_along(genes), function(j) {
    ok_fit <- TRUE
    dm <- t(vapply(mats, function(d) unlist(d[genes[j], c("logFC", "SE")], use.names = FALSE),
                   numeric(2)))
    rownames(dm) <- sets
    colnames(dm) <- c("logFC", "SE")
    out <- tryCatch(meta_gene(genes[j], as.data.frame(dm)),
                    error = function(e) NULL)
    if (is.null(out)) return(NULL)
    cbind(gene = genes[j], out)
  }, mc.cores = NCPU)
  # an mclapply worker that throws returns a try-error object, not NULL - it must never be
  # recycled into the result table (that silently produced error-message "genes" once)
  res <- res[vapply(res, function(z) is.data.frame(z) && nrow(z) == 1, logical(1))]
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res)
  for (v in c("k", "logFC", "SE", "zval", "pval", "ci.lb", "ci.ub", "QE", "QEp", "I2", "tau2",
              "used_dl")) {
    out[[v]] <- as.numeric(out[[v]])
  }
  # per-cohort logFC columns, so every pooled estimate can be traced back to its cohorts
  for (s in sets) {
    out[[paste0("logFC_", s)]] <- mats[[s]][out$gene, "logFC"]
    out[[paste0("adjP_", s)]] <- mats[[s]][out$gene, "adj.P.Val"]
  }
  ncol_lfc <- paste0("logFC_", sets)
  lfc_mat <- as.matrix(out[, ncol_lfc, drop = FALSE])
  out$consistent_direction <- apply(sign(lfc_mat), 1, function(z) all(z == z[1], na.rm = TRUE))
  out$p.adj <- p.adjust(out$pval, method = "BH")
  out$n_not_estimable <- sum(!is.finite(out$pval))
  out$estimable <- is.finite(out$pval) & is.finite(out$logFC) & is.finite(out$p.adj) &
    !is.na(out$consistent_direction)
  out$meta_deg <- out$estimable & abs(out$logFC) > LFC & out$p.adj < FDR &
    out$consistent_direction
  out$i2_high <- !is.na(out$I2) & out$I2 > 75
  out$cohorts <- paste(sets, collapse = "+")
  # a k = 1 "meta" is a relabelled per-cohort estimate, not a pooled one; say so in the table
  out$pooling <- if (max(out$k, na.rm = TRUE) == 1) "single cohort (no pooling)" else
    "random effects (REML)"
  out <- out[order(out$pval), ]
  cbind(compartment = comp, contrast = ct, out)
}

meta_rows <- list()
val_rows <- list()
sens_rows <- list()

for (comp in c("Glom", "TI")) {

  primary_sets <- tiers$set[tiers$compartment == comp & tiers$tier == "primary"]
  val_sets <- tiers$set[tiers$compartment == comp & tiers$tier == "validation"]

  for (ct in contrasts) {

    cat("\n--", comp, ct, "| primary:", paste(primary_sets, collapse = ", "), "\n")

    m <- run_meta(comp, ct, primary_sets)

    if (is.null(m)) { cat("   nothing to pool\n"); next }

    f <- file.path(meta_dir, paste0("meta_", comp, "_", ct, ".csv"))
    write.csv(m, f, row.names = FALSE)

    sig <- m[m$meta_deg, ]
    write.csv(sig, file.path(meta_dir, paste0("meta_", comp, "_", ct, "_metaDEG.csv")),
              row.names = FALSE)

    meta_rows[[paste(comp, ct)]] <- data.frame(
      compartment = comp, contrast = ct, cohorts = m$cohorts[1],
      k_max = max(m$k, na.rm = TRUE), pooling = m$pooling[1],
      genes_pooled = nrow(m), meta_deg_up = sum(m$meta_deg & m$logFC > 0),
      meta_deg_down = sum(m$meta_deg & m$logFC < 0), meta_deg_total = sum(m$meta_deg),
      meta_deg_i2_high = sum(m$meta_deg & m$i2_high),
      n_dl_fallback = sum(m$used_dl, na.rm = TRUE),
      median_i2_meta_deg = if (all(is.na(m$I2[m$meta_deg]))) NA_real_ else
        round(median(m$I2[m$meta_deg], na.rm = TRUE), 1),
      median_tau2_meta_deg = if (all(is.na(m$tau2[m$meta_deg]))) NA_real_ else
        round(median(m$tau2[m$meta_deg], na.rm = TRUE), 4),
      stringsAsFactors = FALSE)

    cat("   pooled", nrow(m), "genes | meta-DEG", sum(m$meta_deg),
        "( up", sum(m$meta_deg & m$logFC > 0), "/ down", sum(m$meta_deg & m$logFC < 0),
        ") | k_max", max(m$k), "\n")

    # --- independent validation against the ERCB cohort ----------------------
    vs <- val_sets[vapply(val_sets, function(s) !is.null(tabs[[paste(s, ct)]]), logical(1))]

    for (s in vs) {

      v <- tabs[[paste(s, ct)]]
      common <- intersect(m$gene, v$gene)
      mm <- m[match(common, m$gene), ]
      vv <- v[match(common, v$gene), ]

      deg <- mm$meta_deg

      val_rows[[paste(comp, ct, s)]] <- data.frame(
        compartment = comp, contrast = ct, validation_set = s,
        genes_in_common = length(common),
        meta_degs_tested = sum(deg),
        replicated_adjP05_samesign = sum(deg & vv$adj.P.Val < FDR & sign(vv$logFC) == sign(mm$logFC)),
        replicated_adjP05_samesign_pct = round(100 * sum(deg & vv$adj.P.Val < FDR &
                                                          sign(vv$logFC) == sign(mm$logFC)) /
                                                 max(1, sum(deg)), 1),
        sign_concordance_all_genes_pct = round(100 * mean(sign(vv$logFC) == sign(mm$logFC)), 1),
        lfc_pearson_all_genes = round(cor(mm$logFC, vv$logFC), 3),
        stringsAsFactors = FALSE)
    }

    # --- sensitivity: add the ERCB cohort to the meta ------------------------
    m2 <- run_meta(comp, ct, c(primary_sets, vs))

    if (!is.null(m2) && !is.null(vs)) {

      write.csv(m2, file.path(meta_dir, paste0("meta_allcohorts_", comp, "_", ct, ".csv")),
                row.names = FALSE)

      # how many primary meta-DEGs stay significant once ERCB is pooled in?
      keep <- m$gene[m$meta_deg]

      sens_rows[[paste(comp, ct)]] <- data.frame(
        compartment = comp, contrast = ct, cohorts = m2$cohorts[1],
        k_max = max(m2$k, na.rm = TRUE),
        primary_meta_degs = length(keep),
        still_significant_allcohorts = sum(m2$meta_deg[match(keep, m2$gene)], na.rm = TRUE),
        allcohorts_meta_degs = sum(m2$meta_deg),
        stringsAsFactors = FALSE)
    }

    # --- volcano -------------------------------------------------------------
    mm <- m
    mm$status <- ifelse(mm$meta_deg, ifelse(mm$logFC > 0, "up in MN", "down in MN"), "ns")
    mm$status <- factor(mm$status, levels = c("up in MN", "down in MN", "ns"))

    lab <- mm[mm$meta_deg, ]
    lab <- lab[order(-abs(lab$logFC)), ]
    lab <- head(lab, 12)

    p <- ggplot(mm, aes(logFC, -log10(pval), colour = status)) +
      geom_point(size = 1.3, alpha = 0.7) +
      geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey40") +
      geom_hline(yintercept = -log10(fdr_line(mm)), linetype = "dashed", colour = "grey40") +
      geom_text_repel(data = lab, aes(label = gene), size = 2.4, max.overlaps = 25,
                      show.legend = FALSE) +
      scale_colour_manual(values = c("up in MN" = "#C0392B", "down in MN" = "#2471A3",
                                     "ns" = "grey65")) +
      labs(title = paste0(comp, " - meta-analysis ", gsub("_", " ", ct),
                          " (", paste(primary_sets, collapse = " + "), ")"),
           subtitle = paste0("dashed horizontal line = nominal p corresponding to BH < 0.05 ",
                             "(significance is judged on BH-adjusted p, not on this line)"),
           x = "random-effects pooled log2FC", y = "-log10 pooled P") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.title = element_blank())

    save_plot(p, file.path(fig_dir, paste0("04c_", comp, "_", ct, "_volcano")))

    # --- forest plots for the strongest meta-DEGs ----------------------------
    forest_genes <- head(sig$gene[order(sig$pval)], 5)

    if (length(forest_genes)) {

      sets_used <- strsplit(m$cohorts[1], "\\+")[[1]]

      for (g in forest_genes) {

        dm <- do.call(rbind, lapply(sets_used, function(s) {
          d <- tabs[[paste(s, ct)]]
          if (is.null(d) || !(g %in% d$gene)) return(NULL)
          r <- d[d$gene == g, ]
          data.frame(cohort = s, logFC = r$logFC, SE = r$SE)
        }))

        fit <- tryCatch(rma(yi = dm$logFC, sei = dm$SE, method = "REML", slab = dm$cohort),
                        error = function(e) NULL)
        if (is.null(fit)) {
          fit <- tryCatch(rma(yi = dm$logFC, sei = dm$SE, method = "DL", slab = dm$cohort),
                          error = function(e) NULL)
        }
        if (is.null(fit)) next

        save_device_plot(function() forest(fit,
                                           main = paste0(g, " - ", comp, " ", gsub("_", " ", ct)),
                                           xlab = "log2 fold change (MN vs comparator)",
                                           col = "#C0392B"),
                         file.path(fig_dir, paste0("04c_", comp, "_", ct, "_forest_", g)),
                         width = 7, height = 3.4)
      }
    }
  }
}

meta_sum <- do.call(rbind, meta_rows)
val_sum <- do.call(rbind, val_rows)
sens_sum <- if (length(sens_rows)) do.call(rbind, sens_rows) else NULL

write.csv(meta_sum, file.path(tab_dir, "meta_summary.csv"), row.names = FALSE)
write.csv(val_sum, file.path(tab_dir, "meta_validation.csv"), row.names = FALSE)

if (!is.null(sens_sum)) {
  write.csv(sens_sum, file.path(tab_dir, "meta_sensitivity_allcohorts.csv"), row.names = FALSE)
}

# --- cohort-level logFC heatmap of the strongest meta-DEGs -------------------

for (comp in c("Glom", "TI")) {

  f <- file.path(meta_dir, paste0("meta_", comp, "_MN_vs_nonMN.csv"))

  if (!file.exists(f)) next

  m <- read.csv(f, stringsAsFactors = FALSE)
  sig <- head(m[m$meta_deg, ][order(m$pval[m$meta_deg]), ], 40)

  if (nrow(sig) < 3) next

  sets_used <- strsplit(m$cohorts[1], "\\+")[[1]]
  lfc_cols <- paste0("logFC_", sets_used)
  hm <- as.matrix(sig[, lfc_cols, drop = FALSE])
  rownames(hm) <- sig$gene
  colnames(hm) <- sets_used

  lim <- max(abs(hm), na.rm = TRUE)

  save_device_plot(function() pheatmap(hm, cluster_cols = FALSE, cluster_rows = TRUE,
                                       color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(101),
                                       breaks = seq(-lim, lim, length.out = 102),
                                       fontsize_row = 7, angle_col = 45,
                                       main = paste0(comp, " top 40 meta-DEGs (MN vs non-MN): per-cohort log2FC")),
                   file.path(fig_dir, paste0("04c_", comp, "_metaDEG_heatmap")),
                   width = 6, height = 8)
}

cat("\n================ META SUMMARY ================\n")

print(meta_sum, row.names = FALSE)

cat("\n---- independent validation (ERCB) ----\n")

print(val_sum, row.names = FALSE)

cat("\n---- sensitivity: ERCB pooled into the meta ----\n")

if (!is.null(sens_sum)) print(sens_sum, row.names = FALSE)

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 04c_meta_analysis.R"),
             capture.output(print(meta_sum, row.names = FALSE)),
             capture.output(print(val_sum, row.names = FALSE))),
           file.path(PROJ, "logs", "meta_analysis.log"))

cat("\nDone.\n")
