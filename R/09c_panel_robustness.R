#!/usr/bin/env Rscript
# Paper 1 - Phase E step 3: does the locked panel survive the questions a reviewer will ask?
#
# 09b established discrimination (ERCB 0.881, RNA-seq 0.894) and flagged cohort-dependent
# miscalibration. This step tests five things the current data can actually answer:
#
#   A. RECALIBRATION  - split ERCB in half, fit an intercept+slope recalibration on half A, score
#                       half B before/after. Shows miscalibration is fixable without fooling
#                       yourself with the whole validation set.
#   B. THIRD COHORT   - pla2r_lcm (GSE216841, LCM glomeruli, 12 MN vs 14 MCD), a sample source no
#                       other arm shares. Small, so it is reported as an independent but
#                       underpowered check.
#   C. COMPARTMENT    - apply the glomerular panel to the tubulointerstitial cohorts. A glomerular
#                       signature should NOT transfer into the interstitium if the compartment
#                       claim in the manuscript is real.
#   D. ANNOTATION     - GSE99339 is the same ERCB cohort processed with a different Brainarray CDF
#                       revision. Same biology, different probe annotation: if the panel survives
#                       that, it is not a probe-set artefact.
#   E. CROSS-PLATFORM TRAINING - fit the five genes on the RNA-seq cohort (z-scored within cohort)
#                       and test on ERCB. Tests whether the coefficients themselves transfer.
#
# Outputs
#   output/tables/09c_panel_robustness.csv        per cohort: endpoint, n, AUC, CI, provenance
#   output/tables/09c_panel_recalibration.csv     before/after calibration on the held-out half
#   output/tables/09c_crossplatform_coefficients.csv
#   output/models_panel/09c_{robustness_roc,recalibration}.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(pROC); library(ggplot2); library(org.Hs.eg.db) })

set_dir <- file.path(PROJ, "work", "sets")
rma_dir <- file.path(PROJ, "work", "rma")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_panel")
data_dir <- DATA

SEED <- 20261019L
set.seed(SEED)

save_plot <- function(p, stem, width = 8, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

panel <- read.csv(file.path(tab_dir, "09a_panel_consensus.csv"), stringsAsFactors = FALSE)$gene
cat("locked panel:", paste(panel, collapse = ", "), "\n")

load_set <- function(set, keep_control = FALSE) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(ms), ann$gsm))
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  grp[!(grp %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  keep <- if (keep_control) grp != "UNKNOWN" else !(grp %in% c("Control", "UNKNOWN"))
  list(X = t(ms[, keep, drop = FALSE]),
       y = factor(ifelse(grp[keep] == "MN", "MN", "rest"), levels = c("rest", "MN")),
       group = factor(grp[keep], levels = c("Control", "FSGS", "MCD", "Other", "MN")))
}

# --- locked model from discovery -------------------------------------------
disc <- load_set("disc_glom_array")
mu   <- colMeans(disc$X[, panel, drop = FALSE])
sdv  <- apply(disc$X[, panel, drop = FALSE], 2, sd)
sdv[!is.finite(sdv) | sdv == 0] <- 1

# a panel gene absent from a cohort (SLAMF8 is missing on the TI RNA-seq platform) is scored at the
# discovery mean, i.e. z = 0, instead of crashing the subset
zscore <- function(X) {
  out <- matrix(0, nrow(X), length(panel), dimnames = list(NULL, panel))
  have <- intersect(panel, colnames(X))
  if (length(have)) {
    m <- X[, have, drop = FALSE]
    zz <- sweep(sweep(m, 2, mu[have]), 2, sdv[have], "/")
    zz[!is.finite(zz)] <- 0
    out[, have] <- zz
  }
  out
}

fit <- glm(disc$y ~ ., data = data.frame(zscore(disc$X), y = disc$y), family = binomial)
score <- function(X) as.numeric(predict(fit, newdata = as.data.frame(zscore(X)), type = "response"))

eval_auc <- function(obj, label, note) {
  s <- score(obj$X)
  r <- roc(obj$y, s, levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  data.frame(cohort = label, endpoint = "MN vs other GN", n = length(obj$y),
             n_MN = sum(obj$y == "MN"), n_other = sum(obj$y == "rest"),
             auc = round(as.numeric(auc(r)), 4),
             ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
             note = note, stringsAsFactors = FALSE)
}

rows <- list()
rows[["ercb"]]  <- eval_auc(load_set("val_glom_array"), "ERCB glomerular (reference)",
                            "same result as 09b, repeated for the table")
rows[["ti"]]    <- eval_auc(load_set("val_ti_array"), "ERCB tubulointerstitial",
                            "compartment specificity: glomerular panel applied to TI")
rows[["rnati"]] <- eval_auc(load_set("cross_ti_rnaseq"), "RNA-seq tubulointerstitial",
                            "SLAMF8 absent on this platform -> scored at the discovery mean")
rows[["lcm"]]   <- eval_auc(load_set("pla2r_lcm"), "LCM glomeruli GSE216841 (third cohort)",
                            "12 MN vs 14 MCD; PLA2R status NOT in the metadata (see B2/P-45)")

# --- D. annotation-revision robustness (GSE99339, same ERCB cohort, other CDF) ---
read_sm_titles <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE)
  ttl <- ln[startsWith(ln, "!Sample_title")]
  gsm <- ln[startsWith(ln, "!Sample_geo_accession")]
  parse <- function(x) gsub('^"|"$', "", strsplit(x, "\t")[[1]][-1])
  data.frame(gsm = parse(gsm), title = parse(ttl), stringsAsFactors = FALSE)
}

annot_gse99339 <- local({
  # filenames are <GSE>-<GPL>_series_matrix.txt.gz, so the GPL token comes BEFORE "series_matrix"
  # the release is split over two platform files: GPL19109 (Plus2) and GPL19184 (133A). An earlier
  # pattern of "GPL1910[94]" matched only the first, which silently excluded the 21 MGN samples.
  fs <- list.files(file.path(data_dir, "GSE99339"), pattern = "GPL191(09|84)_series_matrix",
                   full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, function(f) tryCatch(read_sm_titles(f), error = function(e) NULL)))
})

robust_row <- NULL
if (!is.null(annot_gse99339)) {
  g39 <- readRDS(file.path(rma_dir, "GSE99339_gene_rma_merged.rds"))
  # Entrez-keyed: map to symbols, keep the best probe-free row per symbol
  map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = rownames(g39),
                                                keytype = "ENTREZID", columns = "SYMBOL"))
  map <- map[!is.na(map$SYMBOL) & map$SYMBOL %in% panel, ]
  if (nrow(map)) {
    ens <- sapply(split(map$ENTREZID, map$SYMBOL), function(ids) {
      ids <- intersect(ids, rownames(g39))
      if (!length(ids)) return(NA_character_)
      ids[which.max(rowMeans(g39[ids, , drop = FALSE], na.rm = TRUE))]
    })
    ens <- ens[!is.na(ens)]
    sub <- g39[unname(ens), , drop = FALSE]
    rownames(sub) <- names(ens)

    a <- annot_gse99339
    a <- a[a$gsm %in% colnames(sub), ]
    # GEO titles carry the diagnosis; there is no participant id, so this is a cohort-level check
    lab <- ifelse(grepl("-MGN[0-9]+|-MN[0-9]+|MN$", a$title), "MN",
           ifelse(grepl("-LD[0-9]+|Living", a$title, ignore.case = TRUE), "Control",
           ifelse(grepl("FSGS", a$title) & !grepl("MCD", a$title), "FSGS",
           ifelse(grepl("-MCD|MCD[0-9]+", a$title), "MCD", "Other"))))
    lab[grepl("-FSGS&MCD", a$title)] <- "Other"
    keep <- lab %in% c("MN", "FSGS", "MCD", "Other")
    X <- t(sub[, a$gsm[keep], drop = FALSE])
    miss <- setdiff(panel, colnames(X))
    if (length(miss)) X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
    obj <- list(X = X[, panel, drop = FALSE],
                y = factor(ifelse(lab[keep] == "MN", "MN", "rest"), levels = c("rest", "MN")))
    cat(sprintf("GSE99339: %d glomerular samples (%d MN), genes mapped %d/%d, missing %s\n",
                nrow(X), sum(obj$y == "MN"), ncol(X) - length(miss), length(panel),
                if (length(miss)) paste(miss, collapse = ",") else "none"))
    robust_row <- eval_auc(obj, "GSE99339 (ERCB patients, other CDF revision)",
                           "annotation-robustness: same cohort, different Brainarray revision")
  }
}
if (!is.null(robust_row)) rows[["annot"]] <- robust_row

rob <- do.call(rbind, rows)
write.csv(rob, file.path(tab_dir, "09c_panel_robustness.csv"), row.names = FALSE)
cat("\n--- robustness / transfer ---\n"); print(rob, row.names = FALSE)

# --- A. recalibration on a held-out half of ERCB ---------------------------
ercb <- load_set("val_glom_array")
set.seed(SEED)
idx <- unlist(lapply(split(seq_along(ercb$y), ercb$y), function(i) sample(i, floor(length(i) / 2))))
fit_half <- setdiff(seq_along(ercb$y), idx)

s_ercb <- score(ercb$X)
logit <- function(p) log(pmax(pmin(p, 1 - 1e-6), 1e-6) / (1 - pmax(pmin(p, 1 - 1e-6), 1e-6)))

cal_fit <- glm(ercb$y[fit_half] ~ logit(s_ercb[fit_half]), family = binomial)
ab <- coef(cal_fit)
s_cal <- plogis(ab[1] + ab[2] * logit(s_ercb[idx]))

brier <- function(y, s) mean((as.integer(y == "MN") - s)^2)
cal_rows <- rbind(
  data.frame(half = "held-out ERCB", version = "before recalibration",
             calibration_in_large = round(mean(s_ercb[idx]) - mean(ercb$y[idx] == "MN"), 4),
             calibration_slope = round(coef(glm(ercb$y[idx] ~ s_ercb[idx], family = binomial))[2], 4),
             brier = round(brier(ercb$y[idx], s_ercb[idx]), 4),
             auc = round(as.numeric(auc(roc(ercb$y[idx], s_ercb[idx], levels = c("rest", "MN"),
                                           direction = "<", quiet = TRUE))), 4)),
  data.frame(half = "held-out ERCB", version = "after recalibration",
             calibration_in_large = round(mean(s_cal) - mean(ercb$y[idx] == "MN"), 4),
             calibration_slope = round(coef(glm(ercb$y[idx] ~ s_cal, family = binomial))[2], 4),
             brier = round(brier(ercb$y[idx], s_cal), 4),
             auc = round(as.numeric(auc(roc(ercb$y[idx], s_cal, levels = c("rest", "MN"),
                                           direction = "<", quiet = TRUE))), 4)))
write.csv(cal_rows, file.path(tab_dir, "09c_panel_recalibration.csv"), row.names = FALSE)
cat("\n--- recalibration (intercept ", round(ab[1], 3), ", slope ", round(ab[2], 3), ") ---\n", sep = "")
print(cal_rows, row.names = FALSE)

# --- E. cross-platform training -------------------------------------------
rna <- load_set("cross_glom_rnaseq")
zself <- function(X) {
  out <- matrix(0, nrow(X), length(panel), dimnames = list(NULL, panel))
  have <- intersect(panel, colnames(X))
  if (length(have)) {
    zz <- scale(X[, have, drop = FALSE])
    zz[!is.finite(zz)] <- 0
    out[, have] <- zz
  }
  out
}
fit_x <- glm(rna$y ~ ., data = data.frame(zself(rna$X), y = rna$y), family = binomial)
s_x <- as.numeric(predict(fit_x, newdata = as.data.frame(zself(ercb$X)), type = "response"))
auc_x <- as.numeric(auc(roc(ercb$y, s_x, levels = c("rest", "MN"), direction = "<", quiet = TRUE)))
co_tab <- data.frame(gene = panel,
                     coef_discovery_fit = round(coef(fit)[panel], 4),
                     coef_rnaseq_fit = round(coef(fit_x)[panel], 4),
                     sign_concordant = sign(coef(fit)[panel]) == sign(coef(fit_x)[panel]),
                     stringsAsFactors = FALSE)
write.csv(co_tab, file.path(tab_dir, "09c_crossplatform_coefficients.csv"), row.names = FALSE)
cat(sprintf("\ncross-platform training (RNA-seq -> ERCB): AUC %.4f | coefficient sign concordance %d/%d\n",
            auc_x, sum(co_tab$sign_concordant), nrow(co_tab)))

# --- figures ---------------------------------------------------------------
roc_rows <- list()
for (nm in c("ERCB glomerular (reference)", "ERCB tubulointerstitial",
             "LCM glomeruli GSE216841 (third cohort)")) {
  obj <- switch(nm,
                "ERCB glomerular (reference)" = ercb,
                "ERCB tubulointerstitial" = load_set("val_ti_array"),
                "LCM glomeruli GSE216841 (third cohort)" = load_set("pla2r_lcm"))
  r <- roc(obj$y, score(obj$X), levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  roc_rows[[nm]] <- data.frame(cohort = nm, specificity = r$specificities,
                               sensitivity = r$sensitivities)
}
rocs <- do.call(rbind, roc_rows)
labs <- paste0(rob$cohort, ": ", sprintf("%.3f", rob$auc))
names(labs) <- rob$cohort

save_plot(ggplot(rocs, aes(x = 1 - specificity, y = sensitivity, colour = cohort)) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
            geom_line(linewidth = 0.9) + coord_equal() +
            scale_colour_brewer(palette = "Dark2", labels = labs[names(roc_rows)]) +
            labs(title = "Panel transfer across compartments and cohorts",
                 subtitle = sprintf("locked %d-gene model; glomerular signal is expected to weaken in the interstitium",
                                    length(panel)),
                 x = "1 - specificity", y = "sensitivity", colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), legend.position = c(0.66, 0.2),
                  plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09c_robustness_roc"), width = 6.8, height = 5.6)

cal_df <- data.frame(score = c(s_ercb[idx], s_cal),
                     obs = rep(as.integer(ercb$y[idx] == "MN"), 2),
                     version = rep(c("before", "after"), each = length(idx)))
cal_sum <- do.call(rbind, lapply(split(cal_df, cal_df$version), function(d) {
  br <- unique(quantile(d$score, probs = seq(0, 1, length.out = 5)))
  if (length(br) < 3) return(NULL)
  b <- cut(d$score, br, include.lowest = TRUE)
  data.frame(version = d$version[1], pred = tapply(d$score, b, mean),
             obs = tapply(d$obs, b, mean), n = tapply(d$obs, b, length))
}))
cal_sum$se <- sqrt(cal_sum$obs * (1 - cal_sum$obs) / cal_sum$n)

save_plot(ggplot(cal_sum, aes(x = pred, y = obs, colour = version)) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
            geom_errorbar(aes(ymin = pmax(0, obs - 1.96 * se), ymax = pmin(1, obs + 1.96 * se)),
                          width = 0.012) +
            geom_point(size = 2) + geom_line(linewidth = 0.7) +
            coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
            scale_colour_manual(values = c(before = "firebrick", after = "seagreen")) +
            labs(title = "Recalibration on a held-out half of ERCB",
                 subtitle = sprintf("intercept %.2f, slope %.2f fitted on the other half", ab[1], ab[2]),
                 x = "mean predicted probability", y = "observed proportion MN", colour = NULL) +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09c_recalibration"), width = 6, height = 5.2)

cat("\n09c done\n")
