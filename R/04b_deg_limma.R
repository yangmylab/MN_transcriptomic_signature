#!/usr/bin/env Rscript
# Paper 1 - Phase B step 7: differential expression with limma, per cohort, per compartment.
#
# Contrasts (the plan asks for MN vs healthy donors AND MN vs each other glomerulopathy
# separately, per compartment):
#   MN_vs_Control   MN vs healthy living/tumour-nephrectomy donors
#   MN_vs_FSGS      MN vs FSGS
#   MN_vs_MCD       MN vs minimal-change disease
#   MN_vs_Other     MN vs every remaining diagnosis (IgAN, SLE, DN, HT, RPGN, TMD, FSGS_MCD, "Other")
#   MN_vs_nonMN     binary model, MN vs every non-MN sample including controls
#
# The discovery tubulointerstitial matrix is the ComBat-adjusted one from 04a; a sensitivity table
# compares its logFC to the un-adjusted matrix so the correction can be defended.
# Thresholds: |log2FC| > 1 and BH-adjusted P < 0.05.
# QC policy: every sample is retained (R/03d_qc_policy.R); outlier stage is not applied here.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(limma) })

set_dir <- file.path(PROJ, "work", "sets")
combat_dir <- file.path(PROJ, "work", "combat")
deg_dir <- file.path(OUT, "tables", "deg")
tab_dir <- file.path(OUT, "tables")

dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

LFC <- 1
FDR <- 0.05
MIN_N <- 5
DETECT_FRAC <- 0.25

harmonise <- function(x) {
  out <- ifelse(x %in% c("Control-LD", "Control-TN"), "Control", x)
  out[!(out %in% c("MN", "FSGS", "MCD", "Control"))] <- "Other"
  out
}

# detection filter: log2 intensity > 5 for arrays, log2 CPM > 1 for RNA-seq, in >= 25 % of samples
detected <- function(ms, assay) {
  thr <- if (identical(assay, "microarray")) 5 else 1
  frac <- rowMeans(ms > thr)
  keep <- frac >= DETECT_FRAC
  # releases that combine two array types carry NA for genes absent from one of them
  # (P-10 in the handoff): such genes cannot enter limma and are dropped.
  complete <- rowSums(is.na(ms)) == 0
  keep & complete & apply(ms, 1, sd, na.rm = TRUE) > 0
}

fit_contrast <- function(fit0, design, num, den) {
  cm <- makeContrasts(contrasts = paste(num, "-", den), levels = design)
  eBayes(contrasts.fit(fit0, cm), trend = TRUE, robust = TRUE)
}

parse_contrast <- function(fit) {
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, 1]
  tt$SE <- se[rownames(tt)]
  tt$gene <- rownames(tt)
  tt[, c("gene", "logFC", "SE", "t", "P.Value", "adj.P.Val", "B", "AveExpr")]
}

config <- data.frame(
  set = c("disc_glom_array", "disc_ti_array", "val_glom_array", "val_ti_array",
          "cross_glom_rnaseq", "cross_ti_rnaseq", "pla2r_lcm"),
  compartment = c("Glom", "TI", "Glom", "TI", "Glom", "TI", "Glom"),
  role = c("discovery", "discovery", "independent validation", "independent validation",
           "independent RNA-seq cohort", "independent RNA-seq cohort", "PLA2R sub-analysis"),
  stringsAsFactors = FALSE)

# one pass: fit every cohort and every contrast
deg_rows <- list()
sum_rows <- list()
sensitivity <- list()

for (i in seq_len(nrow(config))) {

  s <- config$set[i]
  cat("\n=====", s, "(", config$role[i], ") =====\n")

  ann <- read.csv(file.path(set_dir, paste0(s, "_samples.csv")), stringsAsFactors = FALSE)

  src <- if (s == "disc_ti_array") file.path(combat_dir, paste0(s, "_expr_combat_symbol.rds")) else
    file.path(set_dir, paste0(s, "_expr_symbol.rds"))

  ms <- readRDS(src)

  stopifnot(identical(colnames(ms), ann$gsm))

  keep <- detected(ms, ann$assay[1])

  ms <- ms[keep, , drop = FALSE]

  grp <- harmonise(ann$disease)
  ok <- !is.na(grp) & ann$disease != "UNKNOWN"

  if (any(!ok)) cat("  dropped samples with unknown diagnosis:", sum(!ok), "\n")

  ms <- ms[, ok, drop = FALSE]
  grp <- grp[ok]

  cat("  samples:", ncol(ms), "| genes after detection filter:", nrow(ms), "\n")
  print(table(grp))

  lev <- c("Control", "FSGS", "MCD", "Other", "MN")

  design <- model.matrix(~ 0 + factor(grp, levels = lev))
  colnames(design) <- lev

  fit <- lmFit(ms, design)

  # keep absent levels (cross_ti_rnaseq has no healthy donors) so the availability test is defined
  avail <- table(factor(grp, levels = lev))

  pairs <- list(c("MN_vs_Control", "MN", "Control"),
                c("MN_vs_FSGS", "MN", "FSGS"),
                c("MN_vs_MCD", "MN", "MCD"),
                c("MN_vs_Other", "MN", "Other"))

  for (p in pairs) {

    nm <- p[1]

    if (avail[p[2]] < MIN_N || avail[p[3]] < MIN_N) {
      cat("  skip", nm, "- group smaller than", MIN_N, "\n")
      next
    }

    tt <- parse_contrast(fit_contrast(fit, design, p[2], p[3]))
    tt$deg <- abs(tt$logFC) > LFC & tt$adj.P.Val < FDR
    tt$direction <- ifelse(tt$deg, ifelse(tt$logFC > 0, "up", "down"), "ns")

    write.csv(tt, file.path(deg_dir, paste0(s, "__", nm, ".csv")), row.names = FALSE)

    sig <- tt[tt$deg, , drop = FALSE]

    if (nrow(sig)) {
      deg_rows[[paste(s, nm)]] <- data.frame(set = s, compartment = config$compartment[i],
                                             contrast = nm, sig, stringsAsFactors = FALSE)
    }

    sum_rows[[paste(s, nm)]] <- data.frame(
      set = s, compartment = config$compartment[i], role = config$role[i], contrast = nm,
      n_case = as.integer(avail[p[2]]), n_ref = as.integer(avail[p[3]]),
      genes_tested = nrow(tt), deg_up = sum(tt$deg & tt$logFC > 0),
      deg_down = sum(tt$deg & tt$logFC < 0), deg_total = sum(tt$deg),
      stringsAsFactors = FALSE)

    cat("  ", nm, ": tested", nrow(tt), "| DEG", sum(tt$deg),
        "( up", sum(tt$deg & tt$logFC > 0), "/ down", sum(tt$deg & tt$logFC < 0), ")\n")
  }

  # binary model for MN vs all non-MN
  if (avail["MN"] >= MIN_N && (ncol(ms) - avail["MN"]) >= MIN_N) {

    is_mn <- factor(ifelse(grp == "MN", "MN", "nonMN"), levels = c("nonMN", "MN"))
    d2 <- model.matrix(~ 0 + is_mn)
    colnames(d2) <- levels(is_mn)

    fit2 <- lmFit(ms, d2)

    tt <- parse_contrast(fit_contrast(fit2, d2, "MN", "nonMN"))
    tt$deg <- abs(tt$logFC) > LFC & tt$adj.P.Val < FDR
    tt$direction <- ifelse(tt$deg, ifelse(tt$logFC > 0, "up", "down"), "ns")

    write.csv(tt, file.path(deg_dir, paste0(s, "__MN_vs_nonMN.csv")), row.names = FALSE)

    sig <- tt[tt$deg, , drop = FALSE]

    if (nrow(sig)) {
      deg_rows[[paste(s, "MN_vs_nonMN")]] <- data.frame(
        set = s, compartment = config$compartment[i], contrast = "MN_vs_nonMN", sig,
        stringsAsFactors = FALSE)
    }

    sum_rows[[paste(s, "MN_vs_nonMN")]] <- data.frame(
      set = s, compartment = config$compartment[i], role = config$role[i], contrast = "MN_vs_nonMN",
      n_case = as.integer(avail["MN"]), n_ref = as.integer(ncol(ms) - avail["MN"]),
      genes_tested = nrow(tt), deg_up = sum(tt$deg & tt$logFC > 0),
      deg_down = sum(tt$deg & tt$logFC < 0), deg_total = sum(tt$deg),
      stringsAsFactors = FALSE)

    cat("   MN_vs_nonMN: tested", nrow(tt), "| DEG", sum(tt$deg), "\n")
  }

  # sensitivity: ComBat vs raw discovery TI, for every contrast already written
  if (s == "disc_ti_array") {

    raw <- readRDS(file.path(set_dir, paste0(s, "_expr_symbol.rds")))
    raw <- raw[keep, ok, drop = FALSE]

    design_raw <- model.matrix(~ 0 + factor(grp, levels = lev))
    colnames(design_raw) <- lev

    fit_raw <- lmFit(raw, design_raw)

    for (nm in c("MN_vs_Control", "MN_vs_FSGS", "MN_vs_MCD", "MN_vs_Other")) {

      p <- switch(nm, MN_vs_Control = c("MN", "Control"), MN_vs_FSGS = c("MN", "FSGS"),
                  MN_vs_MCD = c("MN", "MCD"), MN_vs_Other = c("MN", "Other"))

      a <- read.csv(file.path(deg_dir, paste0(s, "__", nm, ".csv")), stringsAsFactors = FALSE)
      b <- parse_contrast(fit_contrast(fit_raw, design_raw, p[1], p[2]))
      b$deg <- abs(b$logFC) > LFC & b$adj.P.Val < FDR

      common <- intersect(a$gene, b$gene)
      a <- a[match(common, a$gene), ]
      b <- b[match(common, b$gene), ]

      sensitivity[[nm]] <- data.frame(
        contrast = nm, genes = length(common),
        lfc_pearson = round(cor(a$logFC, b$logFC), 4),
        lfc_spearman = round(cor(a$logFC, b$logFC, method = "spearman"), 4),
        deg_combat = sum(a$deg), deg_raw = sum(b$deg),
        deg_overlap = sum(a$deg & b$deg),
        deg_jaccard = round(sum(a$deg & b$deg) / max(1, sum(a$deg | b$deg)), 4),
        stringsAsFactors = FALSE)
    }
  }
}

deg_all <- do.call(rbind, deg_rows)
deg_sum <- do.call(rbind, sum_rows)

write.csv(deg_sum, file.path(tab_dir, "deg_summary.csv"), row.names = FALSE)
write.csv(deg_all, file.path(tab_dir, "deg_all_significant.csv"), row.names = FALSE)

if (length(sensitivity)) {
  write.csv(do.call(rbind, sensitivity), file.path(tab_dir, "deg_combat_sensitivity.csv"),
            row.names = FALSE)
}

cat("\n================ DEG SUMMARY ================\n")

print(deg_sum, row.names = FALSE)

if (length(sensitivity)) {

  cat("\n---- ComBat sensitivity ----\n")

  print(do.call(rbind, sensitivity), row.names = FALSE)
}

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 04b_deg_limma.R"),
             capture.output(print(deg_sum, row.names = FALSE))),
           file.path(PROJ, "logs", "deg_limma.log"))

cat("\nDone. Tables in", deg_dir, "\n")
