#!/usr/bin/env Rscript
# Paper 1 - Phase E step 6: PLA2R-stratified analysis of the panel (resolves handoff B12).
#
# B12 recorded PLA2R status as unavailable. It is in fact present in the GSE216841 series matrix as a
# per-sample characteristic ("phenotype: PLA2R_Positive" / "PLA2R_Negative"), spread over two
# !Sample_characteristics_ch1 rows; 01a did not harvest it. This script parses it, writes a reusable
# table, and asks whether the panel tracks PLA2R-associated disease.
#
# Outputs
#   metadata/GSE216841_pla2r.tsv                     gsm -> PLA2R status (reusable annotation)
#   output/tables/09f_pla2r_strata.csv               panel AUC by stratum
#   output/tables/09f_pla2r_expression.csv           the five genes by stratum
#   output/models_panel/09f_pla2r_strata.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(pROC); library(ggplot2) })

set_dir <- file.path(PROJ, "work", "sets")
tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "models_panel")
meta_dir <- file.path(PROJ, "metadata")
SEED <- 20261019L
set.seed(SEED)

panel <- read.csv(file.path(tab_dir, "09a_panel_consensus.csv"), stringsAsFactors = FALSE)$gene

# --- parse the PLA2R characteristic ---------------------------------------
sm_f <- file.path(DATA, "GSE216841", "GSE216841_series_matrix.txt.gz")
stopifnot(file.exists(sm_f))
ln <- readLines(gzfile(sm_f), warn = FALSE)
get_row <- function(key) {
  x <- ln[startsWith(ln, key)]
  lapply(x, function(z) gsub('^"|"$', "", strsplit(z, "\t")[[1]][-1]))
}
gsm <- get_row("!Sample_geo_accession")[[1]]
# each !Sample_characteristics_ch1 LINE is one vector with an entry per sample, and a sample may
# carry tissue / cell-type / phenotype on different lines - so iterate over lines, not over tokens
char_rows <- get_row("!Sample_characteristics_ch1")
pla2r <- rep(NA_character_, length(gsm))
for (r in char_rows) {
  if (length(r) != length(gsm)) next
  hit <- startsWith(r, "phenotype:")
  if (any(hit)) pla2r[hit] <- sub("^phenotype:\\s*", "", r[hit])
}
# the processed count matrix is keyed by the submitters' sample ids ("MEMBR1", "MCD3", "N_CTRL1"),
# not by GSM, so the series matrix descriptions are the join key:
#   "MEMBR1.MN.PLA2Rpos" -> MEMBR1 | "MCD3.MCD" -> MCD3 | "G3.N_CTRL" -> N_CTRL3
desc <- get_row("!Sample_description")[[1]]
sample_key <- ifelse(grepl("N_CTRL$", desc),
                     paste0("N_CTRL", sub("^G([0-9]+)\\..*$", "\\1", desc)),
                     sub("\\..*$", "", desc))
pla2r_df <- data.frame(gsm = gsm, sample_key = sample_key, PLA2R = pla2r,
                       stringsAsFactors = FALSE)
cat("PLA2R parsed:", sum(!is.na(pla2r_df$PLA2R)), "of", nrow(pla2r_df), "samples\n")
print(table(pla2r_df$PLA2R, useNA = "ifany"))

dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)
write.table(pla2r_df, file.path(meta_dir, "GSE216841_pla2r.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# --- score the cohort ------------------------------------------------------
load_set <- function(set) {
  ms <- readRDS(file.path(set_dir, paste0(set, "_expr_symbol.rds")))
  ann <- read.csv(file.path(set_dir, paste0(set, "_samples.csv")), stringsAsFactors = FALSE)
  stopifnot(identical(colnames(ms), ann$gsm))
  grp <- ifelse(ann$disease %in% c("Control-LD", "Control-TN"), "Control", ann$disease)
  keep <- !(grp %in% c("Control", "UNKNOWN"))
  list(ms = ms[, keep, drop = FALSE],
       y = factor(ifelse(grp[keep] == "MN", "MN", "rest"), levels = c("rest", "MN")),
       gsm = ann$gsm[keep], raw = grp[keep])
}
disc <- load_set("disc_glom_array")
lcm  <- load_set("pla2r_lcm")

mu  <- rowMeans(disc$ms[panel, , drop = FALSE])
sdv <- apply(disc$ms[panel, , drop = FALSE], 1, sd)
sdv[!is.finite(sdv) | sdv == 0] <- 1
z <- function(ms) {
  out <- matrix(0, ncol(ms), length(panel), dimnames = list(NULL, panel))
  have <- intersect(panel, rownames(ms))
  zz <- sweep(sweep(t(ms[have, , drop = FALSE]), 2, mu[have]), 2, sdv[have], "/")
  zz[!is.finite(zz)] <- 0
  out[, have] <- zz
  out
}
fit <- glm(disc$y ~ ., data = data.frame(z(disc$ms), y = disc$y), family = binomial)

s <- as.numeric(predict(fit, newdata = as.data.frame(z(lcm$ms)), type = "response"))
status <- pla2r_df$PLA2R[match(lcm$gsm, pla2r_df$sample_key)]
cat("\nLCM cohort:", length(s), "samples |", sum(is.na(status)), "without a PLA2R call\n")
print(table(lcm$raw, status, useNA = "ifany"))

# --- AUC by stratum --------------------------------------------------------
strata <- list(
  list(nm = "MN vs MCD (all)", idx = lcm$y %in% c("MN", "rest")),
  list(nm = "MN PLA2R+ vs MCD", idx = (lcm$raw == "MN" & status == "PLA2R_Positive") |
                                         lcm$raw == "MCD"),
  list(nm = "MN PLA2R- vs MCD", idx = (lcm$raw == "MN" & status == "PLA2R_Negative") |
                                         lcm$raw == "MCD"))

rows <- list()
for (st in strata) {
  idx <- st$idx & !is.na(st$idx)
  if (sum(idx) < 6) next
  yy <- lcm$y[idx]
  yy <- factor(ifelse(yy == "MN", "MN", "rest"), levels = c("rest", "MN"))
  if (length(unique(yy)) < 2) next
  r <- roc(yy, s[idx], levels = c("rest", "MN"), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  rows[[st$nm]] <- data.frame(stratum = st$nm, n = sum(idx), n_MN = sum(yy == "MN"),
                              n_comparator = sum(yy == "rest"),
                              auc = round(as.numeric(auc(r)), 4),
                              ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
                              median_score_MN = round(median(s[idx][yy == "MN"]), 4),
                              median_score_other = round(median(s[idx][yy == "rest"]), 4),
                              stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)

# direct test of PLA2R+ versus PLA2R- within MN
mn_idx <- lcm$raw == "MN" & !is.na(status)
if (length(unique(status[mn_idx])) > 1) {
  yy <- factor(status[mn_idx], levels = c("PLA2R_Negative", "PLA2R_Positive"))
  r <- roc(yy, s[mn_idx], levels = c("PLA2R_Negative", "PLA2R_Positive"), direction = "<",
           quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  res <- rbind(res, data.frame(
    stratum = "within MN: PLA2R+ vs PLA2R-", n = sum(mn_idx),
    n_MN = sum(status[mn_idx] == "PLA2R_Positive"),
    n_comparator = sum(status[mn_idx] == "PLA2R_Negative"),
    auc = round(as.numeric(auc(r)), 4), ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
    median_score_MN = round(median(s[mn_idx][status[mn_idx] == "PLA2R_Positive"]), 4),
    median_score_other = round(median(s[mn_idx][status[mn_idx] == "PLA2R_Negative"]), 4),
    stringsAsFactors = FALSE))
}
write.csv(res, file.path(tab_dir, "09f_pla2r_strata.csv"), row.names = FALSE)
cat("\n--- panel AUC by PLA2R stratum ---\n"); print(res, row.names = FALSE)

# --- the five genes by stratum --------------------------------------------
expr <- do.call(rbind, lapply(panel, function(g) {
  v <- as.numeric(lcm$ms[g, ])
  grp <- ifelse(lcm$raw == "MN" & status == "PLA2R_Positive", "MN PLA2R+",
         ifelse(lcm$raw == "MN" & status == "PLA2R_Negative", "MN PLA2R-",
         ifelse(lcm$raw == "MCD", "MCD", "Control")))
  data.frame(gene = g, group = grp, value = v, stringsAsFactors = FALSE)
}))
expr <- expr[expr$group != "Control" | TRUE, ]
agg <- aggregate(value ~ gene + group, data = expr, FUN = median)
names(agg)[3] <- "median_expr"
write.csv(agg, file.path(tab_dir, "09f_pla2r_expression.csv"), row.names = FALSE)

save_plot <- function(p, stem, width = 8, height = 5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

res$label <- factor(res$stratum, levels = res$stratum[order(res$auc)])
save_plot(ggplot(res, aes(x = label, y = auc, fill = grepl("within MN", label))) +
            geom_col(width = 0.66, show.legend = FALSE) +
            geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
            geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
            geom_text(aes(label = paste0("n=", n)), vjust = -0.4, size = 3, colour = "grey30") +
            coord_flip(ylim = c(0.3, 1.05)) +
            scale_fill_manual(values = c(`FALSE` = "grey45", `TRUE` = "steelblue")) +
            labs(title = "Panel performance by PLA2R stratum (LCM cohort)",
                 subtitle = "GSE216841 phenotype characteristic; MN PLA2R+ / PLA2R- against MCD, and within MN",
                 x = NULL, y = "AUC (95% CI)") +
            theme_bw(base_size = 10) +
            theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")),
          file.path(fig_dir, "09f_pla2r_strata"), width = 7.4, height = 4.4)

cat("\n09f done\n")
