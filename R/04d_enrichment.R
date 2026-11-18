#!/usr/bin/env Rscript
# Paper 1 - Phase B step 9: functional profiling of the meta-DEGs.
#
#   ORA   : GO (BP/CC/MF, org.Hs.eg.db) and Reactome / KEGG-legacy (MSigDB v2025.1.Hs GMT files)
#           on the meta-DEG symbols of every contrast that yields >= 10 meta-DEGs, with the genes
#           actually tested in the meta as the universe.
#   GSEA  : fgsea on the pooled z-score of every gene in every primary meta, for Hallmark (H) and
#           the immunologic collection (C7).
#
# Gene sets are the MSigDB v2025.1.Hs symbol GMT files in work/msigdb (reference annotation, no
# package-level download at run time).
#
# Outputs
#   output/tables/enrichment/ora_<compartment>_<contrast>_{GO,Reactome,KEGG}.csv
#   output/tables/enrichment/gsea_<compartment>_<contrast>_{Hallmark,C7}.csv
#   output/tables/enrichment_summary.csv
#   output/tables/enrichment_expected_signals.csv
#   output/enrichment/04d_*.{pdf,jpeg}

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(fgsea);
                                 library(ggplot2); library(enrichplot) })

meta_dir <- file.path(OUT, "tables", "meta")
enr_dir <- file.path(OUT, "tables", "enrichment")
fig_dir <- file.path(OUT, "enrichment")
gmt_dir <- file.path(PROJ, "work", "msigdb")

dir.create(enr_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

MIN_DEG <- 10
PADJ_CUT <- 0.2

gmt <- c(Hallmark = "h.all.v2025.1.Hs.symbols.gmt",
         C7 = "c7.all.v2025.1.Hs.symbols.gmt",
         Reactome = "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
         KEGG = "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt")

# long format (one row per term-gene pair) so the same object feeds enricher() and fgsea()
read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE)
  ln <- ln[nzchar(ln)]
  parts <- strsplit(ln, "\t")
  bad <- !vapply(parts, function(z) length(z) >= 3, logical(1))
  if (any(bad)) stop("malformed GMT lines in ", f, ": ", sum(bad))
  out <- do.call(rbind, lapply(parts, function(z) {
    data.frame(term = z[1], gene = unique(z[-1]), stringsAsFactors = FALSE)
  }))
  out[!duplicated(out), ]
}

term2gene <- lapply(gmt, function(f) read_gmt(file.path(gmt_dir, f)))
names(term2gene) <- names(gmt)

for (n in names(term2gene)) cat("gene set collection", n, ":",
                                length(unique(term2gene[[n]]$term)), "sets\n")

pathways <- lapply(term2gene, function(d) split(d$gene, d$term))

save_plot <- function(p, stem, width = 8, height = 5.5) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

meta_sum <- read.csv(file.path(OUT, "tables", "meta_summary.csv"), stringsAsFactors = FALSE)

ora_rows <- list()
gsea_rows <- list()
signal_rows <- list()

for (i in seq_len(nrow(meta_sum))) {

  comp <- meta_sum$compartment[i]
  ct <- meta_sum$contrast[i]
  f <- file.path(meta_dir, paste0("meta_", comp, "_", ct, ".csv"))

  if (!file.exists(f)) next

  m <- read.csv(f, stringsAsFactors = FALSE)
  m <- m[m$estimable, ]

  universe <- unique(m$gene)
  degs <- unique(m$gene[m$meta_deg])

  cat("\n=====", comp, ct, "| universe", length(universe), "| meta-DEG", length(degs), "\n")

  # ---------------- ORA ----------------
  if (length(degs) >= MIN_DEG) {

    for (ont in c("BP", "CC", "MF")) {

      res <- enrichGO(gene = degs, universe = universe, OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
                      pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

      d <- as.data.frame(res)

      if (nrow(d)) {
        d$ontology <- ont
        write.csv(d, file.path(enr_dir, paste0("ora_", comp, "_", ct, "_GO_", ont, ".csv")),
                  row.names = FALSE)
        ora_rows[[paste(comp, ct, "GO", ont)]] <- data.frame(
          compartment = comp, contrast = ct, collection = paste0("GO_", ont),
          n_terms = nrow(d), n_significant = sum(d$p.adjust < PADJ_CUT),
          top_term = d$Description[which.min(d$p.adjust)],
          top_padj = signif(min(d$p.adjust), 3), stringsAsFactors = FALSE)
      }
    }

    # GO BP dotplot is the readable summary plate
    res_bp <- enrichGO(gene = degs, universe = universe, OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff = PADJ_CUT, qvalueCutoff = PADJ_CUT, readable = TRUE)

    if (!is.null(res_bp) && nrow(as.data.frame(res_bp))) {
      p <- dotplot(res_bp, showCategory = 15, font.size = 9) +
        labs(title = paste0(comp, " ", sub("_", " ", ct), " - GO biological process")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      save_plot(p, file.path(fig_dir, paste0("04d_", comp, "_", ct, "_GO_BP_dotplot")),
                width = 9, height = 6)
    }

    for (coll in c("Reactome", "KEGG")) {

      res <- enricher(gene = degs, universe = universe, TERM2GENE = term2gene[[coll]],
                      pvalueCutoff = 1, qvalueCutoff = 1)

      d <- as.data.frame(res)

      if (nrow(d)) {
        write.csv(d, file.path(enr_dir, paste0("ora_", comp, "_", ct, "_", coll, ".csv")),
                  row.names = FALSE)
        ora_rows[[paste(comp, ct, coll)]] <- data.frame(
          compartment = comp, contrast = ct, collection = coll, n_terms = nrow(d),
          n_significant = sum(d$p.adjust < PADJ_CUT),
          top_term = d$ID[which.min(d$p.adjust)],
          top_padj = signif(min(d$p.adjust), 3), stringsAsFactors = FALSE)

        dd <- d[d$p.adjust < PADJ_CUT, ]

        if (nrow(dd)) {
          dd <- head(dd[order(dd$p.adjust), ], 15)
          dd$Description <- factor(dd$ID, levels = rev(dd$ID))
          p <- ggplot(dd, aes(x = -log10(p.adjust), y = Description, size = Count,
                              colour = -log10(p.adjust))) +
            geom_point() +
            scale_colour_gradient(low = "#7FB3D5", high = "#C0392B") +
            labs(title = paste0(comp, " ", sub("_", " ", ct), " - ", coll),
                 x = "-log10 BH-adjusted P", y = NULL) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
          save_plot(p, file.path(fig_dir, paste0("04d_", comp, "_", ct, "_", coll, "_dotplot")),
                    width = 9, height = 5.5)
        }
      }
    }

  }

  # ---------------- GSEA ----------------
  if (sum(m$estimable) >= 2000) {

    stats <- m$zval
    names(stats) <- m$gene
    stats <- stats[is.finite(stats)]
    stats <- stats[order(-stats)]

    for (coll in c("Hallmark", "C7")) {

      set.seed(20261018)

      g <- fgsea(pathways = pathways[[coll]], stats = stats, minSize = 15, maxSize = 500,
                 eps = 0, nproc = 1)

      if (is.null(g) || !nrow(g)) next

      g <- as.data.frame(g)
      g$leadingEdge <- vapply(g$leadingEdge, paste, "", collapse = "/")
      g <- g[order(g$padj), ]
      write.csv(g, file.path(enr_dir, paste0("gsea_", comp, "_", ct, "_", coll, ".csv")),
                row.names = FALSE)

      gsea_rows[[paste(comp, ct, coll)]] <- data.frame(
        compartment = comp, contrast = ct, collection = coll, n_sets = nrow(g),
        n_significant = sum(g$padj < 0.05, na.rm = TRUE),
        n_up_in_MN = sum(g$padj < 0.05 & g$NES > 0, na.rm = TRUE),
        n_down_in_MN = sum(g$padj < 0.05 & g$NES < 0, na.rm = TRUE),
        top_set = g$pathway[1], top_padj = signif(g$padj[1], 3),
        stringsAsFactors = FALSE)

      if (ct %in% c("MN_vs_Control", "MN_vs_nonMN")) {

        sig <- g[!is.na(g$padj) & g$padj < 0.05, ]

        if (nrow(sig)) {

          dd <- head(sig[order(-abs(sig$NES)), ], 12)
          dd$pathway <- factor(dd$pathway, levels = dd$pathway[order(dd$NES)])

          p <- ggplot(dd, aes(NES, pathway, fill = NES > 0)) +
            geom_col() +
            scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2471A3"),
                              labels = c("down in MN", "up in MN"), name = NULL) +
            labs(title = paste0(comp, " ", sub("_", " ", ct), " - GSEA ", coll),
                 x = "normalised enrichment score (NES)", y = NULL) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

          save_plot(p, file.path(fig_dir, paste0("04d_", comp, "_", ct, "_GSEA_", coll, "_bar")),
                    width = 9.5, height = 5.5)

          if (coll == "Hallmark") {
            top3 <- head(sig$pathway[order(sig$padj)], 3)
            for (tp in top3) {
              pp <- plotEnrichment(pathways[[coll]][[tp]], stats) +
                labs(title = paste0(comp, " ", sub("_", " ", ct), " - ", tp),
                     x = "gene rank", y = "enrichment score") +
                theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              save_plot(pp, file.path(fig_dir, paste0("04d_", comp, "_", ct, "_GSEA_",
                                                      gsub("[^A-Za-z0-9]+", "_", tp))),
                        width = 9, height = 6)
            }
          }
        }
      }
    }
  }
}

# --- targeted check of the signals the plan expects to see -------------------
# Read back the written tables rather than the in-loop objects: the file name carries the
# compartment / contrast / collection, and matching is on the term label itself.
signals <- c(complement = "complement", fc_gamma_receptor = "Fc.?gamma",
             antigen_presentation = "antigen", ecm_collagen = "extracellular matrix|collagen",
             hypoxia = "hypoxia", phagocytosis = "phagocyt", immunoglobulin = "immunoglobulin")

ora_files <- list.files(enr_dir, pattern = "^ora_", full.names = TRUE)

for (f in ora_files) {

  nm <- sub("\\.csv$", "", basename(f))
  m <- regmatches(nm, regexec("^ora_(Glom|TI)_(MN_vs_[A-Za-z]+)_(GO_BP|GO_CC|GO_MF|Reactome|KEGG)$",
                              nm))[[1]]

  if (length(m) != 4) next

  d <- read.csv(f, stringsAsFactors = FALSE)
  lab <- if ("Description" %in% names(d)) d$Description else d$ID
  padj <- d$p.adjust

  for (k in names(signals)) {

    idx <- grep(signals[k], lab, ignore.case = TRUE, perl = TRUE)

    if (!length(idx)) next

    best <- idx[which.min(padj[idx])]

    signal_rows[[paste(m[2], m[3], m[4], k)]] <- data.frame(
      compartment = m[2], contrast = m[3], collection = m[4], keyword = k,
      best_term = lab[best], p_adj = signif(padj[best], 3),
      significant = padj[best] < PADJ_CUT, source = "ORA", stringsAsFactors = FALSE)
  }
}

gsea_files <- list.files(enr_dir, pattern = "^gsea_", full.names = TRUE)

for (f in gsea_files) {

  nm <- sub("\\.csv$", "", basename(f))
  m <- regmatches(nm, regexec("^gsea_(Glom|TI)_(MN_vs_[A-Za-z]+)_(Hallmark|C7)$", nm))[[1]]

  if (length(m) != 4) next

  d <- read.csv(f, stringsAsFactors = FALSE)

  for (k in names(signals)) {

    idx <- grep(signals[k], d$pathway, ignore.case = TRUE, perl = TRUE)

    if (!length(idx)) next

    best <- idx[which.min(d$padj[idx])]

    signal_rows[[paste(m[2], m[3], m[4], k)]] <- data.frame(
      compartment = m[2], contrast = m[3], collection = m[4], keyword = k,
      best_term = d$pathway[best], p_adj = signif(d$padj[best], 3),
      significant = !is.na(d$padj[best]) & d$padj[best] < 0.05, source = "GSEA",
      stringsAsFactors = FALSE)
  }
}

ora_sum <- if (length(ora_rows)) do.call(rbind, ora_rows) else NULL
gsea_sum <- if (length(gsea_rows)) do.call(rbind, gsea_rows) else NULL
sig_sum <- if (length(signal_rows)) do.call(rbind, signal_rows) else NULL

if (!is.null(ora_sum)) write.csv(ora_sum, file.path(OUT, "tables", "enrichment_summary_ora.csv"),
                                row.names = FALSE)

if (!is.null(gsea_sum)) write.csv(gsea_sum, file.path(OUT, "tables", "enrichment_summary_gsea.csv"),
                                  row.names = FALSE)

if (!is.null(sig_sum)) write.csv(sig_sum, file.path(OUT, "tables", "enrichment_expected_signals.csv"),
                                 row.names = FALSE)

cat("\n================ ORA SUMMARY ================\n")

if (!is.null(ora_sum)) print(ora_sum, row.names = FALSE)

cat("\n================ GSEA SUMMARY ================\n")

if (!is.null(gsea_sum)) print(gsea_sum, row.names = FALSE)

cat("\n================ EXPECTED SIGNALS ================\n")

if (!is.null(sig_sum)) print(sig_sum, row.names = FALSE)

writeLines(c(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 04d_enrichment.R"),
             capture.output(print(ora_sum, row.names = FALSE)),
             capture.output(print(gsea_sum, row.names = FALSE))),
           file.path(PROJ, "logs", "enrichment.log"))

cat("\nDone.\n")
