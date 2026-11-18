#!/usr/bin/env Rscript
# Paper 1 - Phase D follow-up, step 08f: what IS each co-expression module?
#
# Phase D located the MN signal in two small modules (glom pink, 84 genes; glom red, 144 genes) but
# left them unnamed. A colour is not a figure panel, and a panel gene ("ECM1 belongs to pink") is
# only interpretable if pink has a biological identity. This step runs ORA on every non-grey module
# against the same collections Phase B used, then assigns each module a short programme label from
# its strongest term - the join key that Phase E carries into the panel table.
#
# Collections: GO BP/CC/MF (org.Hs.eg.db), Reactome and KEGG-legacy (MSigDB v2025.1.Hs GMTs in
# work/msigdb), Hallmark (H). Universe = every gene in that compartment's network, i.e. exactly the
# genes that had a chance of being in a module.
#
# Outputs
#   output/tables/08f_module_enrichment_<comp>.csv        every module x collection x term
#   output/tables/08f_module_enrichment_top.csv           top 5 terms per module (paper table)
#   output/tables/08f_module_programme_labels.csv         one row per module: label + MN association
#   output/wgcna/08f_<comp>_module_dotplot.{pdf,jpeg}     top terms per module
#   output/wgcna/08f_<comp>_mn_modules.{pdf,jpeg}         the MN-associated modules in detail

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2) })

tab_dir <- file.path(OUT, "tables")
fig_dir <- file.path(OUT, "wgcna")
gmt_dir <- file.path(PROJ, "work", "msigdb")

MIN_MODULE <- 30L
PADJ_CUT   <- 0.05
N_TOP      <- 5L
TRAIT_CELLS <- c("Podocyte", "PT", "TAL", "DCT", "CD_PC", "PEC", "Endothelium",
                 "Fibroblast", "Mesangium", "Myeloid", "B", "T_NK", "Plasma", "Mast")

save_png <- function(p, stem, width = 9, height = 6) {
  if (interactive()) print(p)
  ggsave(paste0(stem, ".pdf"), plot = p, device = "pdf", width = width, height = height)
  ggsave(paste0(stem, ".jpeg"), plot = p, device = "jpeg", quality = 95, dpi = 300,
         bg = "white", width = width, height = height)
}

gmt <- c(Hallmark = "h.all.v2025.1.Hs.symbols.gmt",
         Reactome = "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
         KEGG = "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt")

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

ALL_F <- file.path(tab_dir, "08f_module_enrichment_all.csv")
RESUME <- file.exists(ALL_F) && !("--force" %in% commandArgs(TRUE))
if (RESUME) {
  cat("resuming from", ALL_F, "- delete it or pass --force to recompute the ORA pass\n")
  long <- read.csv(ALL_F, stringsAsFactors = FALSE)
  labels <- read.csv(file.path(tab_dir, "08f_module_programme_labels.csv"),
                     stringsAsFactors = FALSE)
} else {
long_rows <- list()
label_rows <- list()

for (comp in c("glom", "ti")) {
  cat("\n================ ", toupper(comp), " ================\n", sep = "")

  mods <- read.csv(file.path(tab_dir, paste0("wgcna_modules_", comp, ".csv")), stringsAsFactors = FALSE)
  mt <- read.csv(file.path(tab_dir, paste0("wgcna_module_trait_", comp, ".csv")), stringsAsFactors = FALSE)
  adj <- read.csv(file.path(tab_dir, paste0("wgcna_module_trait_adjusted_", comp, ".csv")),
                  stringsAsFactors = FALSE)

  universe <- unique(mods$gene)
  cat("universe:", length(universe), "genes |", length(unique(mods$module)) - 1, "non-grey modules\n")

  mod_sizes <- table(mods$module)
  test_mods <- names(mod_sizes)[names(mod_sizes) != "grey" & mod_sizes >= MIN_MODULE]
  cat("modules tested:", paste(test_mods, collapse = ","), "\n")

  for (m in test_mods) {
    genes <- unique(mods$gene[mods$module == m])
    rows <- list()

    # --- GO (three ontologies) ---------------------------------------------
    for (ont in c("BP", "CC", "MF")) {
      res <- tryCatch(enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
                               pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
                      error = function(e) NULL)
      d <- if (is.null(res)) NULL else as.data.frame(res)
      if (!is.null(d) && nrow(d)) {
        rows[[length(rows) + 1]] <- data.frame(
          compartment = comp, module = m, n_module_genes = length(genes),
          collection = paste0("GO_", ont), term_id = d$ID, term = d$Description,
          gene_ratio = d$GeneRatio, bg_ratio = d$BgRatio,
          p = d$pvalue, padj = p.adjust(d$pvalue, "BH"),
          n_overlap = as.integer(sub("/.*", "", d$GeneRatio)),
          genes = d$geneID, stringsAsFactors = FALSE)
      }
    }

    # --- MSigDB collections ------------------------------------------------
    for (nm in names(term2gene)) {
      res <- tryCatch(enricher(gene = genes, universe = universe, TERM2GENE = term2gene[[nm]],
                               pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
                      error = function(e) NULL)
      d <- if (is.null(res)) NULL else as.data.frame(res)
      if (!is.null(d) && nrow(d)) {
        rows[[length(rows) + 1]] <- data.frame(
          compartment = comp, module = m, n_module_genes = length(genes),
          collection = nm, term_id = d$ID, term = d$ID,
          gene_ratio = d$GeneRatio, bg_ratio = d$BgRatio,
          p = d$pvalue, padj = p.adjust(d$pvalue, "BH"),
          n_overlap = as.integer(sub("/.*", "", d$GeneRatio)),
          genes = d$geneID, stringsAsFactors = FALSE)
      }
    }

    if (!length(rows)) next
    d <- do.call(rbind, rows)
    d <- d[order(d$padj), ]
    long_rows[[paste(comp, m, sep = "_")]] <- d

    top <- head(d[d$padj < PADJ_CUT, ], N_TOP)
    if (!nrow(top)) top <- head(d, 1)

    # 08b stores the numeric module label in `module` and the colour in `module_color`
    pick1 <- function(x) if (length(x)) x[1] else NA_real_
    mn_cor    <- pick1(mt$cor[mt$module_color == m & mt$trait == "MN"])
    mn_padj   <- pick1(mt$padj[mt$module_color == m & mt$trait == "MN"])
    adj_padj  <- pick1(adj$padj[adj$module_color == m])

    label_rows[[paste(comp, m, sep = "_")]] <- data.frame(
      compartment = comp, module = m, n_genes = length(genes),
      programme_label = top$term[1], label_collection = top$collection[1],
      label_padj = signif(top$padj[1], 3),
      second_term = if (nrow(top) > 1) top$term[2] else NA_character_,
      third_term = if (nrow(top) > 2) top$term[3] else NA_character_,
      n_significant_terms = sum(d$padj < PADJ_CUT),
      MN_bicor = round(mn_cor, 4), MN_padj = signif(mn_padj, 3),
      MN_adjusted_padj = signif(adj_padj, 3),
      stringsAsFactors = FALSE)

    cat(sprintf("  %-9s n=%4d  top: %-52s (%s, padj %.2g)  sig terms %d\n",
                m, length(genes), substr(top$term[1], 1, 52), top$collection[1],
                top$padj[1], sum(d$padj < PADJ_CUT)))
  }
}

}

if (!RESUME) {
long <- do.call(rbind, long_rows)
rownames(long) <- NULL
write.csv(long, file.path(tab_dir, "08f_module_enrichment_all.csv"), row.names = FALSE)
cat("\nwrote 08f_module_enrichment_all.csv (", nrow(long), "rows )\n")

top <- do.call(rbind, lapply(split(long, paste(long$compartment, long$module)),
                             function(d) head(d[order(d$padj), ], N_TOP)))
top <- top[order(top$compartment, top$module, top$padj), ]
write.csv(top, file.path(tab_dir, "08f_module_enrichment_top.csv"), row.names = FALSE)
cat("wrote 08f_module_enrichment_top.csv (", nrow(top), "rows )\n")

labels <- do.call(rbind, label_rows)
rownames(labels) <- NULL
labels <- labels[order(labels$compartment, -labels$n_genes), ]
write.csv(labels, file.path(tab_dir, "08f_module_programme_labels.csv"), row.names = FALSE)
cat("wrote 08f_module_programme_labels.csv (", nrow(labels), "rows )\n")
}

# the resume path must rebuild the top table from the saved long table
if (!exists("top")) {
  top <- do.call(rbind, lapply(split(long, paste(long$compartment, long$module)),
                               function(d) head(d[order(d$padj), ], N_TOP)))
  top <- top[order(top$compartment, top$module, top$padj), ]
}

# --- figures ---------------------------------------------------------------
plot_df <- top[top$padj < PADJ_CUT, ]
if (nrow(plot_df)) {
  plot_df$term_short <- ifelse(nchar(plot_df$term) > 58,
                               paste0(substr(plot_df$term, 1, 55), "..."), plot_df$term)
  plot_df$module <- as.character(plot_df$module)
  lvl <- unique(labels$module[order(labels$compartment, -labels$n_genes)])
  lvl <- lvl[!duplicated(lvl)]
  plot_df$module <- factor(plot_df$module, levels = lvl[lvl %in% unique(plot_df$module)])
  plot_df$score <- -log10(plot_df$padj)

  for (comp in c("glom", "ti")) {
    d <- plot_df[plot_df$compartment == comp, ]
    if (!nrow(d)) next
    p <- ggplot(d, aes(x = score, y = reorder(term_short, score), size = n_overlap, colour = collection)) +
      geom_point() +
      facet_wrap(~ module, scales = "free_y", ncol = 2) +
      scale_colour_brewer(palette = "Dark2") +
      labs(title = paste0("Module functional identity - ", comp),
           subtitle = "top BH-significant terms per module (ORA, universe = network genes)",
           x = expression(-log[10]~"BH padj"), y = NULL, size = "genes", colour = "collection") +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
            plot.title = element_text(face = "bold"))
    save_png(p, file.path(fig_dir, paste0("08f_", comp, "_module_dotplot")),
             width = 12, height = max(5, 2.2 * length(unique(d$module)) * 0.45))
  }

  # detail plate for the modules that actually track MN
  mn_mods <- labels$module[labels$compartment == "glom" & !is.na(labels$MN_padj) & labels$MN_padj < 0.05]
  d2 <- plot_df[plot_df$compartment == "glom" &
                  as.character(plot_df$module) %in% mn_mods, ]
  if (nrow(d2)) {
    p2 <- ggplot(d2, aes(x = score, y = reorder(term_short, score), fill = collection)) +
      geom_col(width = 0.72) +
      facet_wrap(~ module, scales = "free_y", ncol = 1) +
      scale_fill_brewer(palette = "Dark2") +
      labs(title = "Functional identity of the MN-associated glomerular modules",
           subtitle = paste0("modules with BH padj < 0.05 for MN: ", paste(mn_mods, collapse = ", ")),
           x = expression(-log[10]~"BH padj"), y = NULL, fill = "collection") +
      theme_bw(base_size = 9) +
      theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
            plot.title = element_text(face = "bold"))
    save_png(p2, file.path(fig_dir, "08f_glom_mn_modules"),
             width = 10, height = max(4, 1.5 * length(mn_mods) + 2))
  }
}

cat("\n08f done\n")
