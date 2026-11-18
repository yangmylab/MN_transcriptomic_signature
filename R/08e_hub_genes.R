#!/usr/bin/env Rscript
# Paper 1 - Phase D, step 08e: hub genes of preserved modules, and the Phase E candidate space.
#
# Protocol: hub genes are |kME| > 0.8 inside modules that PASSED the ERCB preservation test.
# A hub of an unpreserved module is a hub of a discovery-cohort artefact, so it is excluded here
# rather than quietly carried into Phase E.
#
# The candidate space for Phase E is the union of four independent evidence lines, and every gene
# keeps its provenance so the panel can be justified line by line:
#   1. hub           |kME| > 0.8 in a preserved module
#   2. meta-DEG      member of the random-effects meta-DEGs for this compartment (any contrast)
#   3. gate feature  selected by the elastic-net separability gate at some K
#   4. celltype      bulk DEG attributed to a cell type with tau > 0.5 (specific, not ubiquitous)
#
# Outputs
#   output/tables/wgcna_hub_genes.csv              every hub gene with its module and preservation Z
#   output/tables/wgcna_hub_gene_overlap.csv       how many genes sit in which evidence intersection
#   output/tables/phaseE_candidate_space.csv       ranked candidate space, all evidence lines
#   output/tables/phaseE_candidate_core.csv        genes with >= 2 evidence lines (Phase E input)
#   output/wgcna/08e_<comp>_hub_network.{pdf,jpeg} top-3 preserved modules, top-25 hubs by |kME|
#   output/wgcna/08e_<comp>_hub_kme.{pdf,jpeg}     kME distribution inside preserved modules

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

suppressPackageStartupMessages({ library(WGCNA); library(igraph) })

tab_dir <- file.path(OUT, "tables")
set_dir <- file.path(PROJ, "work", "sets")
wk_dir  <- file.path(PROJ, "work", "wgcna")
fig_dir <- file.path(OUT, "wgcna")

KME_CUT  <- 0.8
TAU_CUT  <- 0.5
# the pDC (97 atlas cells) and Neutrophil (535 cells) columns carry ambient signal and are
# unreliable (handoff B8); they are excluded from the trait set and equally from the
# cell-type-attribution evidence line
UNRELIABLE_CELLTYPES <- c("Neutrophil", "pDC")
TOP_N_NET <- 25L

fig <- function(stem, width, height, draw) {
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()
  jpeg(paste0(stem, ".jpeg"), width = width, height = height, units = "in", res = 300,
       quality = 95, bg = "white")
  draw()
  dev.off()
}

write_tab <- function(x, name) {
  f <- file.path(tab_dir, name)
  write.csv(x, f, row.names = FALSE)
  cat("  wrote", f, "(", nrow(x), "rows )\n")
}

attr_tab <- read.csv(file.path(tab_dir, "bulk_deg_celltype_attribution.csv"), stringsAsFactors = FALSE)
gate_tab <- read.csv(file.path(tab_dir, "separability_selected_genes.csv"), stringsAsFactors = FALSE)
gate_genes <- unique(gate_tab$gene)

meta_files <- list.files(file.path(tab_dir, "meta"), pattern = "_metaDEG\\.csv$", full.names = TRUE)
meta_all <- do.call(rbind, lapply(meta_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d[, c("compartment", "contrast", "gene", "logFC", "p.adj", "cohorts")]
}))
meta_all$compartment <- ifelse(meta_all$compartment == "Glom", "glom", "ti")
meta_all <- meta_all[!duplicated(meta_all[, c("compartment", "gene", "contrast")]), ]

cat("meta-DEG rows:", nrow(meta_all), "| contrasts:", length(unique(meta_all$contrast)),
    "| gate genes:", length(gate_genes), "| attributed genes:", length(unique(attr_tab$gene)), "\n")

hubs_all <- list()
cand_all <- list()

for (comp in c("glom", "ti")) {
  cat("\n================ ", toupper(comp), " ================\n", sep = "")

  mods <- read.csv(file.path(tab_dir, paste0("wgcna_modules_", comp, ".csv")), stringsAsFactors = FALSE)
  pres <- read.csv(file.path(tab_dir, paste0("wgcna_preservation_", comp, ".csv")), stringsAsFactors = FALSE)
  preserve_set <- pres$module[pres$verdict == "preserved"]
  weak_set     <- pres$module[pres$verdict == "weak"]
  cat("preserved modules:", if (length(preserve_set)) paste(preserve_set, collapse = ",") else "none",
      "| weak:", if (length(weak_set)) paste(weak_set, collapse = ",") else "none", "\n")

  h <- mods[!is.na(mods$kME) & abs(mods$kME) > KME_CUT & mods$module %in% preserve_set, ]
  h <- h[order(-abs(h$kME)), ]
  h$compartment <- comp
  h$module_preservation_Z <- pres$Zsummary[match(h$module, pres$module)]
  h$module_preservation_q <- pres$qsummary[match(h$module, pres$module)]
  h$module_verdict <- "preserved"
  hubs_all[[comp]] <- h
  cat("hub genes in preserved modules:", nrow(h), "across", length(unique(h$module)), "modules\n")

  # module-level summary, including unpreserved modules so nothing is hidden
  ms <- do.call(rbind, lapply(split(mods, mods$module), function(d) {
    data.frame(compartment = comp, module = d$module[1], n_genes = nrow(d),
               n_hub = sum(abs(d$kME) > KME_CUT, na.rm = TRUE),
               median_abs_kME = round(median(abs(d$kME), na.rm = TRUE), 4),
               max_abs_kME = round(max(abs(d$kME), na.rm = TRUE), 4),
               preservation_Z = pres$Zsummary[match(d$module[1], pres$module)],
               preservation_verdict = pres$verdict[match(d$module[1], pres$module)],
               stringsAsFactors = FALSE)
  }))
  ms <- ms[order(-ms$n_genes), ]

  # ---- candidate space ----------------------------------------------------
  md <- meta_all[meta_all$compartment == comp, ]
  md_best <- do.call(rbind, lapply(split(md, md$gene), function(d) {
    d[which.min(d$p.adj), c("gene", "contrast", "logFC", "p.adj"), drop = FALSE][1, ]
  }))
  names(md_best)[2:4] <- c("meta_contrast", "meta_logFC", "meta_padj")

  # attribution rows are keyed "meta_<Glom|TI>_<contrast>", the meta tables carry the bare
  # contrast - build the full key or the join silently returns nothing
  md$set_key <- paste0("meta_", ifelse(comp == "glom", "Glom", "TI"), "_", md$contrast)
  at <- attr_tab[attr_tab$set %in% unique(md$set_key), ]
  at <- at[!duplicated(at$gene), c("gene", "top_celltype", "tau", "second_celltype")]
  at$reliable <- !(at$top_celltype %in% UNRELIABLE_CELLTYPES)
  cat("  attributed genes:", nrow(at), "| reliable cell type:", sum(at$reliable),
      "| of those tau >", TAU_CUT, ":", sum(at$reliable & at$tau > TAU_CUT), "\n")

  # B9: the gate's K = "all" row is the full 13,379-gene elastic-net model, not a marker panel, so
  # a gene selected only there does not count as a gate evidence line
  gate_reduced <- gate_tab$gene[gate_tab$K != "all"]

  genes <- Reduce(union, list(h$gene, md_best$gene, gate_reduced,
                              at$gene[at$reliable & at$tau > TAU_CUT]))
  cand <- data.frame(gene = genes, compartment = comp, stringsAsFactors = FALSE)
  cand$is_hub         <- cand$gene %in% h$gene
  cand$hub_module     <- h$module[match(cand$gene, h$gene)]
  cand$hub_kME        <- round(h$kME[match(cand$gene, h$gene)], 4)
  cand$module_Z       <- h$module_preservation_Z[match(cand$gene, h$gene)]
  cand$is_metaDEG     <- cand$gene %in% md_best$gene
  cand$meta_contrast  <- md_best$meta_contrast[match(cand$gene, md_best$gene)]
  cand$meta_logFC     <- round(md_best$meta_logFC[match(cand$gene, md_best$gene)], 4)
  cand$meta_padj      <- signif(md_best$meta_padj[match(cand$gene, md_best$gene)], 4)
  cand$in_gate        <- cand$gene %in% gate_genes
  cand$gate_min_K     <- vapply(cand$gene, function(g) {
    k <- gate_tab$K[gate_tab$gene == g]
    if (!length(k)) return(NA_character_)
    k <- k[k != "all"]
    if (!length(k)) "all" else as.character(min(as.integer(k)))
  }, character(1))
  cand$top_celltype   <- at$top_celltype[match(cand$gene, at$gene)]
  cand$celltype_tau   <- round(at$tau[match(cand$gene, at$gene)], 4)
  cand$celltype_reliable <- at$reliable[match(cand$gene, at$gene)]
  cand$is_specific    <- !is.na(cand$celltype_tau) & cand$celltype_tau > TAU_CUT &
                          !is.na(cand$celltype_reliable) & cand$celltype_reliable
  cand$in_gate_reduced <- cand$gene %in% gate_reduced
  cand$n_evidence <- cand$is_hub + cand$is_metaDEG + cand$in_gate_reduced + cand$is_specific
  cand <- cand[order(-cand$n_evidence, -abs(cand$hub_kME), cand$meta_padj), ]
  cand_all[[comp]] <- cand

  ov <- as.data.frame(table(n_evidence = cand$n_evidence), stringsAsFactors = FALSE)
  names(ov) <- c("n_evidence", "n_genes")
  ov$compartment <- comp
  cat("candidate genes:", nrow(cand), "| with >=2 evidence lines:", sum(cand$n_evidence >= 2),
      "| all four:", sum(cand$n_evidence == 4), "\n")

  # ---- figures ------------------------------------------------------------
  pres_ok <- pres[pres$verdict == "preserved", ]
  pres_ok <- pres_ok[order(-pres_ok$Zsummary), ]
  net_mods <- head(pres_ok$module, 3)
  dat <- readRDS(file.path(wk_dir, paste0("datExpr_", comp, ".rds")))
  pdf_ok <- length(net_mods) > 0
  if (pdf_ok) {
    fig(file.path(fig_dir, paste0("08e_", comp, "_hub_network")), 11, 4.2, function() {
      par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
      for (m in net_mods) {
        top <- head(h$gene[h$module == m], TOP_N_NET)
        if (length(top) < 3) { plot.new(); title(main = paste0(m, ": <3 hubs")); next }
        cm <- cor(dat[, top, drop = FALSE], use = "p")
        g <- graph_from_adjacency_matrix(abs(cm), mode = "undirected", weighted = TRUE, diag = FALSE)
        E(g)$weight <- abs(cm)[lower.tri(abs(cm))]
        g <- delete_edges(g, E(g)[weight < quantile(E(g)$weight, 0.4)])
        V(g)$size <- 4 + 10 * abs(h$kME[match(V(g)$name, h$gene)])
        cols <- ifelse(cor(dat[, V(g)$name, drop = FALSE], dat[, top[1], drop = FALSE])[, 1] > 0,
                       "firebrick", "steelblue")
        set.seed(20261019)
        plot(g, layout = layout_with_fr(g), vertex.color = cols, vertex.frame.color = "white",
             vertex.label = V(g)$name, vertex.label.cex = 0.55, vertex.label.color = "grey15",
             edge.color = "grey80", main = paste0(m, " (Z=", round(pres$Zsummary[pres$module == m], 1),
                                                  ", ", nrow(h[h$module == m, ]), " hubs)"))
      }
    })
  }

  fig(file.path(fig_dir, paste0("08e_", comp, "_hub_kme")), 7.5, 4.6, function() {
    par(mar = c(4.5, 4.5, 3, 1))
    keep <- mods[mods$module %in% unique(c(preserve_set, weak_set)), ]
    boxplot(abs(kME) ~ module, data = keep[order(-abs(keep$kME)), ],
            col = ifelse(levels(factor(keep$module)) %in% preserve_set, "seagreen", "orange"),
            border = "grey30", las = 2, ylab = "|kME|", xlab = NULL,
            main = paste0("Module membership - ", comp, " (green = preserved, orange = weak)"))
    abline(h = KME_CUT, lty = 2, col = "firebrick")
  })

  hubs_all[[comp]] <- h
}

hubs <- do.call(rbind, hubs_all)
write_tab(hubs[, c("compartment", "gene", "module", "module_label", "kME", "p.kME",
                   "module_preservation_Z", "module_preservation_q", "module_verdict")],
          "wgcna_hub_genes.csv")

cand <- do.call(rbind, cand_all)
cand <- cand[order(cand$compartment, -cand$n_evidence, -abs(cand$hub_kME), cand$meta_padj), ]
write_tab(cand, "phaseE_candidate_space.csv")

# the actual Phase E starting point: genes supported by at least two independent evidence lines
core <- cand[cand$n_evidence >= 2, ]
write_tab(core, "phaseE_candidate_core.csv")
cat("core candidates (>=2 evidence lines):", nrow(core), "\n")

ov_all <- do.call(rbind, lapply(cand_all, function(x) {
  o <- as.data.frame(table(n_evidence = x$n_evidence), stringsAsFactors = FALSE)
  names(o) <- c("n_evidence", "n_genes")
  data.frame(compartment = x$compartment[1], o, stringsAsFactors = FALSE)
}))
write_tab(ov_all, "wgcna_hub_gene_overlap.csv")

cat("\n--- hubs by module ---\n")
print(as.data.frame(table(hubs$compartment, hubs$module)))
cat("\n08e done\n")
