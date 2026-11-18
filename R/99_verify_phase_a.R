#!/usr/bin/env Rscript
# Paper 1 - Phase A verification: assert the invariants the analysis depends on, so a regression in
# any upstream step fails loudly instead of silently corrupting a downstream figure.
# Exits 0 when every check passes, 1 otherwise.

local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[startsWith(a, "--file=")])
       source(file.path(if (length(f)) dirname(normalizePath(f[1])) else getwd(), "_lib.R")) })

fails <- 0
check <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %-58s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label, detail))
  if (!isTRUE(ok)) fails <<- fails + 1
}

tab_dir <- file.path(OUT, "tables")
set_dir <- file.path(PROJ, "work", "sets")

# ---- 1. catalogue: 178 unique MN patients, 4 releases fully redundant ---------------------------

md <- read.csv(file.path(tab_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
md$pk <- paste(md$cohort, md$patient_id)
mn <- md[md$disease == "MN", ]

check("metadata has 2165 samples (sum of all 13 accessions)", nrow(md) == 2165, nrow(md))
check("metadata has 453 MN samples", nrow(mn) == 453, nrow(mn))
check("metadata spans 13 accessions", length(unique(md$accession)) == 13, length(unique(md$accession)))
check("178 unique MN patients", length(unique(mn$pk)) == 178, length(unique(mn$pk)))

# Sample-key provenance (BAP audit 2026-11-18). `gsm` doubles as the matrix column key, and only the
# array series carry a real GEO accession: GSE197307/GSE182380 titles are local indices
# (Sample_1..n), GSE216841/GSE73953 use submitter ids. 134 local indices collide between the two
# RNA-seq series, so no downstream step may treat a non-GSM key as a participant identifier (P-2).
check("sample-key provenance recorded for every row",
      all(c("id_type", "id_is_globally_unique") %in% names(md)) && !any(is.na(md$id_type)),
      paste(names(table(md$id_type)), table(md$id_type), sep = "=", collapse = " "))
check("only real GEO accessions are marked globally unique",
      all(md$id_is_globally_unique == (md$id_type == "geo_gsm")), "")
check("local-index keys occur only in the two RNA-seq series",
      setequal(unique(md$accession[md$id_type == "local_index"]), c("GSE197307", "GSE182380")),
      paste(unique(md$accession[md$id_type == "local_index"]), collapse = ","))

sub <- read.csv(file.path(tab_dir, "release_subset_relations.csv"), stringsAsFactors = FALSE)
check("4 redundancy groups among the 13 releases",
      sum(sub$a_in_b) == 7 && length(unique(sub$accession[sub$a_in_b])) == 6,
      paste(sum(sub$a_in_b), "containment relations /", length(unique(sub$accession[sub$a_in_b])), "accessions"))

# ---- 2. every cohort set matches its declared composition ---------------------------------------

summ <- read.csv(file.path(tab_dir, "cohort_sets_summary.csv"), stringsAsFactors = FALSE)
check("7 cohort sets declared", nrow(summ) == 7, nrow(summ))

for (i in seq_len(nrow(summ))) {
  s <- summ$set[i]
  m <- readRDS(file.path(set_dir, paste0(s, "_expr.rds")))
  a <- read.csv(file.path(set_dir, paste0(s, "_samples.csv")), stringsAsFactors = FALSE)
  check(paste0(s, ": samples match summary"), ncol(m) == summ$samples[i] && nrow(a) == summ$samples[i],
        paste(ncol(m), "vs", summ$samples[i]))
  check(paste0(s, ": MN count matches summary"), sum(a$disease == "MN") == summ$MN[i],
        paste(sum(a$disease == "MN"), "vs", summ$MN[i]))
  check(paste0(s, ": one sample per patient"), !any(duplicated(paste(a$cohort, a$patient_id))),
        length(unique(paste(a$cohort, a$patient_id))))
  check(paste0(s, ": sample order aligns with matrix"), identical(a$gsm, colnames(m)))
}

# ---- 3. the leak checks the design depends on ----------------------------------------------------

disc <- read.csv(file.path(set_dir, "disc_glom_array_samples.csv"), stringsAsFactors = FALSE)
disct <- read.csv(file.path(set_dir, "disc_ti_array_samples.csv"), stringsAsFactors = FALSE)
disc_pat <- unique(c(paste(disc$cohort, disc$patient_id), paste(disct$cohort, disct$patient_id)))

for (v in c("val_glom_array", "val_ti_array")) {
  a <- read.csv(file.path(set_dir, paste0(v, "_samples.csv")), stringsAsFactors = FALSE)
  ov <- sum(paste(a$cohort, a$patient_id) %in% disc_pat)
  check(paste0(v, ": zero patient overlap with discovery (independent)"), ov == 0, ov)
}

for (c_ in c("cross_glom_rnaseq", "cross_ti_rnaseq")) {
  a <- read.csv(file.path(set_dir, paste0(c_, "_samples.csv")), stringsAsFactors = FALSE)
  # GSE197307 / GSE182380 carry no participant id (sample titles are local indices), so patient-level
  # overlap with the array cohorts cannot be computed; the sets must say so instead of reporting a
  # number produced by identifier collisions.
  check(paste0(c_, ": patient linkage marked unavailable"),
        is.na(summ$MN_overlapping_discovery[summ$set == c_]))
  check(paste0(c_, ": role labelled independent RNA-seq cohort"),
        summ$role[summ$set == c_] == "independent RNA-seq cohort")
}

check("LCM cohort is independent (0 overlap)",
      sum(paste(read.csv(file.path(set_dir, "pla2r_lcm_samples.csv"))$cohort,
                read.csv(file.path(set_dir, "pla2r_lcm_samples.csv"))$patient_id) %in% disc_pat) == 0)

# ---- 4. normalised matrices and QC artefacts exist -----------------------------------------------

rma <- list.files(file.path(PROJ, "work", "rma"), pattern = "gene_rma.*\\.rds$")
rna <- list.files(file.path(PROJ, "work", "rnaseq"), pattern = "_log2expr\\.rds$")
check("7 array cohorts normalised (GSE108109/108112 dropped by design)", length(rma) == 7, length(rma))
check("3 RNA-seq cohorts normalised", length(rna) == 3, length(rna))

pdfs <- list.files(file.path(OUT, "qc"), pattern = "^03c_.*_pca\\.pdf$")
jpegs <- list.files(file.path(OUT, "qc"), pattern = "^03c_.*_pca\\.jpeg$")
check("7 QC PCA figures as PDF", length(pdfs) == 7, length(pdfs))
check("7 QC PCA figures as JPEG", length(jpegs) == 7, length(jpegs))
check("QC summary covers every set", nrow(read.csv(file.path(tab_dir, "qc_summary.csv"))) == 7)

# ---- 5. no NA in the discovery matrices (merged U133 sets may carry NA by design) ----------------

for (s in c("disc_glom_array", "disc_ti_array", "cross_glom_rnaseq", "cross_ti_rnaseq", "pla2r_lcm")) {
  m <- readRDS(file.path(set_dir, paste0(s, "_expr.rds")))
  check(paste0(s, ": no missing values"), !anyNA(m))
}

cat(sprintf("\n%s - %d check(s) failed\n", if (fails == 0) "ALL CHECKS PASSED" else "FAILURES", fails))
quit(status = if (fails == 0) 0L else 1L)
