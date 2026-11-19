#!/usr/bin/env python3
r"""Compose the manuscript figures as VECTOR PDFs from the saved single-panel vector PDFs.

Each figure is a LaTeX grid of \includegraphics{panel.pdf} calls, so every source panel stays
vector end-to-end (no rasterisation), then pdftoppm produces a JPEG preview at 300 dpi.

Usage:  python3 tools/make_figures.py [--only fig4]
Output: output/figures/<figure>.pdf and .jpeg  (+ output/figures/supplementary/)
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_cohort_design import build_design_panels

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "output" / "figures"
SUPP = OUT / "supplementary"
TEXDIR = PROJ / "work" / "figtex"

O = PROJ / "output"

# The Phase E panel size is data-driven (09a's consensus rule), so figure wording must follow the
# locked panel rather than a hard-coded "five-gene" — it became 4 genes after the 2026-10-20 audit.
try:
    with open(O / "tables" / "09a_panel_consensus.csv", encoding="utf-8") as fh:
        N_PANEL = max(0, sum(1 for line in fh.read().splitlines()[1:] if line.strip()))
except OSError:
    N_PANEL = 0
PANEL_WORD = {3: "three", 4: "four", 5: "five", 6: "six"}.get(N_PANEL, str(N_PANEL))

# (panel pdf, sub-label, width fraction of the text block)
FIGURES = {
    "Fig1_study_design_and_QC": {
        "label": "Figure 1.",
        "title": "Cohort design and analytical independence",
        "cols": 2,
        "panels": [
            (O / "qc/Fig1_study_flow.pdf", "A", 1.0),
            (O / "qc/Fig1_patient_deduplication.pdf", "B", 1.0),
            (O / "qc/Fig1_cohort_allocation.pdf", "C", 1.0),
            (O / "qc/04a_disc_ti_array_pca_before_after.pdf", "D", 1.0),
            (O / "qc/Fig1_combat_quantitative.pdf", "E", 1.0),
        ],
    },
    "Fig2_discovery_and_meta": {
        "label": "Figure 2.",
        "title": "Differential expression and random-effects meta-analysis",
        "cols": 2,
        "panels": [
            (O / "deg/04c_Glom_MN_vs_Control_volcano.pdf", "A", 1.0),
            (O / "deg/04c_Glom_MN_vs_Control_forest_CCNL1.pdf", "B", 1.0),
            (O / "deg/04c_Glom_MN_vs_Control_forest_HBA2.pdf", "C", 1.0),
            (O / "deg/04c_Glom_MN_vs_FSGS_volcano.pdf", "D", 1.0),
            (O / "deg/04c_Glom_MN_vs_nonMN_volcano.pdf", "E", 1.0),
            (O / "deg/04c_Glom_metaDEG_heatmap.pdf", "F", 1.0),
        ],
    },
    "Fig3_separability_gate": {
        "label": "Figure 3.",
        "title": "Separability of MN from its mimics, and why a DEG count hides it",
        "cols": 2,
        "panels": [
            (O / "models_gate/05a_gate_roc.pdf", "A", 1.0),
            (O / "models_gate/05a_gate_auc_vs_k.pdf", "B", 1.0),
            (O / "models_gate/05b_gate_permutation_null.pdf", "C", 1.0),
            (O / "models_gate/05b_gate_auc_by_comparator.pdf", "D", 1.0),
            (O / "models_gate/05a_gate_confusion.pdf", "E", 1.0),
        ],
    },
    "Fig4_WGCNA_modules": {
        "label": "Figure 4.",
        "title": "Glomerular co-expression modules, trait association and ERCB preservation",
        "cols": 2,
        "panels": [
            (O / "wgcna/08a_glom_module_dendro.pdf", "A", 1.0),
            (O / "wgcna/08b_glom_module_trait.pdf", "B", 1.0),
            (O / "wgcna/08b_glom_me_box.pdf", "C", 1.0),
            (O / "wgcna/08d_glom_preservation.pdf", "D", 1.0),
            (O / "wgcna/08f_glom_mn_modules.pdf", "E", 1.0),
        ],
    },
    "Fig5_biomarker_panel": {
        "label": "Figure 5.",
        "title": f"A {PANEL_WORD}-gene glomerular panel for MN",
        "cols": 2,
        "panels": [
            (O / "models_panel/09a_frequency.pdf", "A", 1.0),
            (O / "models_panel/09b_roc.pdf", "B", 1.0),
            (O / "models_panel/09b_forest.pdf", "C", 1.0),
            (O / "models_panel/09b_score_box.pdf", "D", 1.0),
            (O / "models_panel/09e_by_comparator.pdf", "E", 1.0),
        ],
    },
    "Fig6_atlas_and_deconvolution": {
        "label": "Figure 6.",
        "title": "Kidney single-cell atlas and bulk deconvolution",
        "cols": 2,
        "panels": [
            (O / "scrna/07c_umap_celltype.pdf", "A", 1.0),
            (O / "scrna/07d_composition.pdf", "B", 1.0),
            (O / "scrna/07f_deconvolution_stacked.pdf", "C", 1.0),
            (O / "scrna/07f_deconvolution_boxplots.pdf", "D", 1.0),
            (O / "scrna/07e_bulk_vs_pseudobulk.pdf", "E", 1.0),
        ],
    },
    # ---------------- supplementary ----------------
    "FigS1_bulk_cohort_qc": {
        "label": "Figure S1.",
        "title": "Bulk cohort expression quality control",
        "cols": 2,
        "panels": [
            (O / "qc/03c_disc_glom_array_pca.pdf", "A", 1.0),
            (O / "qc/03c_val_glom_array_pca.pdf", "B", 1.0),
            (O / "qc/03c_cross_glom_rnaseq_pca.pdf", "C", 1.0),
            (O / "qc/03c_pla2r_lcm_pca.pdf", "D", 1.0),
        ],
    },
    "FigS2_cohort_and_singlecell_qc": {
        "label": "Figure S2.",
        "title": "Cohort QC, gate stability and single-cell annotation",
        "cols": 2,
        "panels": [
            (O / "qc/03c_val_ti_array_pca.pdf", "A", 1.0),
            (O / "qc/03c_cross_ti_rnaseq_pca.pdf", "B", 1.0),
            (O / "models_gate/05c_gate_seed_stability.pdf", "C", 1.0),
            (O / "scrna/07c_marker_dotplot.pdf", "D", 1.0),
            (O / "scrna/07b_sample_qc.pdf", "E", 1.0),
        ],
    },
    "FigS3_WGCNA_stability_and_hubs": {
        "label": "Figure S3.",
        "title": "Module reproducibility, hub genes and annotation detail",
        "cols": 2,
        "panels": [
            (O / "wgcna/08a_glom_softpower.pdf", "A", 1.0),
            (O / "wgcna/08c_glom_stability.pdf", "B", 1.0),
            (O / "wgcna/08e_glom_hub_network.pdf", "C", 1.0),
            (O / "wgcna/08e_glom_hub_kme.pdf", "D", 1.0),
            (O / "wgcna/08f_glom_module_dotplot.pdf", "E", 1.0),
        ],
    },
    "FigS4_functional_profiling": {
        "label": "Figure S4.",
        "title": "Functional profiling of the glomerular MN signature",
        "cols": 2,
        "panels": [
            (O / "enrichment/04d_Glom_MN_vs_Control_GO_BP_dotplot.pdf", "A", 1.0),
            (O / "enrichment/04d_Glom_MN_vs_Control_GSEA_C7_bar.pdf", "B", 1.0),
            (O / "enrichment/04d_Glom_MN_vs_Control_GSEA_HALLMARK_APICAL_JUNCTION.pdf", "C", 1.0),
            (O / "enrichment/04d_Glom_MN_vs_Control_GSEA_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf", "D", 1.0),
            (O / "enrichment/04d_Glom_MN_vs_nonMN_GO_BP_dotplot.pdf", "E", 1.0),
        ],
    },
    "FigS5_panel_robustness": {
        "label": "Figure S5.",
        "title": "Panel robustness, calibration and transfer",
        "cols": 2,
        "panels": [
            (O / "models_panel/09c_robustness_roc.pdf", "A", 1.0),
            (O / "models_panel/09b_calibration.pdf", "B", 1.0),
            (O / "models_panel/09c_recalibration.pdf", "C", 1.0),
            (O / "models_panel/09b_dca.pdf", "D", 1.0),
            (O / "models_panel/09d_leave_one_out.pdf", "E", 1.0),
        ],
    },
    "FigS6_localisation_and_pla2r": {
        "label": "Figure S6.",
        "title": "Where the panel genes sit, panel sensitivity and PLA2R stratum",
        "cols": 2,
        "panels": [
            (O / "models_panel/09d_gene_celltype.pdf", "A", 1.0),
            (O / "models_panel/09d_pseudobulk.pdf", "B", 1.0),
            (O / "models_panel/09e_no_slamf8_roc.pdf", "C", 1.0),
            (O / "models_panel/09e_authors_check.pdf", "D", 1.0),
            (O / "models_panel/09f_pla2r_strata.pdf", "E", 1.0),
        ],
    },
    "FigS7_subtypes_negative": {
        "label": "Figure S7.",
        "title": "Exploratory molecular subtyping is not reproducible",
        "cols": 2,
        "panels": [
            (O / "subtypes/06a_array_consensus_k2.pdf", "A", 1.0),
            (O / "subtypes/06a_stability_pac.pdf", "B", 1.0),
            (O / "subtypes/06a_centroid_concordance.pdf", "C", 1.0),
            (O / "subtypes/06a_array_subtype_heatmap.pdf", "D", 1.0),
            (O / "subtypes/06a_rnaseq_consensus_k3.pdf", "E", 1.0),
        ],
    },
}

PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[paperwidth=190mm,paperheight=%(height).0fmm,margin=8mm]{geometry}
\usepackage{graphicx}
\usepackage{helvet}
\renewcommand{\familydefault}{\sfdefault}
\usepackage{xcolor}
\usepackage{caption}
\pagestyle{empty}
\setlength{\parindent}{0pt}
\newcommand{\panel}[4]{%% #1 file, #2 label, #3 width fraction, #4 height in mm
  \begin{minipage}[t]{#3\linewidth}\centering
    {\bfseries\small (#2)}\par\vspace{1mm}
    \includegraphics[width=\linewidth,height=#4mm,keepaspectratio]{#1}
  \end{minipage}}
\begin{document}
"""

# Panel height is a fixed physical size rather than 0.40\textheight: with a content-driven page
# height the two would be circular (page height -> textheight -> panel height -> page height), and
# the fixed 250 mm page left a 30-35 % empty band at the bottom of every plate (audit F37).
PANEL_H_MM = 82


def _ink_height_mm(pdf_path):
    """Height of the drawn content in mm, measured with Ghostscript's bbox device.

    Ghostscript cannot read this exfat workspace (P-49), so the file is staged through /tmp.
    Returns None when the measurement fails, in which case the caller keeps its estimate.
    """
    stage = Path("/tmp") / f"_bbox_{pdf_path.name}"
    shutil.copyfile(pdf_path, stage)
    r = subprocess.run(["gs", "-dNOPAUSE", "-dBATCH", "-sDEVICE=bbox", str(stage)],
                       capture_output=True, text=True)
    stage.unlink(missing_ok=True)
    ys = []
    for line in (r.stdout + r.stderr).splitlines():
        if "HiResBoundingBox" in line:
            x0, y0, x1, y1 = (float(v) for v in line.split(":")[1].split())
            ys.append((y0, y1))
    if not ys:
        return None
    return (max(y1 for _, y1 in ys) - min(y0 for y0, _ in ys)) * 25.4 / 72.0


def _page_count(pdf_path):
    """Number of pages in a PDF (a plate must be exactly one)."""
    r = subprocess.run(["pdfinfo", str(pdf_path)], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split()[1])
    return 1


def build(name, spec, dest):
    missing = [str(p) for p, _, _ in spec["panels"] if not p.exists()]
    if missing:
        print(f"  SKIP {name}: missing {len(missing)} panel(s)")
        for m in missing:
            print("      ", m)
        return False

    cols = spec["cols"]
    # Multi-panel plates (>= 5 panels) use 62 mm panels so that three rows still fit an A4 page;
    # the old 82 mm is kept for compact 4-panel plates.
    ph = spec.get("panel_h", 62 if len(spec["panels"]) >= 5 else PANEL_H_MM)
    frac = 0.98 / cols
    rows = (len(spec["panels"]) + cols - 1) // cols
    # generous first guess, then shrink to the measured content height + 2 x 8 mm margins
    height_mm = 16 + rows * (ph + 16)

    def write_tex(page_mm):
        body = [PREAMBLE % {"height": page_mm},
                r"{\bfseries\large %s %s.}\par\vspace{2mm}" % (spec["label"], spec["title"]), "\n"]
        for i, (path, label, w) in enumerate(spec["panels"]):
            body.append(r"\panel{%s}{%s}{%.3f}{%d}" % (path, label, frac * w, ph))
            body.append("\n\\vspace{3mm}\n")
            if (i + 1) % cols == 0 and i != len(spec["panels"]) - 1:
                body.append("\n\\par\\vspace{2mm}\n")
        body.append("\n\\end{document}\n")
        TEXDIR.mkdir(parents=True, exist_ok=True)
        dest.mkdir(parents=True, exist_ok=True)
        tex = TEXDIR / f"{name}.tex"
        tex.write_text("".join(body))
        return tex

    for _ in range(3):
        tex = write_tex(height_mm)
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            f"-output-directory={TEXDIR}", str(tex)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  FAILED {name} (pdflatex)")
            print(r.stdout[-1500:])
            return False
        ink = _ink_height_mm(TEXDIR / f"{name}.pdf")
        if ink is None:
            break
        target = round(ink + 16)  # 8 mm top + 8 mm bottom margin
        if abs(target - height_mm) <= 2:
            break
        height_mm = target

    pdf_src = TEXDIR / f"{name}.pdf"

    # A plate must be exactly ONE page. The bbox-derived height can undershoot by a few mm, which
    # pushes the last row onto a second (blank-ish) page: journals reject that, and pdftoppm leaves
    # a stray "-2.jpg" evidence file behind. Grow the page until LaTeX stops breaking.
    guard = 0
    base_mm = 16 + rows * (ph + 16)
    height_mm = max(height_mm, base_mm)   # the bbox loop shrinks; a multi-page plate needs the row budget
    while _page_count(pdf_src) > 1 and guard < 6:
        guard += 1
        height_mm += 10
        print(f"  {name}: {_page_count(pdf_src)} pages -> retry at {height_mm} mm")
        tex = write_tex(height_mm)
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            f"-output-directory={TEXDIR}", str(tex)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  FAILED {name} (pdflatex, page-height guard)")
            print(r.stdout[-1500:])
            return False
    n_pages = _page_count(pdf_src)
    if n_pages > 1:
        print(f"  WARN {name}: still {n_pages} pages after {guard} height increases")
    pdf_dst = dest / f"{name}.pdf"

    # Submission-quality pass: R's default pdf() device leaves Helvetica un-embedded, and journals
    # reject that. Ghostscript rewrites the composite with every font embedded and subset.
    # ghostscript on this box cannot read the exfat workspace, so the pass is staged via /tmp
    stage_in = Path("/tmp") / f"_fig_in_{name}.pdf"
    stage_out = Path("/tmp") / f"_fig_out_{name}.pdf"
    shutil.copyfile(pdf_src, stage_in)
    gs = subprocess.run(
        ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dNOSAFER", "-sDEVICE=pdfwrite",
         "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true", "-dSubsetFonts=true",
         "-dCompatibilityLevel=1.4", f"-sOutputFile={stage_out}", str(stage_in)],
        capture_output=True, text=True)
    pdf_dst.write_bytes(stage_out.read_bytes() if gs.returncode == 0 and stage_out.exists()
                        else pdf_src.read_bytes())
    if gs.returncode != 0:
        print(f"  note: ghostscript embedding failed for {name}; using the pdflatex output")
    for f in (stage_in, stage_out):
        f.unlink(missing_ok=True)
    subprocess.run(["pdftoppm", "-jpeg", "-r", "300", "-jpegopt", "quality=92",
                    str(pdf_dst), str(dest / name)], check=True)
    jpg = dest / f"{name}-1.jpg"
    final_jpg = dest / f"{name}.jpeg"
    if jpg.exists():
        jpg.replace(final_jpg)
    for extra in dest.glob(f"{name}-*.jpg"):   # never leave pdftoppm page-preview leftovers
        extra.unlink(missing_ok=True)
    print(f"  {name}: {len(spec['panels'])} panels, {_page_count(pdf_dst)} page -> PDF + JPEG")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    if args.only is None or "fig1" in args.only.lower():
        build_design_panels()

    ok = 0
    total = 0
    for name, spec in FIGURES.items():
        if args.only and args.only.lower() not in name.lower():
            continue
        total += 1
        dest = OUT if name.startswith("Fig") and not name.startswith("FigS") else SUPP
        print(f"[{name}]")
        if build(name, spec, dest):
            ok += 1
    print(f"\n{ok}/{total} figures built -> {OUT}")
    return 0 if ok == total else 1


if __name__ == "__main__":
    sys.exit(main())
