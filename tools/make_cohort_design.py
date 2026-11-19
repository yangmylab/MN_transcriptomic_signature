#!/usr/bin/env python3
"""Build the Figure 1 design panels from saved cohort and QC result tables."""
import csv
import shutil
import subprocess
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "output"
TEXDIR = PROJ / "work" / "figtex" / "cohort_design"
DEST = OUT / "qc"


def rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def esc(value):
    return (str(value).replace("\\", r"\textbackslash{}").replace("&", r"\&")
            .replace("%", r"\%").replace("_", r"\_").replace("#", r"\#"))


def compile_panel(name, body, height_mm=115):
    TEXDIR.mkdir(parents=True, exist_ok=True)
    DEST.mkdir(parents=True, exist_ok=True)
    tex = TEXDIR / f"{name}.tex"
    tex.write_text(r"""\documentclass[10pt]{article}
\usepackage[paperwidth=170mm,paperheight=""" + str(height_mm) + r"""mm,margin=5mm]{geometry}
\usepackage{tikz,xcolor,helvet}
\usetikzlibrary{arrows.meta,positioning,calc}
\definecolor{navy}{HTML}{1F4E79}
\definecolor{teal}{HTML}{2A9D8F}
\definecolor{orange}{HTML}{E76F51}
\definecolor{lightblue}{HTML}{DCE6F1}
\definecolor{lightgrey}{HTML}{E9ECEF}
\definecolor{midgrey}{HTML}{6C757D}
\pagestyle{empty}
\setlength{\parindent}{0pt}
\begin{document}
\sffamily\fontsize{10}{12}\selectfont
\hyphenpenalty=10000\exhyphenpenalty=10000
""" + body + r"""
\end{document}
""", encoding="utf-8")
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                        f"-output-directory={TEXDIR}", str(tex)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"{name}: pdflatex failed\n{r.stdout[-3000:]}")
    stage_in = Path("/tmp") / f"_{name}_in.pdf"
    stage_out = Path("/tmp") / f"_{name}_out.pdf"
    shutil.copyfile(TEXDIR / f"{name}.pdf", stage_in)
    gs = subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dNOSAFER", "-sDEVICE=pdfwrite",
                         "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true", "-dSubsetFonts=true",
                         "-dCompatibilityLevel=1.4", f"-sOutputFile={stage_out}", str(stage_in)],
                        capture_output=True, text=True)
    if gs.returncode:
        raise RuntimeError(f"{name}: Ghostscript failed: {gs.stdout}{gs.stderr}")
    pdf = DEST / f"{name}.pdf"
    pdf.write_bytes(stage_out.read_bytes())
    for stage in (stage_in, stage_out):
        stage.unlink(missing_ok=True)
    subprocess.run(["pdftoppm", "-f", "1", "-singlefile", "-jpeg", "-r", "300",
                    "-jpegopt", "quality=95", str(pdf), str(DEST / name)], check=True)
    (DEST / f"{name}.jpg").replace(DEST / f"{name}.jpeg")
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True, check=True).stdout
    if "Pages:           1" not in info:
        raise RuntimeError(f"{name}: expected one page\n{info}")
    print(f"  {name}: PDF + JPEG")


def build_design_panels():
    meta = rows(OUT / "tables" / "sample_metadata.csv")
    patients = rows(OUT / "tables" / "patient_map.csv")
    cohorts = rows(OUT / "tables" / "cohort_sets_summary.csv")
    dedup = rows(OUT / "tables" / "dedup_summary.csv")
    relations = rows(OUT / "tables" / "release_subset_relations.csv")
    combat = rows(OUT / "tables" / "combat_summary.csv")[0]

    n_accessions = len({r["accession"] for r in meta})
    n_samples = len(meta)
    n_mn = sum(r["disease"] == "MN" for r in meta)
    n_patients = len(patients)
    n_sets = len(cohorts)
    contained = [r for r in relations if r["a_in_b"] == "TRUE"]
    n_groups = 4

    def node(x, y, top, bottom, colour="navy"):
        return (f"\\node[draw={colour}, rounded corners=2pt, very thick, fill=lightblue, "
                f"minimum width=3.65cm, minimum height=1.8cm, text width=3.3cm, align=center] at ({x},{y}) "
                f"{{\\textbf{{{esc(top)}}}\\\\[2pt]{esc(bottom)}}};\n")

    body = ["\\begin{tikzpicture}[x=1cm,y=1cm]\n",
            "\\node[anchor=west,text width=15.5cm,font=\\bfseries\\fontsize{13}{15}\\selectfont] at (0,9.7) {Cohort assembly and patient-level de-duplication};\n"]
    xs = [1.75, 5.7, 9.65, 13.6]
    tops = [f"{n_accessions} GEO accessions", f"{n_samples:,} samples",
            f"{n_mn} MN samples / {n_patients} patients", f"{n_sets} non-redundant sets"]
    bottoms = ["public transcriptomic releases", "all compartments and diagnoses",
               f"{n_groups} redundancy groups resolved", "discovery and independent validation"]
    for i, x in enumerate(xs):
        body.append(node(x, 7.65, tops[i], bottoms[i], "teal" if i == 3 else "navy"))
        if i < 3:
            body.append(f"\\draw[-{{Stealth[length=3mm]}},very thick,midgrey] ({x + 1.9},7.65) -- ({xs[i + 1] - 1.9},7.65);\n")
    roles = [
        ("Discovery array", "153 samples / 51 MN"),
        ("ERCB external", "196 samples / 21 MN"),
        ("RNA-seq", "282 samples / 62 MN"),
        ("Laser capture", "34 samples / 12 MN"),
    ]
    for i, (role, count) in enumerate(roles):
        x = xs[i]
        body.append(f"\\node[draw=lightblue,fill=lightgrey,rounded corners=2pt,minimum width=3.65cm,minimum height=1.35cm,text width=3.3cm,align=center] at ({x},4.65) {{\\textbf{{{esc(role)}}}\\\\{esc(count)}}};\n")
    body.extend([
        "\\node[anchor=west,text=navy] at (0.2,2.45) {\\textbf{Independence safeguards}};\n",
        "\\node[anchor=west,text width=15.5cm] at (0.2,1.45) {ERCB: 0 shared participants with discovery \\hspace{5mm} Laser capture: 0 overlap};\n",
        "\\node[anchor=west,text width=15.5cm] at (0.2,0.45) {RNA-seq: participant linkage unavailable; analysed as an independent cohort, never inferred from sample titles.};\n",
        "\\end{tikzpicture}\n",
    ])
    compile_panel("Fig1_study_flow", "".join(body))

    family = {r["cohort"]: r for r in dedup if r["cohort"] in ("NEPTUNE-numbered", "ERCB-numbered")}
    body = ["\\begin{tikzpicture}[x=1cm,y=1cm]\n",
            "\\node[anchor=west,text width=15.5cm,font=\\bfseries\\fontsize{13}{15}\\selectfont] at (0,9.7) {Patient redundancy and release containment};\n",
            "\\node[anchor=west,font=\\bfseries] at (0.2,8.25) {MN samples per unique patient};\n"]
    for i, label in enumerate(("NEPTUNE-numbered", "ERCB-numbered")):
        r = family[label]
        y = 6.75 - i * 1.65
        redundancy = float(r["redundancy"])
        width = redundancy * 1.85
        short = "NEPTUNE" if label.startswith("NEPTUNE") else "ERCB"
        body.append(f"\\node[anchor=east] at (2.15,{y + 0.32}) {{{short}}};\n")
        body.append(f"\\fill[lightblue] (2.35,{y}) rectangle ({2.35 + width:.3f},{y + 0.65});\n")
        body.append(f"\\draw[navy,very thick] (2.35,{y}) rectangle ({2.35 + width:.3f},{y + 0.65});\n")
        body.append(f"\\node[anchor=west,text=navy,font=\\bfseries] at ({2.55 + width:.3f},{y + 0.32}) {{{redundancy:.2f}$\\times$}};\n")
        body.append(f"\\node[anchor=west,text=midgrey] at (2.35,{y - 0.42}) {{{r['MN_samples']} samples / {r['unique_MN_patients']} patients}};\n")
    body.append("\\node[anchor=west,text width=6.5cm,font=\\bfseries] at (9.5,8.25) {Fully contained release groups excluded from primary matrices};\n")
    for i, accession in enumerate(("GSE108112", "GSE108109", "GSE104948", "GSE104954")):
        targets = sorted({r["other"] for r in contained if r["accession"] == accession})
        reciprocal = all(any(r["accession"] == t and r["other"] == accession for r in contained) for t in targets)
        symbol = "$\\equiv$" if reciprocal else "$\\subset$"
        target = " and ".join(targets)
        body.append(f"\\node[anchor=west,text width=6.5cm] at (9.5,{6.85 - i * 1.12}) {{{accession} {symbol} {target}}};\n")
    body.extend([
        "\\node[anchor=west,text width=15.5cm,text=midgrey] at (0.2,2.25) {Redundancy was resolved across every disease class, compartment and accession before matrix assembly.};\n",
        "\\node[anchor=west,text width=15.5cm,text=midgrey] at (0.2,1.25) {Alternative releases were retained only for declared annotation or normalisation sensitivity analyses.};\n",
        "\\end{tikzpicture}\n",
    ])
    compile_panel("Fig1_patient_deduplication", "".join(body))

    selected = [
        ("disc_glom_array", "Discovery array"),
        ("val_glom_array", "ERCB external"),
        ("cross_glom_rnaseq", "RNA-seq independent"),
        ("pla2r_lcm", "Laser capture"),
    ]
    lookup = {r["set"]: r for r in cohorts}
    body = ["\\begin{tikzpicture}[x=1cm,y=1cm]\n",
            "\\node[anchor=west,text width=15.5cm,font=\\bfseries\\fontsize{13}{15}\\selectfont] at (0,9.7) {Glomerular cohort allocation and independence};\n",
            "\\fill[teal] (10.7,8.65) rectangle (11.15,9.05); \\node[anchor=west,font=\\fontsize{8}{9}\\selectfont] at (11.3,8.85) {MN};\n",
            "\\fill[lightgrey,draw=midgrey] (12.65,8.65) rectangle (13.1,9.05); \\node[anchor=west,font=\\fontsize{8}{9}\\selectfont] at (13.25,8.85) {Other};\n"]
    for i, (key, label) in enumerate(selected):
        r = lookup[key]
        total = int(r["samples"])
        mn = int(r["MN"])
        other = total - mn
        y = 7.35 - i * 1.35
        scale = 0.032
        body.append(f"\\node[anchor=east] at (3.05,{y + 0.25}) {{{label}}};\n")
        body.append(f"\\fill[teal] (3.25,{y}) rectangle ({3.25 + mn * scale:.3f},{y + 0.52});\n")
        body.append(f"\\fill[lightgrey,draw=midgrey] ({3.25 + mn * scale:.3f},{y}) rectangle ({3.25 + total * scale:.3f},{y + 0.52});\n")
        body.append(f"\\node[anchor=west,text width=3cm,font=\\fontsize{8}{9}\\selectfont] at (12.9,{y + 0.25}) {{{total} total $|$ {mn} MN}};\n")
    body.extend([
        "\\draw[-{Stealth[length=2.5mm]},very thick,teal] (1.0,1.8) -- (1.65,1.8);\n",
        "\\node[anchor=west,text width=14cm] at (1.85,1.8) {ERCB and laser capture have 0 shared participants with discovery.};\n",
        "\\draw[-{Stealth[length=2.5mm]},very thick,orange] (1.0,0.8) -- (1.65,0.8);\n",
        "\\node[anchor=west,text width=14cm] at (1.85,0.8) {RNA-seq lacks participant identifiers and is treated as independent, not as repeated measures.};\n",
        "\\end{tikzpicture}\n",
    ])
    compile_panel("Fig1_cohort_allocation", "".join(body))

    r2_before = float(combat["batch_r2_pc1to5_before"])
    r2_after = float(combat["batch_r2_pc1to5_after"])
    batch_before = float(combat["batch_silhouette_before"])
    batch_after = float(combat["batch_silhouette_after"])
    disease_before = float(combat["disease_silhouette_before"])
    disease_after = float(combat["disease_silhouette_after"])
    reduction = 100 * (1 - r2_after / r2_before)
    cards = [
        (2.55, "Batch R$^2$ on PC1--5", f"{r2_before:.4f} $\\rightarrow$ {r2_after:.4f}", f"{reduction:.1f}\\% reduction", "teal"),
        (7.95, "Batch silhouette", f"{batch_before:.4f} $\\rightarrow$ {batch_after:.4f}", "series structure removed", "teal"),
        (13.35, "Disease silhouette", f"{disease_before:.4f} $\\rightarrow$ {disease_after:.4f}", f"$\\Delta$ = {disease_after - disease_before:+.4f}; preserved", "navy"),
    ]
    body = ["\\begin{tikzpicture}[x=1cm,y=1cm]\n",
            "\\node[anchor=west,text width=15.5cm,font=\\bfseries\\fontsize{13}{15}\\selectfont] at (0,9.7) {ComBat removes series structure without removing disease structure};\n"]
    for x, title, change, note, colour in cards:
        body.append(f"\\node[draw={colour},rounded corners=2pt,very thick,fill=white,minimum width=4.75cm,minimum height=4.2cm,text width=4.3cm,align=center] at ({x},4.35) {{\\textbf{{{title}}}\\\\[12pt]{{\\fontsize{{13}}{{15}}\\selectfont {change}}}\\\\[12pt]\\textcolor{{{colour}}}{{{note}}}}};\n")
    body.extend([
        f"\\node[anchor=west,text width=15.5cm,text=midgrey] at (0.2,1.35) {{ComBat: accession as batch and diagnosis as protected covariate; PC1 variance {combat['pc1_var_before']}\\% $\\rightarrow$ {combat['pc1_var_after']}\\%.}};\n",
        "\\node[anchor=west,text width=15.5cm,text=midgrey] at (0.2,0.35) {Quantitative batch correction is analytical validation, distinct from routine sample-level QC shown in the supplement.};\n",
        "\\end{tikzpicture}\n",
    ])
    compile_panel("Fig1_combat_quantitative", "".join(body))


if __name__ == "__main__":
    build_design_panels()
