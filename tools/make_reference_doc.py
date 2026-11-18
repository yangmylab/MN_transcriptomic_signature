#!/usr/bin/env python3
"""Build the project's pandoc reference document.

Typography (standing user preference):
  * CJK text        -> 仿宋 (FangSong, eastAsia)
  * Latin text      -> Times New Roman
  * code            -> red, monospaced (Consolas)

The template is pandoc's default reference.docx with styles.xml patched, so every
`pandoc in.md --reference-doc=reference/tnr_fangsong_redcode.docx -o out.docx`
gets the fonts without per-document flags.

Usage: python3 tools/make_reference_doc.py
"""

import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reference" / "tnr_fangsong_redcode.docx"

LATIN = "Times New Roman"
CJK = "FangSong"
MONO = "Consolas"
CODE_RED = "FF0000"

THEME_FONTS = (
    '<w:rFonts w:ascii="{latin}" w:hAnsi="{latin}" w:eastAsia="{cjk}" w:cs="{latin}"/>'
    .format(latin=LATIN, cjk=CJK)
)


def patch_styles(xml: str) -> str:
    # 1. every theme-driven rFonts becomes an explicit Times New Roman / 仿宋 pair
    def repl(match: re.Match) -> str:
        tag = match.group(0)
        if "Theme" not in tag and "theme" not in tag:
            return tag
        return THEME_FONTS

    xml = re.sub(r"<w:rFonts\b[^>]*/>", repl, xml, flags=re.S)

    # 2. document defaults, so unstyled runs are covered as well
    if "<w:rPrDefault>" in xml and "<w:rFonts" not in xml.split("<w:rPrDefault>")[1].split("</w:rPrDefault>")[0]:
        xml = xml.replace("<w:rPrDefault>", "<w:rPrDefault>" + THEME_FONTS, 1)

    # 3. Normal style: explicit fonts for anything inheriting from it
    normal = re.search(r'<w:style [^>]*w:styleId="Normal".*?</w:style>', xml, re.S)
    if normal is None:
        sys.exit("Normal style not found in reference.docx")
    block = normal.group(0)
    if "<w:rPr>" in block:
        patched = block.replace("<w:rPr>", "<w:rPr>" + THEME_FONTS, 1)
    else:
        patched = block.replace("</w:style>", "<w:rPr>" + THEME_FONTS + "</w:rPr></w:style>")
    xml = xml.replace(block, patched, 1)

    # 4. inline code (VerbatimChar) and code blocks (SourceCode) in red monospace
    verbatim = re.search(r'<w:style [^>]*w:styleId="VerbatimChar".*?</w:style>', xml, re.S)
    if verbatim is None:
        sys.exit("VerbatimChar style not found in reference.docx")
    vb = verbatim.group(0)
    if "<w:color " not in vb:
        vb = vb.replace("<w:rPr>", '<w:rPr><w:color w:val="%s"/>' % CODE_RED, 1)
    vb = re.sub(r'<w:rFonts\b[^>]*/>',
                '<w:rFonts w:ascii="%s" w:hAnsi="%s" w:eastAsia="%s"/>' % (MONO, MONO, CJK),
                vb)
    xml = xml.replace(verbatim.group(0), vb, 1)

    source_code = (
        '<w:style w:type="paragraph" w:customStyle="1" w:styleId="SourceCode">\n'
        '  <w:name w:val="Source Code" />\n'
        '  <w:basedOn w:val="Normal" />\n'
        '  <w:pPr><w:shd w:val="clear" w:color="auto" w:fill="F2F2F2" /></w:pPr>\n'
        '  <w:rPr>\n'
        '    <w:rFonts w:ascii="%s" w:hAnsi="%s" w:eastAsia="%s" />\n'
        '    <w:color w:val="%s" />\n'
        '    <w:sz w:val="20" />\n'
        '  </w:rPr>\n'
        '</w:style>\n' % (MONO, MONO, CJK, CODE_RED)
    )
    if 'w:styleId="SourceCode"' not in xml:
        xml = xml.replace("</w:styles>", source_code + "</w:styles>", 1)

    return xml


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        default = tmp / "default.docx"

        with default.open("wb") as fh:
            subprocess.run(["pandoc", "--print-default-data-file", "reference.docx"],
                           check=True, stdout=fh)

        with zipfile.ZipFile(default) as z:
            names = z.namelist()
            z.extractall(tmp / "x")

        styles = tmp / "x" / "word" / "styles.xml"
        styles.write_text(patch_styles(styles.read_text(encoding="utf-8")), encoding="utf-8")

        if OUT.exists():
            OUT.unlink()

        with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
            for name in names:
                z.write(tmp / "x" / name, name)

    print("wrote", OUT, OUT.stat().st_size, "bytes")


if __name__ == "__main__":
    main()
