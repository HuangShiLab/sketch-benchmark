#!/usr/bin/env python
"""Minimal Markdown -> .docx converter for the manuscript (no pandoc / python-docx needed).
Supports: # ## ### headings (first # = Title), paragraphs, pipe tables, **bold**, *italic*,
`code`, $math$ (rendered as italic text), lines starting with "- " or "N. " (kept as
indented paragraphs), [PLACEHOLDER ...] lines (italic).  Usage: md2docx.py in.md out.docx"""
import re, sys, zipfile
from xml.sax.saxutils import escape

def runs(text):
    """inline markdown -> list of <w:r> xml"""
    out = []
    tokens = re.split(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|`[^`]+`|\$[^$]+\$|\$\$[^$]+\$\$)", text)
    for t in tokens:
        if not t: continue
        props = ""; body = t
        if t.startswith("**") and t.endswith("**"): props = "<w:b/>"; body = t[2:-2]
        elif t.startswith("$$") and t.endswith("$$"): props = "<w:i/>"; body = t[2:-2]
        elif t.startswith("$") and t.endswith("$"): props = "<w:i/>"; body = t[1:-1]
        elif t.startswith("`") and t.endswith("`"): props = '<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas"/>'; body = t[1:-1]
        elif t.startswith("*") and t.endswith("*") and len(t) > 2: props = "<w:i/>"; body = t[1:-1]
        body = body.replace("\\|", "|")
        parts = body.split("\t")
        for i, p in enumerate(parts):
            if i: out.append("<w:r><w:tab/></w:r>")
            if p: out.append(f'<w:r>{"<w:rPr>" + props + "</w:rPr>" if props else ""}<w:t xml:space="preserve">{escape(p)}</w:t></w:r>')
    return "".join(out)

def para(text, style=None, indent=False, keep_next=False):
    ppr = ""
    if style or indent or keep_next:
        ppr = "<w:pPr>" + (f'<w:pStyle w:val="{style}"/>' if style else "") + ('<w:ind w:left="567" w:hanging="283"/>' if indent else "") + ("<w:keepNext/>" if keep_next else "") + "</w:pPr>"
    return f"<w:p>{ppr}{runs(text)}</w:p>"

def table(rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    cells = [r for r in cells if not all(re.fullmatch(r":?-{2,}:?", c) for c in r)]
    ncol = max(len(r) for r in cells)
    xml = ['<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="5000" w:type="pct"/><w:tblLook w:val="04A0"/></w:tblPr><w:tblGrid>' + "".join('<w:gridCol/>' for _ in range(ncol)) + "</w:tblGrid>"]
    for i, r in enumerate(cells):
        r = r + [""] * (ncol - len(r))
        xml.append("<w:tr>" + ("<w:trPr><w:tblHeader/></w:trPr>" if i == 0 else ""))
        for c in r:
            txt = f"**{c}**" if i == 0 and c and not c.startswith("**") else c
            xml.append(f'<w:tc><w:tcPr><w:tcW w:w="0" w:type="auto"/></w:tcPr>{para(txt, "TableText")}</w:tc>')
        xml.append("</w:tr>")
    xml.append("</w:tbl>"); xml.append(para(""))
    return "".join(xml)

def convert(md):
    body = []; lines = md.splitlines(); i = 0; first_h1 = True
    while i < len(lines):
        l = lines[i]
        if not l.strip(): i += 1; continue
        if l.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"): rows.append(lines[i]); i += 1
            body.append(table(rows)); continue
        if l.startswith("### "): body.append(para(l[4:], "Heading3", keep_next=True))
        elif l.startswith("## "): body.append(para(l[3:], "Heading2", keep_next=True))
        elif l.startswith("# "):
            body.append(para(l[2:], "Title" if first_h1 else "Heading1", keep_next=True)); first_h1 = False
        elif re.match(r"^(- |\d+\. )", l): body.append(para(l, "Normal", indent=True))
        elif l.startswith("[PLACEHOLDER"): body.append(para(f"*{l}*"))
        elif re.match(r"^(Table \d+\.|Figure \d+:|Box \d+\.)", l): body.append(para(l, "Caption", keep_next=l.startswith("Table") or l.startswith("Box")))
        else: body.append(para(l))
        i += 1
    return "".join(body)

STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:eastAsia="Calibri" w:cs="Calibri"/><w:sz w:val="22"/><w:lang w:val="en-GB"/></w:rPr></w:rPrDefault>
<w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="276" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:basedOn w:val="Normal"/><w:qFormat/><w:pPr><w:spacing w:after="240"/></w:pPr><w:rPr><w:b/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="360" w:after="120"/><w:outlineLvl w:val="0"/></w:pPr><w:rPr><w:b/><w:sz w:val="28"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:basedOn w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="240" w:after="80"/><w:outlineLvl w:val="1"/></w:pPr><w:rPr><w:b/><w:sz w:val="24"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading3"><w:name w:val="heading 3"/><w:basedOn w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="160" w:after="60"/><w:outlineLvl w:val="2"/></w:pPr><w:rPr><w:b/><w:i/><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="caption"/><w:basedOn w:val="Normal"/><w:qFormat/><w:rPr><w:sz w:val="20"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="TableText"><w:name w:val="Table Text"/><w:basedOn w:val="Normal"/><w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr><w:rPr><w:sz w:val="18"/></w:rPr></w:style>
<w:style w:type="table" w:default="1" w:styleId="TableNormal"><w:name w:val="Normal Table"/><w:tblPr><w:tblCellMar><w:top w:w="0" w:type="dxa"/><w:left w:w="108" w:type="dxa"/><w:bottom w:w="0" w:type="dxa"/><w:right w:w="108" w:type="dxa"/></w:tblCellMar></w:tblPr></w:style>
<w:style w:type="table" w:styleId="TableGrid"><w:name w:val="Table Grid"/><w:basedOn w:val="TableNormal"/><w:tblPr><w:tblBorders><w:top w:val="single" w:sz="4" w:color="808080"/><w:left w:val="single" w:sz="4" w:color="808080"/><w:bottom w:val="single" w:sz="4" w:color="808080"/><w:right w:val="single" w:sz="4" w:color="808080"/><w:insideH w:val="single" w:sz="4" w:color="808080"/><w:insideV w:val="single" w:sz="4" w:color="808080"/></w:tblBorders></w:tblPr></w:style>
</w:styles>"""
CT = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/></Types>"""
RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>"""
DOCRELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>"""

def main(src, dst):
    body = convert(open(src, encoding="utf8").read())
    doc = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'
           + body + '<w:sectPr><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0"/></w:sectPr></w:body></w:document>')
    with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", CT); z.writestr("_rels/.rels", RELS)
        z.writestr("word/document.xml", doc); z.writestr("word/styles.xml", STYLES); z.writestr("word/_rels/document.xml.rels", DOCRELS)
    print(f"{dst}: {len(body)} bytes of body xml")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
