#!/usr/bin/env python3
"""Build the HRCP bibliography as a minimal, valid .docx (OOXML) with no deps."""
import zipfile
from xml.sax.saxutils import escape

TITLE = "Bibliography"
SUB = [
    "Muon-Induced Background in the Unanalyzed Years of the "
    "DM-Ice17 Dark Matter Detector",
    "Ben Charette · Harvard College Research Program, Fall 2026",
]

# Each entry: list of (text, italic) segments.
REFS = [
    [("1.\tBernabei, R., et al. (DAMA/LIBRA Collaboration). (2018). First model "
      "independent results from DAMA/LIBRA-phase2. ", False),
     ("Universe", True),
     (", 4(11), 116. doi: 10.3390/universe4110116", False)],

    [("2.\tCherwinka, J., et al. (DM-Ice Collaboration). (2014). First data from "
      "DM-Ice17. ", False),
     ("Physical Review D", True),
     (", 90, 092005. doi: 10.1103/PhysRevD.90.092005", False)],

    [("3.\tAartsen, M. G., et al. (IceCube Collaboration). (2017). The IceCube "
      "Neutrino Observatory: Instrumentation and online systems. ", False),
     ("Journal of Instrumentation", True),
     (", 12, P03012. doi: 10.1088/1748-0221/12/03/P03012", False)],

    [("4.\tCherwinka, J., et al. (DM-Ice Collaboration). (2016). Measurement of muon "
      "annual modulation and muon-induced phosphorescence in NaI(Tl) crystals with "
      "DM-Ice17. ", False),
     ("Physical Review D", True),
     (". arXiv:1509.02486", False)],

    [("5.\tAdhikari, G., et al. (COSINE-100 Collaboration). (2018). An experiment to "
      "search for dark-matter interactions using sodium iodide detectors. ", False),
     ("Nature", True),
     (", 564, 83–86. doi: 10.1038/s41586-018-0739-1", False)],

    [("6.\tHubbard, A. J. F. (2015). ", False),
     ("Muon-Induced Backgrounds in the DM-Ice17 NaI(Tl) Dark Matter Detector", True),
     (" (PhD thesis). University of Wisconsin–Madison.", False)],

    [("7.\tTilav, S., et al. (IceCube Collaboration). (2010). Atmospheric variations "
      "as observed by IceCube. arXiv:1001.0776", False)],
]


def run(text, italic=False, bold=False):
    rpr = ""
    if bold:
        rpr += "<w:b/>"
    if italic:
        rpr += "<w:i/>"
    rpr = f"<w:rPr>{rpr}</w:rPr>" if rpr else ""
    return (f'<w:r>{rpr}<w:t xml:space="preserve">'
            f'{escape(text)}</w:t></w:r>')


def para(segments, *, align=None, hanging=False, after=240, line=480, bold=False):
    ppr = "<w:pPr>"
    if align:
        ppr += f'<w:jc w:val="{align}"/>'
    if hanging:
        ppr += '<w:ind w:left="720" w:hanging="720"/>'
    ppr += f'<w:spacing w:after="{after}" w:line="{line}" w:lineRule="auto"/>'
    ppr += "</w:pPr>"
    runs = "".join(run(t, italic=i, bold=bold) for t, i in segments)
    return f"<w:p>{ppr}{runs}</w:p>"


body = []
body.append(para([(TITLE, False)], align="center", after=120, line=240, bold=True))
for line in SUB:
    body.append(para([(line, False)], align="center", after=120, line=240))
body.append(para([("", False)], after=120, line=240))
for r in REFS:
    body.append(para(r, hanging=True))

DOCUMENT = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    f'<w:body>{"".join(body)}'
    '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>'
    '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440"'
    ' w:header="720" w:footer="720" w:gutter="0"/></w:sectPr>'
    '</w:body></w:document>'
)

STYLES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    '<w:docDefaults><w:rPrDefault><w:rPr>'
    '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"'
    ' w:eastAsia="Times New Roman" w:cs="Times New Roman"/>'
    '<w:sz w:val="24"/><w:szCs w:val="24"/>'
    '</w:rPr></w:rPrDefault>'
    '<w:pPrDefault><w:pPr><w:spacing w:after="240" w:line="480"'
    ' w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>'
    '<w:style w:type="paragraph" w:default="1" w:styleId="Normal">'
    '<w:name w:val="Normal"/><w:qFormat/></w:style>'
    '</w:styles>'
)

CONTENT_TYPES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
    '<Default Extension="rels"'
    ' ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
    '<Default Extension="xml" ContentType="application/xml"/>'
    '<Override PartName="/word/document.xml" ContentType="application/vnd.'
    'openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
    '<Override PartName="/word/styles.xml" ContentType="application/vnd.'
    'openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
    '</Types>'
)

RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/'
    '2006/relationships/officeDocument" Target="word/document.xml"/>'
    '</Relationships>'
)

DOC_RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/'
    '2006/relationships/styles" Target="styles.xml"/>'
    '</Relationships>'
)

import sys
out = sys.argv[1]
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("[Content_Types].xml", CONTENT_TYPES)
    z.writestr("_rels/.rels", RELS)
    z.writestr("word/document.xml", DOCUMENT)
    z.writestr("word/styles.xml", STYLES)
    z.writestr("word/_rels/document.xml.rels", DOC_RELS)
print("wrote", out)
