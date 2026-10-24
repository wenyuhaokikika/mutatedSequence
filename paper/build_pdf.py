#!/usr/bin/env python3
"""Render the Markdown manuscript to a print-ready HTML and PDF."""

from __future__ import annotations

import subprocess
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent
MD = ROOT / "manuscript.md"
HTML = ROOT / "manuscript.html"
PDF = ROOT / "manuscript.pdf"

CSS = """
@page {
  size: A4;
  margin: 20mm 18mm 20mm 18mm;
}
html { font-size: 10.5pt; }
body {
  font-family: "Source Serif 4", "Charter", "Georgia", "Times New Roman", serif;
  color: #1f2a30;
  line-height: 1.48;
  max-width: 100%;
  margin: 0;
}
h1 {
  font-family: "Helvetica Neue", Arial, sans-serif;
  font-size: 20pt;
  line-height: 1.2;
  margin: 0 0 10pt 0;
  color: #22333B;
}
h2 {
  font-family: "Helvetica Neue", Arial, sans-serif;
  font-size: 13.5pt;
  margin: 20pt 0 6pt 0;
  padding-bottom: 2pt;
  border-bottom: 0.6pt solid #cfc8ba;
  color: #22333B;
  page-break-after: avoid;
}
h3 {
  font-family: "Helvetica Neue", Arial, sans-serif;
  font-size: 11.5pt;
  margin: 14pt 0 4pt 0;
  color: #2A6F6A;
  page-break-after: avoid;
}
p { margin: 0 0 7pt 0; text-align: justify; }
strong { color: #16242a; }
em { font-style: italic; }
a { color: #1f6f78; text-decoration: none; }
hr { border: 0; border-top: 0.5pt solid #d8d2c4; margin: 14pt 0; }
code {
  font-family: "SFMono-Regular", Consolas, monospace;
  font-size: 8.8pt;
  background: #f4f1e8;
  padding: 0 2pt;
  border-radius: 2pt;
}
pre {
  font-family: "SFMono-Regular", Consolas, monospace;
  font-size: 8.4pt;
  background: #f8f5ed;
  border: 0.5pt solid #e3ddd0;
  padding: 7pt;
  white-space: pre-wrap;
}
blockquote {
  margin: 8pt 0;
  padding: 6pt 10pt;
  border-left: 2.5pt solid #2A9D8F;
  background: #f3f8f6;
  color: #435158;
}
table {
  border-collapse: collapse;
  width: 100%;
  margin: 8pt 0 12pt 0;
  font-size: 8.3pt;
  page-break-inside: avoid;
}
th, td {
  border: 0.5pt solid #cfc8ba;
  padding: 4pt 5pt;
  vertical-align: top;
  text-align: left;
}
th { background: #f0ece1; font-weight: 600; }
ul, ol { margin: 3pt 0 8pt 0; padding-left: 18pt; }
li { margin-bottom: 2.5pt; }
img {
  max-width: 100%;
  display: block;
  margin: 10pt auto 3pt auto;
  page-break-inside: avoid;
}
figure { margin: 0; }
figcaption { font-size: 8.5pt; color: #4d5a60; }
.small { font-size: 8.5pt; color: #4d5a60; }
"""


def main() -> None:
    text = MD.read_text(encoding="utf-8")
    body = markdown.markdown(
        text,
        extensions=["tables", "sane_lists", "fenced_code", "nl2br"],
    )
    # Keep the draft notice visually distinct.
    body = body.replace(
        "<p><strong>Preprint draft for bioRxiv — not yet submitted</strong></p>",
        '<p class="small"><strong>Preprint draft for bioRxiv — not yet submitted</strong></p>',
    )
    HTML.write_text(
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<title>mutatedSequence — manuscript draft</title>\n"
        f"<style>{CSS}</style>\n</head>\n<body>\n{body}\n</body>\n</html>\n",
        encoding="utf-8",
    )

    edge = Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
    if not edge.exists():
        edge = Path(r"C:\Program Files\Microsoft\Edge\Application\msedge.exe")
    if not edge.exists():
        raise SystemExit("Microsoft Edge was not found; open manuscript.html and print to PDF")
    subprocess.run(
        [
            str(edge),
            "--headless",
            "--disable-gpu",
            "--no-pdf-header-footer",
            "--allow-file-access-from-files",
            f"--print-to-pdf={PDF}",
            HTML.resolve().as_uri(),
        ],
        check=True,
    )
    print(f"wrote {HTML.name} and {PDF.name}")


if __name__ == "__main__":
    main()
