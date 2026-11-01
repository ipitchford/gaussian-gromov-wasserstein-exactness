#!/usr/bin/env bash
# Rebuild the PDF from manuscript.md (Pandoc + XeLaTeX + xdvipdfmx).
set -euo pipefail
cd -- "$(dirname -- "$0")"
pandoc manuscript.md \
  --from=markdown+tex_math_single_backslash \
  --standalone \
  --include-in-header=paper_header.tex \
  --variable=geometry:margin=22mm \
  --variable=papersize:a4 --variable=fontsize:11pt \
  --variable=colorlinks:true \
  --resource-path=. --output=affine_rigidity_gaussian_gw.tex
xelatex -no-pdf -interaction=nonstopmode -halt-on-error affine_rigidity_gaussian_gw.tex > pdf_build.log
xelatex -no-pdf -interaction=nonstopmode -halt-on-error affine_rigidity_gaussian_gw.tex >> pdf_build.log
python3 - <<'PY'
from pathlib import Path
import os
import subprocess
import tempfile

root = Path.cwd()
with tempfile.TemporaryDirectory(prefix="gw-pdf-") as directory:
    temporary_pdf = Path(directory) / "paper.pdf"
    with (root / "pdf_build.log").open("ab") as log:
        subprocess.run(
            ["xdvipdfmx", "-q", "-o", str(temporary_pdf),
             str(root / "affine_rigidity_gaussian_gw.xdv")],
            check=True, stdout=log, stderr=log,
        )
    data = temporary_pdf.read_bytes()
    if not data.startswith(b"%PDF-") or b"%%EOF" not in data[-100:]:
        raise RuntimeError("PDF conversion did not produce a complete file")
    staged = root / "affine_rigidity_gaussian_gw.pdf.new"
    with staged.open("wb") as output:
        output.write(data)
        output.flush()
        os.fsync(output.fileno())
    staged.replace(root / "affine_rigidity_gaussian_gw.pdf")
PY
