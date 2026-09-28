#!/usr/bin/env python3
"""review_panel.py -- N random cells per sample as PNGs, barcode and stats stamped on each.

Renders hapCO's own per-cell plot (<barcode>_co.pdf, the view calls are judged
by) to PNG, stacks multi-page plots into one image, and stamps a banner:
    sample | NN/N | barcode | haploidness | molecules | COs called | label
The cell list is saved per seed, so the SAME panel can be re-rendered after a
parameter change and compared cell for cell.

Usage: review_panel.py SAMPLE N SEED [PER_CELL_DIR] [LABEL] [DPI]
  PER_CELL_DIR  default results/crossovers/SAMPLE/per_cell
  LABEL         default pipeline
  DPI           default 300 (the banner scales with it)
Output: qc/review/SAMPLE/LABEL/NN_BARCODE.png, a copy of the original vector
PDF next to each (NN_BARCODE.pdf, for lossless zooming), and panel.tsv
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

import matplotlib
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

S, N, SEED = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
SRC = sys.argv[4] if len(sys.argv) > 4 else "results/crossovers/%s/per_cell" % S
LABEL = sys.argv[5] if len(sys.argv) > 5 else "pipeline"
OUT = os.path.join("qc/review", S, LABEL)
DPI = int(sys.argv[6]) if len(sys.argv) > 6 else 300
SCALE = DPI / 110.0
os.makedirs(OUT, exist_ok=True)

if shutil.which("pdftoppm"):
    def render(pdf, prefix):
        subprocess.run(["pdftoppm", "-r", str(DPI), "-png", pdf, prefix], check=True)
elif shutil.which("gs"):
    def render(pdf, prefix):
        subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=png16m", "-r%d" % DPI,
                        "-sOutputFile=%s-%%d.png" % prefix, pdf], check=True)
else:
    sys.exit("neither pdftoppm nor gs is on PATH -- one of them is needed to turn the PDFs into PNGs")

cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
panel_file = os.path.join("qc/review", S, "panel_seed%d.txt" % SEED)
if os.path.exists(panel_file):
    pick = [l.strip() for l in open(panel_file) if l.strip()]
    print("reusing the saved panel %s (%d cells)" % (panel_file, len(pick)))
else:
    pick = random.Random(SEED).sample(cells, min(N, len(cells)))
    open(panel_file, "w").write("\n".join(pick) + "\n")
ht = pd.read_csv("qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S, sep="\t",
                 usecols=["barcode", "molecules", "haploidness"]).set_index("barcode")
font_path = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf")
font = ImageFont.truetype(font_path, int(22 * SCALE))

rows = []
tmp = tempfile.mkdtemp()
for i, bc in enumerate(pick, 1):
    pdf = os.path.join(SRC, "%s_co.pdf" % bc)
    calls = os.path.join(SRC, "%s_co_pred.txt" % bc)
    if not os.path.exists(pdf):
        print("  %02d %s: no PDF in %s -- skipped" % (i, bc, SRC))
        continue
    ncos = 0
    if os.path.exists(calls):
        for l in open(calls):
            f = l.split()
            if len(f) >= 3 and f[1].isdigit():
                ncos += 1
    prefix = os.path.join(tmp, bc)
    render(pdf, prefix)
    pages = sorted((p for p in os.listdir(tmp) if p.startswith(bc + "-") and p.endswith(".png")),
                   key=lambda p: int(p[len(bc) + 1:-4]))
    ims = [Image.open(os.path.join(tmp, p)).convert("RGB") for p in pages]
    h = ht.loc[bc] if bc in ht.index else None
    text = "%s   |   %02d/%d   |   %s   |   haploidness %s   |   %s molecules   |   %d COs called   |   %s" % (
        S, i, len(pick), bc, "%.2f" % h.haploidness if h is not None else "n/a",
        format(int(h.molecules), ",") if h is not None else "n/a", ncos, LABEL)
    text_w = int(ImageDraw.Draw(Image.new("RGB", (1, 1))).textlength(text, font=font)) + int(40 * SCALE)
    width = max(max(im.width for im in ims), text_w)
    banner = int(56 * SCALE)
    canvas = Image.new("RGB", (width, banner + sum(im.height for im in ims)), "white")
    d = ImageDraw.Draw(canvas)
    d.rectangle([0, 0, width, banner - int(6 * SCALE)], fill=(255, 244, 214))
    d.text((int(16 * SCALE), int(14 * SCALE)), text, fill=(20, 20, 20), font=font)
    y = banner
    for im in ims:
        canvas.paste(im, (0, y))
        y += im.height
    png = os.path.join(OUT, "%02d_%s.png" % (i, bc))
    canvas.save(png, optimize=True)
    shutil.copy(pdf, os.path.join(OUT, "%02d_%s.pdf" % (i, bc)))
    for p in pages:
        os.remove(os.path.join(tmp, p))
    rows.append((i, bc, h.haploidness if h is not None else None,
                 int(h.molecules) if h is not None else None, ncos, png))
shutil.rmtree(tmp, ignore_errors=True)
pd.DataFrame(rows, columns=["n", "barcode", "haploidness", "molecules", "cos_called", "png"]).to_csv(
    os.path.join(OUT, "panel.tsv"), sep="\t", index=False)
print("wrote %d PNGs at %d dpi (+ the vector PDFs) to %s (cell list: %s)" % (len(rows), DPI, OUT, panel_file))
