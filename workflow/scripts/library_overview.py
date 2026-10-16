#!/usr/bin/env python3
"""library_overview.py -- one figure of the pollen scRNA libraries of D. binata and D. paradoxa: where the reads went
(STAR Log.final.out), the library metrics (STARsolo Summary.csv) and how many nuclei survive each step
(STARsolo's estimate, the pipeline's called and selected nuclei from the supervisor report v5).

Reads (run from the CO_smk root): results/starsolo/SAMPLE/Log.final.out, results/starsolo/SAMPLE/Solo.out/
{GeneFull,Gene}/Summary.csv, reports/supervisor_v5/numbers.txt (optional).
Writes reports/library_overview/library_overview.{png,tsv}. Usage: python workflow/scripts/library_overview.py
[--binata Dbinata_hap1] [--paradoxa Dparadoxa_split]. A few seconds.
"""
import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ap = argparse.ArgumentParser()
ap.add_argument("--binata", default="Dbinata_hap1")
ap.add_argument("--paradoxa", default="Dparadoxa_split")
ap.add_argument("--out", default="reports/library_overview")
A = ap.parse_args()
os.makedirs(A.out, exist_ok=True)
SP = [("binata", "D. binata", A.binata, "#0F6E56"), ("paradoxa", "D. paradoxa", A.paradoxa, "#534AB7")]
FATE = [("uniquely mapped", "#2a78d6"), ("multi-mapped", "#eb6834"), ("unmapped: too short", "#1baf7a"),
        ("unmapped: other", "#eda100")]
INK, MUTED, GRID = "#2C2C2A", "#6B6A64", "#E3E1DA"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def num(s):
    s = str(s).strip().rstrip("%")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def find_sample(prefix, preferred):
    if os.path.exists("results/starsolo/%s/Log.final.out" % preferred):
        return preferred
    hits = sorted(glob.glob("results/starsolo/%s*/Log.final.out" % prefix))
    return hits[0].split("/")[2] if hits else None


V5 = {}
if os.path.exists("reports/supervisor_v5/numbers.txt"):
    for l in open("reports/supervisor_v5/numbers.txt"):
        if "\t" in l:
            k, v = l.rstrip("\n").split("\t", 1)
            V5[k] = v

D = []
for key, name, pref, col in SP:
    s = find_sample("D" + key, pref)
    if s is None:
        print("%s: no results/starsolo/D%s*/Log.final.out, skipped" % (name, key))
        continue
    star = {k.strip(): v.strip() for k, v in (l.split("|", 1) for l in open("results/starsolo/%s/Log.final.out" % s)
                                              if "|" in l)}
    solo, feat = {}, ""
    for f in ("GeneFull", "Gene"):
        p = "results/starsolo/%s/Solo.out/%s/Summary.csv" % (s, f)
        if os.path.exists(p):
            solo = {k.strip(): num(v) for k, v in (l.strip().split(",", 1) for l in open(p) if "," in l)}
            feat = f
            break

    def sv(rx):
        k = next((k for k in solo if re.fullmatch(rx, k)), None)
        return solo.get(k, float("nan"))
    d = dict(key=key, name=name, sample=s, col=col, feat=feat,
             reads=num(star.get("Number of input reads", "nan")), length=num(star.get("Average input read length", "nan")),
             fate=[num(star.get("Uniquely mapped reads %", "nan")),
                   num(star.get("% of reads mapped to multiple loci", "nan")) +
                   num(star.get("% of reads mapped to too many loci", "nan")),
                   num(star.get("% of reads unmapped: too short", "nan")),
                   num(star.get("% of reads unmapped: other", "nan")) +
                   (num(star.get("% of reads unmapped: too many mismatches", "0")) or 0.0)],
             valid=sv(r"Reads With Valid Barcodes"), genome=sv(r"Reads Mapped to Genome: Unique"),
             genes=sv(r"Reads Mapped to Gene\w*: Unique Gene\w*"), incells=sv(r"Fraction of Unique Reads in Cells"),
             sat=sv(r"Sequencing Saturation"), est=sv(r"Estimated Number of Cells"),
             umi=sv(r"Median UMI per Cell"), ngenes=sv(r"Median Gene\w* per Cell"),
             called=num(V5.get("%s.cells_called" % key, "nan")), selected=num(V5.get("%s.cells_selected" % key, "nan")))
    D.append(d)
if not D:
    raise SystemExit("no STARsolo output found under results/starsolo/")

fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw={"width_ratios": [1.25, 1.15, 1]})
# a: where the reads went
ax = axs[0]
for i, d in enumerate(D):
    y, x = len(D) - 1 - i, 0.0
    for (lab, c), v in zip(FATE, d["fate"]):
        if v == v and v > 0:
            ax.barh(y, v, left=x, height=0.55, color=c, edgecolor="white", linewidth=2)
            if v >= 4:
                ax.text(x + v / 2, y, "%.0f%%" % v, ha="center", va="center", fontsize=8.5, color="white",
                        fontweight="bold")
            x += v
ax.set_yticks(range(len(D)))
ax.set_yticklabels(["%s\n%s M reads, %s bp" % (d["name"], format(round(d["reads"] / 1e6), ","), "%.0f" % d["length"])
                    for d in D[::-1]], fontsize=8.5)
ax.set_xlim(0, 100)
ax.set_xlabel("% of reads (STAR, Log.final.out)")
ax.set_title("a  Where the reads went", loc="left", fontsize=10)
ax.legend([Patch(color=c) for _, c in FATE], [l for l, _ in FATE], ncol=2, fontsize=7.5, frameon=False,
          loc="upper center", bbox_to_anchor=(0.5, -0.2))
# b: library metrics
ax = axs[1]
M = [("valid barcode", "valid"), ("mapped to the genome, unique", "genome"), ("mapped to genes, unique", "genes"),
     ("unique reads in cells", "incells"), ("sequencing saturation", "sat")]
h = 0.8 / len(D)
for j, d in enumerate(D):
    for i, (lab, k) in enumerate(M):
        v = d[k]
        if v == v:
            y = len(M) - 1 - i + (len(D) / 2 - 0.5 - j) * h
            ax.barh(y, 100 * v, height=h * 0.9, color=d["col"])
            ax.text(100 * v + 1.2, y, "%.0f%%" % (100 * v), va="center", fontsize=7.5, color=INK)
ax.set_yticks(range(len(M)))
ax.set_yticklabels([m[0] for m in M[::-1]], fontsize=8.5)
ax.set_xlim(0, 112)
ax.set_xlabel("% of reads (STARsolo, Summary.csv)")
ax.set_title("b  Library metrics", loc="left", fontsize=10)
ax.legend([Patch(color=d["col"]) for d in D], [d["name"] for d in D], fontsize=7.5, frameon=False, loc="lower right")
# c: nuclei
ax = axs[2]
C = [("STARsolo's estimate", "est"), ("nuclei called (pipeline)", "called"), ("haploid nuclei analysed", "selected")]
for j, d in enumerate(D):
    for i, (lab, k) in enumerate(C):
        v = d[k]
        if v == v and v > 0:
            y = len(C) - 1 - i + (len(D) / 2 - 0.5 - j) * h
            ax.barh(y, v, height=h * 0.9, color=d["col"])
            ax.text(v * 1.08, y, format(int(v), ","), va="center", fontsize=7.5, color=INK)
ax.set_xscale("log")
ax.set_yticks(range(len(C)))
ax.set_yticklabels([c[0] for c in C[::-1]], fontsize=8.5)
ax.set_xlim(right=max(d[k] for d in D for _, k in C if d[k] == d[k]) * 4)
ax.set_xlabel("nuclei (log scale)")
ax.set_title("c  Nuclei kept at each step", loc="left", fontsize=10)
med = "; ".join("%s: median %s UMIs and %s genes per STARsolo cell" % (
    d["name"], format(int(d["umi"]), ",") if d["umi"] == d["umi"] else "–",
    format(int(d["ngenes"]), ",") if d["ngenes"] == d["ngenes"] else "–") for d in D)
fig.text(0.01, -0.06, med + ". Samples: " + ", ".join("%s = %s (%s)" % (d["name"], d["sample"], d["feat"] or "no Summary.csv")
                                                       for d in D), fontsize=7.5, color=MUTED)
fig.suptitle("Pollen scRNA libraries: alignment and nuclei", x=0.01, ha="left", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
png = os.path.join(A.out, "library_overview.png")
fig.savefig(png, dpi=170, bbox_inches="tight")
rows = [("sample", "sample"), ("input reads", "reads"), ("read length, bp", "length"), ("valid barcode", "valid"),
        ("mapped to genome, unique (STARsolo)", "genome"), ("mapped to genes, unique", "genes"),
        ("unique reads in cells", "incells"), ("sequencing saturation", "sat"), ("STARsolo estimated cells", "est"),
        ("median UMIs per cell", "umi"), ("median genes per cell", "ngenes"), ("nuclei called", "called"),
        ("haploid nuclei analysed", "selected")]
with open(os.path.join(A.out, "library_overview.tsv"), "w") as fo:
    fo.write("metric\t" + "\t".join(d["name"] for d in D) + "\n")
    for i, (lab, _) in enumerate(FATE):
        fo.write("STAR: %s, %%\t%s\n" % (lab, "\t".join("%.2f" % d["fate"][i] for d in D)))
    for lab, k in rows:
        fo.write("%s\t%s\n" % (lab, "\t".join(str(d[k]) for d in D)))
print(open(os.path.join(A.out, "library_overview.tsv")).read())
print("wrote %s and library_overview.tsv" % png)
