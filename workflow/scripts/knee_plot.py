#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

samples = [
    ("Cuscuta hap1 (516M reads, 5' v1)", "#1f77b4",
     "results/starsolo/6303_B_hap1/umi_per_barcode.tsv"),
    ("Dbinata (52M reads, 3' v4)", "#d62728",
     "results/starsolo/Dbinata_10x/umi_per_barcode.tsv"),
    ("Dparadoxa (48M reads, 3' v4)", "#2ca02c",
     "results/starsolo/Dparadoxa_10x/umi_per_barcode.tsv"),
]

fig, ax = plt.subplots(figsize=(9, 6))

for label, color, path in samples:
    if not os.path.exists(path):
        print(f"MISSING: {path}")
        continue
    umi = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split()  # handle tab or space
            try:
                umi.append(int(parts[1]))
            except (IndexError, ValueError):
                print(f"BAD LINE in {path}: {repr(line)}")
                continue
    print(f"{label}: read {len(umi)} barcodes, max UMI {max(umi) if umi else 0}")
    umi = sorted(umi, reverse=True)
    ranks = np.arange(1, len(umi) + 1)
    ax.plot(ranks, umi, label=f"{label}: {len(umi):,} barcodes",
            color=color, linewidth=1.8)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1, 1e6)
ax.set_ylim(0.8, 1e4)
ax.set_xlabel("Barcode rank")
ax.set_ylabel("UMI count per barcode")
ax.set_title("Barcode-rank plot: pollen scRNA libraries")
ax.legend(loc="upper right", fontsize=10)
ax.grid(True, which="both", alpha=0.3)
for thresh, name in [(100, "100 UMI"), (1000, "1000 UMI")]:
    ax.axhline(thresh, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
    ax.text(1.5, thresh * 1.15, name, fontsize=8, color="gray")
plt.tight_layout()
plt.savefig("resources/pollen_qc/knee_plot.png", dpi=150)
print("Wrote knee_plot.png")
