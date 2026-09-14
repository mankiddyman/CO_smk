#!/usr/bin/env python3
"""Shared helpers for the CO_smk QC battery.

Summary  -- collects key/value metrics, writes a machine-readable TSV
Flags    -- collects anomaly flags; each names the follow-up it should trigger
plotting -- consistent styling so figures are comparable across samples
"""
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_RAW   = "#d62728"   # raw / unfiltered
C_KEEP  = "#1f77b4"   # retained
C_REF   = "#2ca02c"   # reference / expectation
C_CUT   = "#000000"   # threshold lines


class Summary:
    """Machine-readable metric collector. One TSV row per metric."""

    def __init__(self, sample, stage):
        self.sample = sample
        self.stage = stage
        self.rows = []

    def add(self, key, value, note=""):
        self.rows.append((key, value, note))
        return value

    def write(self, path):
        with open(path, "w") as f:
            f.write("sample\tstage\tkey\tvalue\tnote\n")
            for k, v, n in self.rows:
                f.write(f"{self.sample}\t{self.stage}\t{k}\t{v}\t{n}\n")
        print(f"\nwrote {path}  ({len(self.rows)} metrics)", file=sys.stderr)

    def echo(self):
        if not self.rows:
            return
        w = max(len(k) for k, _, _ in self.rows)
        print(f"\n=== {self.sample} / {self.stage} ===", file=sys.stderr)
        for k, v, n in self.rows:
            line = f"  {k:<{w}}  {v}"
            if n:
                line += f"    # {n}"
            print(line, file=sys.stderr)


class Flags:
    """Anomaly flags. A flag does not answer a question -- it says one exists."""

    def __init__(self, sample, stage):
        self.sample = sample
        self.stage = stage
        self.flags = []

    def check(self, condition, name, message):
        """Raise `name` if `condition` is truthy. Returns the condition."""
        if condition:
            self.flags.append((name, message))
        return bool(condition)

    def write(self, path):
        with open(path, "w") as f:
            if not self.flags:
                f.write(f"# {self.sample} {self.stage}: no flags raised\n")
            for name, msg in self.flags:
                f.write(f"{name}\t{msg}\n")
        print(f"wrote {path}  ({len(self.flags)} flags)", file=sys.stderr)

    def echo(self):
        print(f"\n=== FLAGS: {self.sample} / {self.stage} ===", file=sys.stderr)
        if not self.flags:
            print("  none", file=sys.stderr)
        for name, msg in self.flags:
            print(f"  [{name}] {msg}", file=sys.stderr)


def run(cmd, **kw):
    """Run a shell command, raise on failure, return stdout."""
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        sys.exit(f"FAILED: {cmd}\n{p.stderr[:2000]}")
    return p.stdout


def stream(cmd):
    """Yield stdout lines from a shell command, raising if it exits non-zero."""
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line
    p.stdout.close()
    rc = p.wait()
    if rc != 0:
        sys.exit(f"FAILED (rc={rc}): {cmd}\n{p.stderr.read()[:2000]}")


def read_fai(path):
    """chrom -> length (bp), preserving file order."""
    lens = {}
    for line in open(path):
        f = line.split("\t")
        lens[f[0]] = int(f[1])
    return lens


def chrom_sort_key(c):
    """Sort chromosomes numerically when possible, else lexically."""
    try:
        return (0, int(c), "")
    except ValueError:
        return (1, 0, c)


def save(fig, path, dpi=130):
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"  wrote {path}", file=sys.stderr)
