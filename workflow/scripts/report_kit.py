#!/usr/bin/env python3
"""report_kit.py -- building blocks for code-generated evidence reports.

Copy this file next to a project's report script (e.g. workflow/scripts/report_kit.py) and import it.
It has no dependencies beyond numpy and matplotlib, so it runs wherever the pipeline runs.

What it gives you
  Report        one self-contained HTML page plus numbers.txt, built section by section
    .have(p)          record a file as read (True) or missing (False)
    .put(k, v)        log a number under a key; returns v, so it can wrap the value where it is used
    .guarded(n, fn)   run fn(); on any error record it, show it in red at the top, and carry on
    .h2/.h3/.p        headings and paragraphs (text is HTML; build numbers into it with f-strings)
    .table(head, rows)
    .figure(fig, name, width)   save a matplotlib figure under fig/ and embed it
    .img(path, width)           embed an existing PNG (a missing one shows as "not built")
    .details(title, path=, text=)  a collapsed block with a file's text (evidence, logs)
    .write()          writes report.html and numbers.txt; prints the numbers, the problems, the missing files
  Formatters    fi (integer with thousands separator), ff (fixed decimals), ok (not None/NaN)
  Provenance    git_commit(), git_tag(pattern)
  Code facts    code_const(path, regex, n): read a constant from pipeline code so the report's definitions
                match the code that produced the data; tool(name): find a binary in .snakemake/conda or PATH
  Statistics    poisson_tail, binom_p, bootstrap_ci, chance_runs (for "is this gap or hotspot more than chance?")
"""
import base64
import datetime
import glob
import html
import math
import os
import re
import shutil
import subprocess
import traceback

import numpy as np


# ------------------------------------------------------------------------------------ formatting
def ok(x):
    return x is not None and not (isinstance(x, float) and (math.isnan(x) or math.isinf(x)))


def fi(x):
    """Integer with thousands separators; '–' when missing."""
    try:
        return "{:,}".format(int(round(float(x)))) if ok(float(x)) else "–"
    except (TypeError, ValueError):
        return "–"


def ff(x, d=1):
    """Fixed decimals; '–' when missing."""
    try:
        return "%.*f" % (d, float(x)) if ok(float(x)) else "–"
    except (TypeError, ValueError):
        return "–"


def natural(s):
    """Sort key: chr2 before chr10."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", str(s))]


# ------------------------------------------------------------------------------------ provenance, code facts
def git_commit():
    return subprocess.run("git rev-parse --short HEAD", shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip() or "?"


def git_tag(pattern="*"):
    return subprocess.run("git tag -l '%s' | tail -n 1" % pattern, shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip()


def code_const(path, pattern, n=1):
    """Read numbers from pipeline code, e.g. code_const('workflow/scripts/x.py', r'MIN_CELLS\\s*=\\s*(\\d+)').
    Returns a list of n floats (NaN when not found), so definitions in the text cannot drift from the code."""
    try:
        m = re.search(pattern, open(path).read())
    except OSError:
        m = None
    return [float(x) for x in m.groups()][:n] if m else [float("nan")] * n


def tool(name):
    """A binary from the pipeline's conda environments first, then PATH."""
    c = sorted(glob.glob(".snakemake/conda/*/bin/%s" % name))
    return os.path.abspath(c[0]) if c else shutil.which(name)


# ------------------------------------------------------------------------------------ statistics
def poisson_tail(k, mu):
    """P(X >= k) for X ~ Poisson(mu), computed in log space."""
    if k <= 0:
        return 1.0
    if mu <= 0:
        return 0.0
    s = sum(math.exp(-mu + i * math.log(mu) - math.lgamma(i + 1)) for i in range(int(k)))
    return max(0.0, 1.0 - s)


def binom_p(k, n, p=0.5):
    """Two-sided exact binomial test of k successes in n against p."""
    if n == 0:
        return 1.0
    pk = [math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(n + 1)]
    return min(1.0, sum(q for q in pk if q <= pk[k] * (1 + 1e-9)))


def bootstrap_ci(values, stat=np.mean, n=1000, seed=1, level=95):
    """Percentile bootstrap interval of stat over the rows of values (resampling units, e.g. cells or genes)."""
    rng = np.random.default_rng(seed)
    v = np.asarray(values)
    reps = [stat(v[rng.integers(0, len(v), len(v))]) for _ in range(n)]
    a = (100 - level) / 2
    return float(stat(v)), float(np.percentile(reps, a)), float(np.percentile(reps, 100 - a))


def chance_runs(n_events, length, window, min_len, reps=200, seed=1):
    """How many event-free stretches of >= min_len appear by chance when n_events fall uniformly on length,
    counted in windows of `window`. Returns (mean, 2.5th, 97.5th percentile). Use it before calling a gap real."""
    rng = np.random.default_rng(seed)
    nb = int(math.ceil(length / window))
    out = []
    for _ in range(reps):
        cnt = np.bincount(np.minimum((rng.uniform(0, length, n_events) // window).astype(int), nb - 1), minlength=nb)
        runs, i = 0, 0
        while i < nb:
            if cnt[i] == 0:
                j = i
                while j + 1 < nb and cnt[j + 1] == 0:
                    j += 1
                if (min((j + 1) * window, length) - i * window) >= min_len:
                    runs += 1
                i = j + 1
            else:
                i += 1
        out.append(runs)
    return float(np.mean(out)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


# ------------------------------------------------------------------------------------ the report
CSS = """body{font-family:Helvetica,Arial,sans-serif;max-width:1080px;margin:24px auto;padding:0 16px;color:#2C2C2A;
font-size:13px;line-height:1.45}h1{font-size:21px}h2{font-size:16px;border-bottom:1px solid #D3D1C7;margin-top:30px}
h3{font-size:13px;margin-top:18px}table{border-collapse:collapse;margin:8px 0 14px;font-size:12px}td,th{border:1px
solid #D3D1C7;padding:3px 8px;text-align:left;vertical-align:top}th{background:#F1EFE8}.meta{color:#5F5E5A;font-size:11px}
.missing{color:#A32D2D}.look{color:#3C3489}.path{color:#888780;font-size:11px}code{font-size:11px}pre{font-size:10.5px;
background:#F7F6F2;padding:8px;overflow-x:auto}img{display:block;margin:8px 0}details{margin:4px 0}
@media print{img{page-break-inside:avoid}table{page-break-inside:avoid}}"""


class Report:
    def __init__(self, title, out_dir, generator):
        self.title, self.out, self.generator = title, out_dir, generator
        self.fig_dir = os.path.join(out_dir, "fig")
        os.makedirs(self.fig_dir, exist_ok=True)
        self.parts, self.used, self.missing, self.num, self.problems = [], [], [], {}, []

    # ---- bookkeeping
    def have(self, p):
        if p and os.path.exists(p):
            if p not in self.used:
                self.used.append(p)
            return True
        if p not in self.missing:
            self.missing.append(p)
        return False

    def put(self, key, val):
        self.num[key] = val
        return val

    def guarded(self, name, fn):
        try:
            return fn()
        except Exception as e:
            self.problems.append((name, "%s: %s" % (type(e).__name__, e)))
            traceback.print_exc()
            return None

    # ---- content
    def h2(self, t):
        self.parts.append("<h2>%s</h2>" % t)

    def h3(self, t):
        self.parts.append("<h3>%s</h3>" % t)

    def p(self, t):
        self.parts.append("<p>%s</p>" % t)

    def look(self, t):
        """'What to look for' line under a figure: what the reader should see if the claim holds."""
        self.parts.append("<p class='look'>What to look for: %s</p>" % t)

    def table(self, head, rows):
        h = "".join("<th>%s</th>" % html.escape(str(x)) for x in head)
        b = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % html.escape(str(x)) for x in r) for r in rows)
        self.parts.append("<table><tr>%s</tr>%s</table>" % (h, b))

    def img(self, path, width=100):
        if not path or not os.path.exists(str(path)):
            self.parts.append("<p class='missing'>figure not built: %s</p>" % html.escape(str(path)))
            return
        if str(path) not in self.used and not str(path).startswith(self.fig_dir):
            self.used.append(str(path))
        b = base64.b64encode(open(path, "rb").read()).decode()
        self.parts.append("<img src='data:image/png;base64,%s' style='width:%d%%'>" % (b, width))

    def figure(self, fig, name, width=100, dpi=150):
        import matplotlib.pyplot as plt
        p = os.path.join(self.fig_dir, name if name.endswith(".png") else name + ".png")
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        self.img(p, width)
        return p

    def details(self, title, path=None, text=None):
        if text is None:
            if not self.have(path):
                self.parts.append("<p class='missing'>missing: %s</p>" % html.escape(str(path)))
                return
            text = open(path, errors="replace").read()
        self.parts.append("<details><summary>%s <span class='path'>%s</span></summary><pre>%s</pre></details>"
                          % (html.escape(title), html.escape(path or ""), html.escape(text)))

    # ---- output
    def write(self):
        head = ["<h1>%s</h1>" % self.title,
                "<p class='meta'>Generated %s from commit %s by %s. Every number, table and figure is computed "
                "from the pipeline's outputs; the files read are listed in the appendix.</p>"
                % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), git_commit(), self.generator)]
        if self.problems:
            head.append("<p class='missing'>Parts that could not be built: %s</p>" % html.escape(
                "; ".join("%s (%s)" % p for p in self.problems)))
        tail = ["<h3>Files read</h3><pre>%s</pre>" % html.escape("\n".join(self.used))]
        if self.missing:
            tail.append("<h3>Files not found</h3><pre>%s</pre>" % html.escape("\n".join(map(str, self.missing))))
        body = "\n".join(head + self.parts + tail)
        with open(os.path.join(self.out, "report.html"), "w") as f:
            f.write("<!doctype html><html><head><meta charset='utf-8'><title>%s</title><style>%s</style></head>"
                    "<body>%s</body></html>" % (html.escape(re.sub("<[^>]+>", "", self.title)), CSS, body))
        lines = ["%s\t%s" % (k, ("%.4g" % v) if isinstance(v, (float, np.floating)) else v) for k, v in sorted(self.num.items())]
        with open(os.path.join(self.out, "numbers.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
        print("\n".join(lines))
        print("\nparts not built: %d%s" % (len(self.problems), "".join("\n  %s: %s" % p for p in self.problems)))
        print("missing files: %d%s" % (len(self.missing), "".join("\n  %s" % m for m in self.missing)))
        print("wrote %s" % os.path.join(self.out, "report.html"))
