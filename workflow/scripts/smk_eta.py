#!/usr/bin/env python3
"""smk_eta.py -- where is a one-job Snakemake run of CO_smk, and when will its landscape be ready?

Reads the run's Snakemake output (the sbatch log, logs/*_<JOBID>.out), which stamps every job's
start and finish, and follows the pipeline's critical path:
  HiFi branch   hifi_align: minimap2's "mapped N sequences" lines against the same reads' total
                in a previous sample's log, plus an allowance for sort + index
                -> calling, per chromosome: bcftools writes in position order, so the last
                position written over the chromosome length is the share done; or by 5 Mb
                region: regions done / planned at the measured pace and parallel width
                -> concat, merge, filter_markers
  scRNA branch  star_index -> starsolo_align (Log.progress.out against a previous sample's read
                total) -> filter_bam_uniq -> cell_calling
  then          cell_snp_counting -> ... -> recombination_landscape
Steps not yet running take the median duration of the same rule in a previous sample's
finished jobs (.snakemake/log and logs/run_*.out); per-chromosome calling takes that sample's
Mb per hour. Without a previous measurement a fixed guess is used and marked "guess".
Estimates assume the current speed holds. Also flags an estimate past the Slurm time limit.

Run from the CO_smk root. Usage: smk_eta.py JOBID SAMPLE PREVIOUS_SAMPLE [PREVIOUS_SAMPLE ...]
"""
import argparse
import datetime
import glob
import os
import re
import statistics
import subprocess
import zlib

import yaml

ap = argparse.ArgumentParser()
ap.add_argument("jobid"); ap.add_argument("sample"); ap.add_argument("previous", nargs="+")
A = ap.parse_args()
S, PREV = A.sample, A.previous
now = datetime.datetime.now()
H = datetime.timedelta(hours=1)
GUESS_H = {"hifi_align": 12, "star_index": 1.5, "starsolo_align": 6, "hifi_variants_concat_regions": 0.1,
           "hifi_variants_merge": 0.2, "filter_markers": 0.5, "filter_bam_uniq": 0.5, "cell_calling": 0.2,
           "cell_snp_counting": 4}
GUESS_OTHER_H, GUESS_MB_PER_H, SORT_H, STAR_TAIL_H = 0.3, 40.0, 0.75, 0.75
HIFI_TAIL = ["hifi_variants_merge", "filter_markers"]
CELL_TAIL = ["filter_bam_uniq", "cell_calling"]
FINAL = ["cell_snp_counting", "cellsnp_to_per_cell", "cell_data_molecules", "cell_haploidness", "select_cells",
         "co_calling", "co_aggregate", "marker_classes", "region_call_check", "recombination_landscape"]
TS = re.compile(r"^\[(\w{3} \w{3}\s+\d+ \d\d:\d\d:\d\d \d{4})\]$")


def sh(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                          universal_newlines=True).stdout.strip()


def fmt(t):
    if t is None:
        return "?"
    return (t.strftime("%a ") if t.date() != now.date() else "") + t.strftime("%H:%M")


def parse_log(path):
    """Snakemake log -> ({jobid: {rule, wc, start, end}}, planned jobs per rule, progress, failed rules, cores)."""
    jobs, plan, progress, errors, t, cur, stats, cores = {}, {}, None, [], None, None, False, None
    for line in open(path, errors="replace"):
        line = line.rstrip("\n")
        m = TS.match(line.strip())
        if m:
            t = datetime.datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y")
            continue
        if line.startswith("Job stats:"):
            stats = True
            continue
        if stats:
            f = line.split()
            if f and f[0] == "total":
                stats = False
            elif len(f) == 2 and f[1].isdigit() and f[0] != "job":
                plan.setdefault(f[0], int(f[1]))
            continue
        m = re.match(r"^(?:local)?(?:rule|checkpoint) (\S+):$", line)
        if m:
            cur = {"rule": m.group(1), "wc": {}, "start": t, "end": None}
            continue
        m = re.match(r"^\s+jobid: (\d+)$", line)
        if m and cur is not None:
            jobs[int(m.group(1))] = cur
        m = re.match(r"^\s+wildcards: (.*)$", line)
        if m and cur is not None:
            cur["wc"] = dict(kv.split("=", 1) for kv in m.group(1).split(", ") if "=" in kv)
        m = re.match(r"^Finished job(?:id:)? (\d+)", line)
        if m and int(m.group(1)) in jobs:
            jobs[int(m.group(1))]["end"] = t
        m = re.search(r"\d+ of \d+ steps \(\d+%\) done", line)
        if m:
            progress = m.group(0)
        m = re.match(r"^Error in rule (\S+):", line)
        if m:
            errors.append(m.group(1))
        m = re.match(r"^Provided cores: (\d+)", line)
        if m:
            cores = int(m.group(1))
    return jobs, plan, progress, errors, cores


def hours(j):
    return (j["end"] - j["start"]).total_seconds() / 3600


def segments(path):
    """minimap2 runs in a log, split at its 'loaded/built the index' line."""
    if not os.path.exists(path):
        return []
    lines = open(path, errors="replace").read().splitlines()
    cut = [i for i, l in enumerate(lines) if "the index for" in l]
    return [lines[a:b] for a, b in zip(cut, cut[1:] + [len(lines)])]


def mm_state(seg):
    """(sequences mapped, seconds spent mapping, finished) for one minimap2 run."""
    t0 = re.search(r"::([\d.]+)\*", seg[0])
    mapped, tl = 0, None
    for l in seg:
        m = re.search(r"worker_pipeline::([\d.]+)\*[\d.]+\] mapped (\d+) sequences", l)
        if m:
            mapped, tl = mapped + int(m.group(2)), float(m.group(1))
    return mapped, (tl - float(t0.group(1))) if (tl is not None and t0) else None, any("Real time:" in l for l in seg)


def vcf_last_pos(path):
    """Last position written to a growing bgzipped VCF (reads only the file's tail)."""
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 4 * 1024 * 1024))
            buf = f.read()
    except OSError:
        return None
    i = buf.find(b"\x1f\x8b\x08\x04")
    while i != -1 and buf[i + 12:i + 14] != b"BC":
        i = buf.find(b"\x1f\x8b\x08\x04", i + 1)
    text = b""
    while 0 <= i < len(buf):
        d = zlib.decompressobj(31)
        try:
            text += d.decompress(buf[i:])
        except zlib.error:
            break
        if not d.unused_data:
            break
        i = len(buf) - len(d.unused_data)
    for line in reversed(text.split(b"\n")[:-1]):
        f = line.split(b"\t")
        if len(f) > 1 and f[1].isdigit():
            return int(f[1])
    return None


# ---- this run
cands = sorted(glob.glob("logs/*_%s.out" % A.jobid))
run_log = next((p for p in cands if "jobid:" in open(p, errors="replace").read()), None)
if run_log is None:
    smk = sorted(glob.glob(".snakemake/log/*.snakemake.log"), key=os.path.getmtime, reverse=True)
    run_log = next((p for p in smk if "sample=%s" % S in open(p, errors="replace").read()), None)
if run_log is None:
    raise SystemExit("no Snakemake log for job %s / %s yet" % (A.jobid, S))
jobs, plan, progress, errors, cores = parse_log(run_log)
J = list(jobs.values())


def of(rule):
    return [j for j in J if j["rule"] == rule]


def finished(rule):
    return rule not in plan or sum(1 for j in of(rule) if j["end"]) >= plan[rule]


# ---- the previous sample's measured durations
prev = []
for p in glob.glob(".snakemake/log/*.snakemake.log") + glob.glob("logs/run_*.out"):
    if p == run_log:
        continue
    try:
        prev += [j for j in parse_log(p)[0].values()
                 if j["end"] and j["start"] and j["wc"].get("sample") in PREV and j["end"] > j["start"]]
    except (OSError, ValueError):
        pass


def fai(s):
    p = "results/reference/%s/genome.fa.fai" % s
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()} if os.path.exists(p) else {}


lens_cur = fai(S)
lens = {}
for s in PREV + [S]:
    lens.update(fai(s))
try:
    cbr = ((yaml.safe_load(open("config/config.yaml")) or {}).get("call_by_region") or {}).get(S) or {}
except (OSError, yaml.YAMLError):
    cbr = {}
region_chroms, region_mb = set(cbr.get("chromosomes") or []), float(cbr.get("mb") or 5)


def dur(rule):
    """(hours, basis) for one job of a rule not yet running."""
    v = [hours(j) for j in prev if j["rule"] == rule]
    return (statistics.median(v), "measured") if v else (GUESS_H.get(rule, GUESS_OTHER_H), "guess")


rates = [lens[j["wc"]["chrom"]] / 1e6 / hours(j) for j in prev
         if j["rule"] == "hifi_variants_per_chrom" and j["wc"].get("chrom") in lens]
rate, rate_basis = (statistics.median(rates), "measured") if rates else (GUESS_MB_PER_H, "guess")

# ---- Slurm
st = sh("squeue -j %s -h -o '%%T|%%S|%%l'" % A.jobid).split("|")
deadline = None
print("RUN ETA  %s  (now %s; Snakemake log %s)" % (S, now.strftime("%a %H:%M"), run_log))
if len(st) == 3:
    state, start, lim = st
    try:
        t0 = datetime.datetime.strptime(start, "%Y-%m-%dT%H:%M:%S")
        d, rest = lim.split("-") if "-" in lim else ("0", lim)
        p = [int(x) for x in rest.split(":")]
        p = [0] * (3 - len(p)) + p
        deadline = t0 + datetime.timedelta(days=int(d), hours=p[0], minutes=p[1], seconds=p[2])
    except ValueError:
        t0 = None
    print("  Slurm job %s: %s since %s, limit %s -> must finish by %s" % (A.jobid, state, fmt(t0), lim, fmt(deadline)))
else:
    print("  Slurm job %s is no longer in the queue (finished, failed or cancelled): sacct -j %s" % (A.jobid, A.jobid))
if errors:
    print("  !! ERROR in rule(s): %s -- the run will stop; see logs/<rule>/" % ", ".join(sorted(set(errors))))
run = [j for j in J if j["start"] and not j["end"]]
by_rule = {}
for j in run:
    by_rule.setdefault(j["rule"], []).append((now - j["start"]).total_seconds() / 3600)
print("  Snakemake: %s; running: %s" % (progress or "no step finished yet", ", ".join(
    "%s%s (%.1f h)" % (r, " x%d" % len(v) if len(v) > 1 else "", max(v)) for r, v in by_rule.items()) or "nothing"))

# ---- HiFi branch
print("  HiFi markers")
al = of("hifi_align")
if finished("hifi_align"):
    t_align = max([j["end"] for j in al if j["end"]] or [now])
    print("    hifi_align          done %s" % fmt(t_align))
elif al:
    total = None
    for s in PREV:
        done_runs = [g for g in segments("logs/hifi_align/%s.log" % s) if mm_state(g)[2]]
        if done_runs:
            total = mm_state(done_runs[-1])[0]
            break
    seg = segments("logs/hifi_align/%s.log" % S)
    mapped, map_s, fin = mm_state(seg[-1]) if seg else (0, None, False)
    if fin:
        t_align = now + SORT_H * H
        print("    hifi_align          mapping done; sort + index ~%.2f h -> ~%s" % (SORT_H, fmt(t_align)))
    elif total and mapped and map_s:
        rem_h = (total - mapped) / (mapped / map_s) / 3600
        t_align = now + (rem_h + SORT_H) * H
        print("    hifi_align          %.2f of %.2f M reads (%.0f%%) after %.1f h of mapping -> mapped ~%s, sorted ~%s"
              % (mapped / 1e6, total / 1e6, 100 * mapped / total, map_s / 3600, fmt(now + rem_h * H), fmt(t_align)))
    else:
        h, b = dur("hifi_align")
        t_align = al[0]["start"] + h * H
        print("    hifi_align          running, no progress lines yet (index building?); %.1f h (%s) -> ~%s"
              % (h, b, fmt(t_align)))
else:
    h, b = dur("hifi_align")
    t_align = now + h * H
    print("    hifi_align          not started; %.1f h (%s) -> ~%s" % (h, b, fmt(t_align)))

ends = []
pc = {j["wc"].get("chrom"): j for j in of("hifi_variants_per_chrom")}
n_pc = plan.get("hifi_variants_per_chrom", 0)
if n_pc:
    print("    calling, one bcftools job per chromosome (%.0f Mb/h, %s)" % (rate, rate_basis))
for c, j in sorted(pc.items()):
    L = lens.get(c)
    if j["end"]:
        ends.append(j["end"])
        print("      %-11s done %s" % (c, fmt(j["end"])))
        continue
    pos = vcf_last_pos("results/markers/%s/by_chrom/%s.vcf.gz" % (S, c))
    el = (now - j["start"]).total_seconds()
    if pos and L and pos / L > 0.01:
        eta = j["start"] + datetime.timedelta(seconds=el * L / pos)
        print("      %-11s %5.1f%% of %3.0f Mb after %4.1f h -> ~%s" % (c, 100 * pos / L, L / 1e6, el / 3600, fmt(eta)))
    else:
        eta = j["start"] + (L / 1e6 / rate if L else dur("hifi_variants_per_chrom")[0]) * H
        print("      %-11s started %.1f h ago, too early to measure -> ~%s" % (c, el / 3600, fmt(eta)))
    ends.append(eta)
for c in [c for c in lens_cur if c not in region_chroms and c not in pc
          and not os.path.exists("results/markers/%s/by_chrom/%s.vcf.gz.tbi" % (S, c))][:max(0, n_pc - len(pc))]:
    ends.append(t_align + lens_cur[c] / 1e6 / rate * H)
    print("      %-11s not started; %.1f h after hifi_align -> ~%s" % (c, lens_cur[c] / 1e6 / rate, fmt(ends[-1])))
R = plan.get("hifi_variants_region", 0)
if R:
    rj = of("hifi_variants_region")
    done = [j for j in rj if j["end"]]
    running = [j for j in rj if not j["end"]]
    v = [hours(j) for j in prev if j["rule"] == "hifi_variants_region"]
    each, basis = ((statistics.median([hours(j) for j in done]), "this run") if done else
                   (statistics.median(v), "measured") if v else (region_mb / rate, "from the per-chromosome rate"))
    pts = sorted([(j["start"], 1) for j in rj] + [(j["end"], -1) for j in done])
    width, cur_w = 0, 0
    for _, step in pts:
        cur_w += step
        width = max(width, cur_w)
    guess_w = not width
    width = width or max(1, (cores or 48) - 22)
    base = now if rj else max(now, t_align)
    t_reg = base + (R - len(done)) * each / width * H if len(done) < R else max(j["end"] for j in done)
    if not finished("hifi_variants_concat_regions"):
        t_reg += dur("hifi_variants_concat_regions")[0] * H
    ends.append(t_reg)
    print("    calling by 5 Mb region: %d of %d done, %d running; %.2f h each (%s), %d at a time%s -> all + concat ~%s"
          % (len(done), R, len(running), each, basis, width, " (guess: cores minus STAR and calling)" if guess_w else "",
             fmt(t_reg)))
t_calls = max(ends) if ends else t_align
t_markers = t_calls + sum(dur(r)[0] for r in HIFI_TAIL if not finished(r)) * H
print("    markers ready ~%s (merge + filter after the last calling job)" % fmt(t_markers))

# ---- scRNA branch
print("  cells (scRNA)")
ix = of("star_index")
if finished("star_index"):
    t_index = max([j["end"] for j in ix if j["end"]] or [now])
    print("    star_index          done %s" % fmt(t_index))
else:
    h, b = dur("star_index")
    t_index = (ix[0]["start"] if ix else now) + h * H
    print("    star_index          %s; %.1f h (%s) -> ~%s" % ("running" if ix else "not started", h, b, fmt(t_index)))
so = of("starsolo_align")
if finished("starsolo_align"):
    t_star = max([j["end"] for j in so if j["end"]] or [now])
    print("    starsolo_align      done %s" % fmt(t_star))
else:
    h, b = dur("starsolo_align")
    total = None
    for s in PREV:
        for f in glob.glob("results/starsolo/%s/Log.final.out" % s):
            for l in open(f):
                if "Number of input reads" in l:
                    total = int(l.split("|")[1])
        if total:
            break
    prog = "results/starsolo/%s/Log.progress.out" % S
    rows = [l.split() for l in open(prog)
            if re.match(r"^\s*[A-Z][a-z]{2}\s+\d+\s+\d\d:\d\d:\d\d", l)] if (so and os.path.exists(prog)) else []
    if rows and total and float(rows[-1][3]) > 0:
        speed, nreads = float(rows[-1][3]), int(rows[-1][4])
        rem = max(total - nreads, 0) / (speed * 1e6)
        t_star = now + (rem + STAR_TAIL_H) * H
        print("    starsolo_align      %.0f of %.0f M reads (%.0f%%) at %.0f M/h -> mapped ~%s, done ~%s"
              % (nreads / 1e6, total / 1e6, 100 * nreads / total, speed, fmt(now + rem * H), fmt(t_star)))
    else:
        t_star = (so[0]["start"] if so else t_index) + h * H
        print("    starsolo_align      %s; %.1f h (%s) -> ~%s" % ("running" if so else "not started", h, b, fmt(t_star)))
t_cells = t_star + sum(dur(r)[0] for r in CELL_TAIL if not finished(r)) * H
print("    cells called ~%s" % fmt(t_cells))

# ---- landscape
rest = [r for r in FINAL if r in plan and not finished(r)]
tail_h = sum(dur(r)[0] for r in rest)
guessed = [r for r in rest if dur(r)[1] == "guess"]
t_land = max(t_markers, t_cells) + tail_h * H
print("  landscape ~%s  (%s branch last, then %.1f h for %d steps from cellsnp to the landscape%s)"
      % (fmt(t_land), "HiFi" if t_markers >= t_cells else "scRNA", tail_h, len(rest),
         "; guessed: " + ", ".join(guessed) if guessed else ""))
if deadline:
    print("  %s" % ("fits the Slurm limit (%s)" % fmt(deadline) if t_land <= deadline else
                    "!! PAST the Slurm limit (%s): the job will be killed first -- resubmit the same sbatch "
                    "afterwards, it resumes where it stopped" % fmt(deadline)))
