#!/usr/bin/env python3
"""run_eta.py -- where is a running sample, and when will its landscape be ready?

Calling (one job per chromosome, named std_call_<chrom>): bcftools writes records in
position order, so the last position in the growing VCF over the chromosome's length
is the fraction done; with the job's elapsed time -> its finish time.
STARsolo (job std_scrna): reads processed so far (Log.progress.out) against the input
reads of the sample's previous STARsolo run -> mapping finish time, plus an allowance
for sorting, Solo counting and cell calling.
Final job: an allowance after both (markers, cellsnp, cells, crossovers, landscape).
Estimates assume the current speed holds.

Usage: run_eta.py SAMPLE PREVIOUS_SAMPLE [--final_h 2] [--after_star_h 0.75]
"""
import argparse
import datetime
import glob
import os
import re
import subprocess

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("previous")
ap.add_argument("--final_h", type=float, default=2.0)
ap.add_argument("--after_star_h", type=float, default=0.75)
A = ap.parse_args()
S, P = A.sample, A.previous
now = datetime.datetime.now()


def sh(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                          universal_newlines=True).stdout


def fmt(t):
    if t is None:
        return "?"
    return ("" if t.date() == now.date() else t.strftime("%a ")) + t.strftime("%H:%M")


jobs = {}
for l in sh("squeue -u $USER -h -o '%j|%T|%S'").splitlines():
    n, st, start = l.split("|")
    try:
        t0 = datetime.datetime.strptime(start, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        t0 = None
    jobs[n] = (st, t0)

lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open("results/reference/%s/genome.fa.fai" % S)}
print("RUN ETA  %s  (now %s)" % (S, now.strftime("%a %H:%M")))
print("  variant calling, one job per chromosome:")
call_eta = []
for c in lens:
    vcf = "results/markers/%s/by_chrom/%s.vcf.gz" % (S, c)
    st, t0 = jobs.get("std_call_" + c, (None, None))
    if st is None:
        done = os.path.exists(vcf + ".tbi")
        print("    %-11s %s" % (c, "done" if done else "NOT RUNNING and not finished"))
        call_eta.append(now if done else None)
        continue
    if st != "RUNNING" or t0 is None:
        print("    %-11s %s" % (c, st))
        call_eta.append(None)
        continue
    last = sh("zcat %s 2>/dev/null | tail -n 2 | head -n 1" % vcf).split("\t")
    pos = int(last[1]) if len(last) > 1 and last[1].isdigit() else 0
    frac = pos / float(lens[c])
    el = (now - t0).total_seconds()
    eta = t0 + datetime.timedelta(seconds=el / frac) if frac > 0.01 else None
    call_eta.append(eta)
    print("    %-11s %5.1f%% of %3.0f Mb after %4.1f h   -> done ~%s" % (c, 100 * frac, lens[c] / 1e6, el / 3600, fmt(eta)))

print("  STARsolo -> cells:")
st, t0 = jobs.get("std_scrna", (None, None))
star_eta = None
prog = glob.glob("results/starsolo/%s/**/*Log.progress.out" % S, recursive=True)
total = None
for f in glob.glob("results/starsolo/%s/**/*Log.final.out" % P, recursive=True):
    for l in open(f):
        if "Number of input reads" in l:
            total = int(l.split("|")[1])
if os.path.exists("results/cells/%s/barcodes_called.tsv" % S):
    print("    done (cells called)")
    star_eta = now
elif st is None:
    print("    NOT RUNNING and not finished")
elif not prog:
    print("    %s, no progress log yet (earlier steps or genome loading)" % st)
else:
    rows = [l.split() for l in open(prog[0]) if re.match(r"^\s*[A-Z][a-z]{2}\s+\d+\s+\d\d:\d\d:\d\d", l)]
    if rows and total:
        r = rows[-1]
        speed, done_reads = float(r[3]), int(r[4])
        rem_h = max(total - done_reads, 0) / (speed * 1e6) if speed > 0 else None
        if rem_h is not None:
            star_eta = now + datetime.timedelta(hours=rem_h + A.after_star_h)
        print("    %.0f of %.0f M reads (%.0f%%) at %.0f M/h   -> mapping done ~%s, cells ~%s"
              % (done_reads / 1e6, total / 1e6, 100.0 * done_reads / total, speed,
                 fmt(now + datetime.timedelta(hours=rem_h)) if rem_h is not None else "?", fmt(star_eta)))
    else:
        print("    running; " + ("no progress rows yet" if not rows else "previous run's read total not found"))

known = [e for e in call_eta if e is not None]
if len(known) == len(call_eta) and star_eta:
    start = max(max(known), star_eta)
    print("  landscape: final job starts ~%s, landscape ~%s (%.1f h allowance for the final job)"
          % (fmt(start), fmt(start + datetime.timedelta(hours=A.final_h)), A.final_h))
else:
    print("  landscape: not estimable yet (a step above has no estimate)")
