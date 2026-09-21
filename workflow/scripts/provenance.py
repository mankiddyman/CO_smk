#!/usr/bin/env python3
"""provenance.py -- the params and shell commands Snakemake RECORDED for a
sample's outputs, read from .snakemake/metadata. What the runs actually used --
authoritative where config is only a claim.

Usage: python3 provenance.py Spondias_hap1 6303_B_hap1
Outputs with missing metadata will not appear; use tool logs for those.
"""
import base64
import collections
import datetime
import json
import os
import re
import sys

CO = "/netscratch/dep_mercier/grp_marques/Aaryan/CO_smk"
MD = os.path.join(CO, ".snakemake", "metadata")


def decode(rel):
    s = rel.replace(os.sep, "")
    try:
        return base64.urlsafe_b64decode(s + "=" * (-len(s) % 4)).decode()
    except Exception:
        return None


def main(samples):
    per = collections.defaultdict(dict)
    for root, _, files in os.walk(MD):
        for f in files:
            full = os.path.join(root, f)
            path = decode(os.path.relpath(full, MD))
            if not path:
                continue
            hit = [s for s in samples
                   if re.search(r"(^|/)%s([/._]|$)" % re.escape(s), path)]
            if not hit:
                continue
            try:
                rec = json.load(open(full))
            except Exception:
                continue
            rule = rec.get("rule") or "?"
            for s in hit:
                prev = per[s].get(rule)
                if prev is None or (rec.get("endtime") or 0) > (prev[1].get("endtime") or 0):
                    per[s][rule] = (path, rec)
    for s in samples:
        print("=" * 96)
        print("%s -- %d rules with recorded metadata" % (s, len(per[s])))
        print("=" * 96)
        for rule in sorted(per[s]):
            path, rec = per[s][rule]
            t = rec.get("endtime")
            when = datetime.datetime.fromtimestamp(t).strftime("%Y-%m-%d") if t else "?"
            print("\nrule %s   (%s; e.g. %s)" % (rule, when, path))
            params = rec.get("params") or []
            if params:
                print("   params   : " + " | ".join(str(p) for p in params)[:700])
            sh = re.sub(r"\s+", " ", rec.get("shellcmd") or "").strip()
            if sh:
                print("   shellcmd : " + sh[:900] + (" ..." if len(sh) > 900 else ""))
        print()


if __name__ == "__main__":
    main(sys.argv[1:])
