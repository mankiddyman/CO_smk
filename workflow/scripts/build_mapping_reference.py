#!/usr/bin/env python3
"""build_mapping_reference.py -- the haploid mapping reference, cut from the dual assembly by name.

Every read in the crossover chain (HiFi for markers, single-cell for
genotypes) is mapped to ONE haploid reference. That is only valid if the
reference holds every piece of the genome exactly once. Hi-C phasing labels
hap1/hap2 consistently within a chromosome but arbitrarily across chromosomes,
so with a heterozygous reciprocal translocation 'all of hap1' can hold one arm
twice and the other never (D. paradoxa, 2026-09-28: chr1/chr2 right arms).

This builds the reference from an explicit list (config:
mapping_reference.<sample>.chromosomes, or .pieces), taken from the dual
assembly and its annotation as named in an existing reference's MANIFEST.txt.
Each entry is either
  CHROM                       a whole chromosome of the dual assembly, or
  NAME=CHROM:START-END        a piece of one (1-based, inclusive), named NAME;
and '@FILE' reads the entries from a split plan (breakpoints_from_paf.py: CSV,
or TSV, after '#' lines; columns piece, source, start, end). Pieces of one
chromosome may not overlap. Writes:
  genome.fa (60 bp lines, listed order), genome.fa.fai, contigs.txt,
  pieces.tsv (piece, source, start, end, length: where every sequence comes from),
  annotation.gff3 (features of the listed sequences, moved to piece coordinates;
  a gene cut by a piece boundary is left out with everything under it),
  MANIFEST.txt (source_fasta etc., so the rest of the pipeline can trace it like
  any published reference)

Usage: build_mapping_reference.py SOURCE_MANIFEST OUTDIR ENTRY[,ENTRY...]|@PLAN SAMPLE
"""
import collections
import csv
import datetime
import os
import re
import subprocess
import sys

man, out, spec, sample = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
meta = {}
for l in open(man):
    if ":" in l:
        k, v = l.split(":", 1)
        meta[k.strip()] = v.strip()
src_fa, src_gff = meta["source_fasta"], meta["source_gff"]

# ---- the entries: (name, source chromosome, start, end or None for the whole chromosome)
entries = []
if spec.startswith("@"):
    plan = spec[1:]
    rows = [l for l in open(plan) if l.strip() and not l.startswith("#")]
    for r in csv.DictReader(rows, delimiter="\t" if "\t" in rows[0] else ","):
        entries.append((r["piece"], r["source"], int(r["start"]), int(r["end"])))
else:
    plan = ""
    for e in spec.split(","):
        m = re.match(r"^([^=]+)=([^:]+):(\d+)-(\d+)$", e)
        entries.append((m.group(1), m.group(2), int(m.group(3)), int(m.group(4))) if m else (e, e, None, None))
names = [e[0] for e in entries]
if len(set(names)) != len(names):
    sys.exit("duplicate names: %s" % names)
need = set(e[1] for e in entries)
os.makedirs(out, exist_ok=True)

seqs, cur, buf = {}, None, []
with open(src_fa) as fh:                         # each chromosome joined as soon as it is read
    for l in fh:
        if l.startswith(">"):
            if cur in need:
                seqs[cur] = "".join(buf)
            cur, buf = l[1:].split()[0], []
            continue
        if cur in need:
            buf.append(l.strip())
    if cur in need:
        seqs[cur] = "".join(buf)
    buf = []
missing = sorted(need - set(seqs))
if missing:
    sys.exit("not in %s: %s" % (src_fa, missing))

# whole chromosomes get their full span; pieces are checked against the source length and each other
spans = collections.defaultdict(list)
full = []
for name, src, a, b in entries:
    n = len(seqs[src])
    if a is None:
        a, b = 1, n
    if not 1 <= a <= b <= n:
        sys.exit("%s: %s:%d-%d outside 1-%d" % (name, src, a, b, n))
    spans[src].append((a, b, name))
    full.append((name, src, a, b))
for src, lst in spans.items():
    lst.sort()
    for (a1, b1, n1), (a2, b2, n2) in zip(lst, lst[1:]):
        if a2 <= b1:
            sys.exit("%s and %s overlap on %s" % (n1, n2, src))
cover = {src: (sum(b - a + 1 for a, b, _ in lst), len(seqs[src])) for src, lst in spans.items()}

fai, total, off = [], 0, 0
with open(os.path.join(out, "genome.fa"), "w") as fo:
    for name, src, a, b in full:
        s = seqs[src][a - 1:b]
        hdr = ">%s\n" % name
        fo.write(hdr)
        off += len(hdr)
        fai.append("%s\t%d\t%d\t60\t61\n" % (name, len(s), off))
        for i in range(0, len(s), 60):
            line = s[i:i + 60] + "\n"
            fo.write(line)
            off += len(line)
        total += len(s)
open(os.path.join(out, "genome.fa.fai"), "w").write("".join(fai))
open(os.path.join(out, "contigs.txt"), "w").write("".join(n + "\n" for n, _, _, _ in full))
with open(os.path.join(out, "pieces.tsv"), "w") as fo:
    fo.write("piece\tsource\tstart\tend\tlength\n")
    for name, src, a, b in full:
        fo.write("%s\t%s\t%d\t%d\t%d\n" % (name, src, a, b, b - a + 1))

# ---- annotation: two passes. 1) which features cross a piece boundary (or lie outside every piece),
# and everything under them via Parent; 2) write the rest in piece coordinates.
by_src = collections.defaultdict(list)
for name, src, a, b in full:
    by_src[src].append((a, b, name))


def piece_of(src, s, e):
    for a, b, name in by_src.get(src, []):
        if a <= s and e <= b:
            return name, a
    return None, None


def attrs(col9):
    d = {}
    for kv in col9.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


cut_ids, parents, ids_on = set(), {}, set()
with open(src_gff) as fi:
    for l in fi:
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in by_src:
            continue
        a = attrs(f[8])
        fid = a.get("ID")
        if fid:
            ids_on.add(fid)
            parents[fid] = a.get("Parent", "").split(",") if a.get("Parent") else []
            if piece_of(f[0], int(f[3]), int(f[4]))[0] is None:
                cut_ids.add(fid)
changed = True
while changed:                                   # everything under a dropped feature goes too
    changed = False
    for fid, ps in parents.items():
        if fid not in cut_ids and any(p in cut_ids for p in ps):
            cut_ids.add(fid)
            changed = True

nfeat, ngene, ndrop, ngene_drop = 0, 0, 0, 0
with open(src_gff) as fi, open(os.path.join(out, "annotation.gff3"), "w") as fo:
    fo.write("##gff-version 3\n")
    for name, src, a, b in full:
        fo.write("##sequence-region %s 1 %d\n" % (name, b - a + 1))
    for l in fi:
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in by_src:
            continue
        a = attrs(f[8])
        fid = a.get("ID")
        pids = a.get("Parent", "").split(",") if a.get("Parent") else []
        pc, start = piece_of(f[0], int(f[3]), int(f[4]))
        if pc is None or (fid and fid in cut_ids) or any(p in cut_ids for p in pids):
            ndrop += 1
            ngene_drop += f[2] == "gene"
            continue
        f[0], f[3], f[4] = pc, str(int(f[3]) - start + 1), str(int(f[4]) - start + 1)
        fo.write("\t".join(f) + "\n")
        nfeat += 1
        ngene += f[2] == "gene"


def git(*a):
    r = subprocess.run(["git"] + list(a), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
    return r.stdout.strip()


dirty = git("status", "--porcelain", "--", "workflow", "config")
open(os.path.join(out, "MANIFEST.txt"), "w").write(
    "species:          %s\n"
    "haplotype:        composite\n"
    "built:            %s (rule mapping_reference, sample %s)\n"
    "source_fasta:     %s\n"
    "source_gff:       %s\n"
    "source_manifest:  %s\n"
    "chromosomes:      %s\n"
    "pieces:           %s\n"
    "split_plan:       %s\n"
    "repo_commit:      %s%s\n"
    "contigs:          %d\n"
    "total_bp:         %d\n"
    "gff_features:     %d\n"
    "gff_genes:        %d\n"
    "gff_dropped:      %d features (%d genes) cut by a piece boundary\n"
    % (meta.get("species", "?"), datetime.datetime.now().isoformat(timespec="seconds"), sample, src_fa, src_gff,
       os.path.abspath(man), ",".join(sorted(need)),
       ",".join("%s=%s:%d-%d" % (n, s, a, b) for n, s, a, b in full if (a, b) != (1, len(seqs[s]))) or "none",
       os.path.abspath(plan) if plan else "none", git("rev-parse", "--short", "HEAD"),
       "  (DIRTY WORKING TREE)" if dirty else "", len(full), total, nfeat, ngene, ndrop, ngene_drop))
print("%s: %d sequences, %.1f Mb, %d GFF features (%d genes; %d features, %d genes dropped at piece boundaries) -> %s"
      % (sample, len(full), total / 1e6, nfeat, ngene, ndrop, ngene_drop, out))
for (name, src, a, b), l in zip(full, fai):
    print("  %-18s %7.1f Mb  %s" % (name, int(l.split("\t")[1]) / 1e6,
                                    "" if (a, b) == (1, len(seqs[src])) else "%s:%d-%d" % (src, a, b)))
for src, (got, n) in sorted(cover.items()):
    if got != n:
        print("  NOTE: %s is %.1f Mb, its pieces cover %.1f Mb" % (src, n / 1e6, got / 1e6))
