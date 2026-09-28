#!/usr/bin/env python3
"""build_mapping_reference.py -- the haploid mapping reference, cut from the dual assembly by name.

Every read in the crossover chain (HiFi for markers, single-cell for
genotypes) is mapped to ONE haploid reference. That is only valid if the
reference holds every piece of the genome exactly once. Hi-C phasing labels
hap1/hap2 consistently within a chromosome but arbitrarily across chromosomes,
so with a heterozygous reciprocal translocation 'all of hap1' can hold one arm
twice and the other never (D. paradoxa, 2026-09-28: chr1/chr2 right arms).

This builds the reference from an explicit list of chromosomes (config:
mapping_reference.<sample>.chromosomes), taken from the dual assembly and its
annotation as named in an existing reference's MANIFEST.txt, and writes:
  genome.fa (60 bp lines, listed order), genome.fa.fai, contigs.txt,
  annotation.gff3 (the listed sequences' features only), MANIFEST.txt
  (source_fasta etc., so the rest of the pipeline can trace it like any
  published reference)

Usage: build_mapping_reference.py SOURCE_MANIFEST OUTDIR CHROM[,CHROM...] SAMPLE
"""
import datetime
import os
import subprocess
import sys

man, out, chroms, sample = sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4]
meta = {}
for l in open(man):
    if ":" in l:
        k, v = l.split(":", 1)
        meta[k.strip()] = v.strip()
src_fa, src_gff = meta["source_fasta"], meta["source_gff"]
if len(set(chroms)) != len(chroms):
    sys.exit("chromosome list has duplicates: %s" % chroms)
os.makedirs(out, exist_ok=True)

seqs, cur = {}, None
with open(src_fa) as fh:
    for l in fh:
        if l.startswith(">"):
            cur = l[1:].split()[0]
            if cur in chroms:
                seqs[cur] = []
            continue
        if cur in seqs:
            seqs[cur].append(l.strip())
missing = [c for c in chroms if c not in seqs]
if missing:
    sys.exit("not in %s: %s" % (src_fa, missing))

fai, total, off = [], 0, 0
with open(os.path.join(out, "genome.fa"), "w") as fo:
    for c in chroms:
        s = "".join(seqs.pop(c))
        hdr = ">%s\n" % c
        fo.write(hdr)
        off += len(hdr)
        fai.append("%s\t%d\t%d\t60\t61\n" % (c, len(s), off))
        for i in range(0, len(s), 60):
            line = s[i:i + 60] + "\n"
            fo.write(line)
            off += len(line)
        total += len(s)
open(os.path.join(out, "genome.fa.fai"), "w").write("".join(fai))
open(os.path.join(out, "contigs.txt"), "w").write("".join(c + "\n" for c in chroms))

keep, nfeat, ngene = set(chroms), 0, 0
with open(src_gff) as fi, open(os.path.join(out, "annotation.gff3"), "w") as fo:
    for l in fi:
        if l.startswith("##sequence-region"):
            if len(l.split()) > 1 and l.split()[1] in keep:
                fo.write(l)
            continue
        if l.startswith("#"):
            fo.write(l)
            continue
        f = l.split("\t", 3)
        if f[0] in keep:
            fo.write(l)
            nfeat += 1
            if len(f) > 2 and f[2] == "gene":
                ngene += 1


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
    "repo_commit:      %s%s\n"
    "contigs:          %d\n"
    "total_bp:         %d\n"
    "gff_features:     %d\n"
    "gff_genes:        %d\n"
    % (meta.get("species", "?"), datetime.datetime.now().isoformat(timespec="seconds"), sample, src_fa, src_gff,
       os.path.abspath(man), ",".join(chroms), git("rev-parse", "--short", "HEAD"),
       "  (DIRTY WORKING TREE)" if dirty else "", len(chroms), total, nfeat, ngene))
print("%s: %d chromosomes, %.1f Mb, %d GFF features (%d genes) -> %s"
      % (sample, len(chroms), total / 1e6, nfeat, ngene, out))
for c, l in zip(chroms, fai):
    print("  %-12s %7.1f Mb" % (c, int(l.split("\t")[1]) / 1e6))
