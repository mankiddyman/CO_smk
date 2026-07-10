#!/usr/bin/env python3
"""
Preprocess the Spondias hap1 Helixer GTF into a STAR-ready GTF.

Two fixes:
  1. Seqnames: 'scaffold_N' -> 'N'  (to match the assembly FASTA headers >1, >2, ...)
  2. transcript_id: the Helixer GTF put a functional DESCRIPTION in transcript_id
     (non-unique). The real transcript identifier depends on feature type:
       - mRNA line:          its own ID attr is the transcript id
       - exon/CDS/UTR lines: the Parent attr names the transcript they belong to
         (their own ID is feature-level, e.g. "...1.exon.4", which must NOT be used)

gene_id is already clean and unique (hap1_scaffold_1_003111) and is left untouched.

Usage:
  python3 fix_spondias_gtf.py IN.gtf OUT.gtf
"""
import sys
import re

def get_attr(attr_str, key):
    """Return the value of key "..." from a GTF attribute string, or None."""
    m = re.search(rf'{key} "([^"]*)"', attr_str)
    return m.group(1) if m else None

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: fix_spondias_gtf.py IN.gtf OUT.gtf")
    inp, outp = sys.argv[1], sys.argv[2]

    n_lines = 0
    n_seqname_fixed = 0
    n_tid_fixed = 0
    n_tid_missing_id = 0
    seqnames_before = set()
    seqnames_after = set()

    with open(inp) as fh, open(outp, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                out.write(line)
                continue

            # Fix 1: seqname scaffold_N -> N
            seqnames_before.add(parts[0])
            if parts[0].startswith("scaffold_"):
                parts[0] = parts[0][len("scaffold_"):]
                n_seqname_fixed += 1
            seqnames_after.add(parts[0])

            # Fix 2: set transcript_id to the correct TRANSCRIPT identifier.
            #   - mRNA line:            its own ID attr IS the transcript id
            #   - exon/CDS/UTR lines:   the Parent attr is the transcript they belong to
            feature = parts[2]
            attr = parts[8]
            tid_val = get_attr(attr, "transcript_id")
            if tid_val is not None:
                if feature == "mRNA":
                    new_tid = get_attr(attr, "ID")
                else:
                    new_tid = get_attr(attr, "Parent")
                if new_tid is None:
                    gid = get_attr(attr, "gene_id")
                    new_tid = f"{gid}.1" if gid else tid_val
                    n_tid_missing_id += 1
                attr = re.sub(r'transcript_id "[^"]*"',
                              f'transcript_id "{new_tid}"', attr)
                n_tid_fixed += 1
            parts[8] = attr

            out.write("\t".join(parts) + "\n")
            n_lines += 1

    sys.stderr.write(f"Lines processed:            {n_lines}\n")
    sys.stderr.write(f"Seqnames renamed:           {n_seqname_fixed}\n")
    sys.stderr.write(f"transcript_id set:          {n_tid_fixed}\n")
    sys.stderr.write(f"transcript_id fallback:     {n_tid_missing_id}\n")
    sys.stderr.write(f"Seqnames before: {sorted(seqnames_before)[:6]}...\n")
    sys.stderr.write(f"Seqnames after:  {sorted(seqnames_after)[:6]}...\n")

if __name__ == "__main__":
    main()
