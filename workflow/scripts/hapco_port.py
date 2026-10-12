#!/usr/bin/env python3
"""hapco_port.py -- hapCO_identification.R's crossover calling, line for line in Python, without the plots.

Why: a settings sweep runs the caller hundreds of times per cell. In hapCO the expensive part, the
two smoothing passes, depends only on baseAF and windowAF; the block rules (marker_num, block_size,
terminal_marker_num) act on the smoothed calls. Here a cell is smoothed once per baseAF and every block
setting reuses it, in one process, with no PDF. The calling is a transcription of the R functions,
quirks included (edge rows vote by count, interior rows by sum; pass 2 updates in place; the block
'acc' column), and is checked against the pipeline's own hapCO output on every cell before a sweep is
trusted (co_sweep.py --check).

As a module:
    rows = read_cell(path)                                  # per-cell table: chrom pos ref rc alt ac
    res = call_cell(rows, chrom_map, dict(base_af=.3, window_af=.4, genotype=.2,
                    block_size=2e6, marker_num=4, terminal_marker_num=8, cell_markers=200))
    res -> None (too few markers) or {chrom: (breakpoints, final_blocks)}
           breakpoints: [start, end, geno_pre, geno_post, acc]; final_blocks: [start, end, genotype, acc]
As a script (same arguments as hapCO_identification.R; writes <prefix>_co_pred.txt and
<prefix>_co_block_pred.txt, no PDF):
    hapco_port.py --input T --prefix BC --chrom_map M --outpath DIR [--cell_markers 200 --block_size 2e6
                  --marker_num 8 --baseAF 0.3 --windowAF 0.4 --genotype 0.2 --terminal_marker_num 0]
"""
import argparse
import csv
import math
import os


def read_cell(path):
    """Rows in file order: (chrom, pos, ref_count, alt_count)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                rows.append((f[0], int(float(f[1])), int(float(f[3])), int(float(f[5]))))
    return rows


def read_chrom_map(path):
    with open(path) as fh:
        return [(r["name"], int(float(r["size"]))) for r in csv.DictReader(fh, delimiter="\t")]


def smooth(raw, base_af, window_af):
    """allele_count_smoother: raw = REF share per row (column 7); returns column 8."""
    n = len(raw)
    smt = [0.0] * n
    for i in range(n):                                   # pass 1: votes from the RAW neighbours
        lo, hi = max(i - 2, 0), min(i + 2, n - 1)
        den = hi - lo
        if i == 0:
            win = sum(1 for j in range(1, hi + 1) if raw[j] > 0) / den
        elif i == n - 1:
            win = sum(1 for j in range(lo, i) if raw[j] > 0) / den
        else:
            win = (math.fsum(raw[lo:i]) + math.fsum(raw[i + 1:hi + 1])) / den
        if win > 1 - base_af:
            smt[i] = 1.0
        elif win < base_af:
            smt[i] = 0.0
        else:
            smt[i] = raw[i]
    for i in range(n):                                   # pass 2: in place, so left neighbours are new
        lo, hi = max(i - 2, 0), min(i + 2, n - 1)
        den = hi - lo
        if i == 0:
            win = sum(1 for j in range(1, hi + 1) if smt[j] > 0) / den
        elif i == n - 1:
            win = sum(1 for j in range(lo, i) if smt[j] > 0) / den
        else:
            win = (sum(1 for j in range(lo, i) if smt[j] > 0) + sum(1 for j in range(i + 1, hi + 1) if smt[j] > 0)) / den
        if window_af <= win <= 1 - window_af:
            smt[i] = raw[i]
    return smt


def _geno(v, genotype):
    return 1 if v >= 1 - genotype else (0 if v <= genotype else -1)


def genotype_blocks(pos, smt, genotype):
    """get_genotype_block: runs of identical smoothed values -> [sta, end, genotype, markernum]."""
    cur = [pos[0], pos[0], _geno(smt[0], genotype), 1]
    blocks = []
    for i in range(1, len(smt)):
        if smt[i] != smt[i - 1]:
            blocks.append(cur)
            cur = [pos[i], pos[i], _geno(smt[i], genotype), 1]
        else:
            cur[1] = pos[i]
            cur[3] += 1
    blocks.append(cur)
    return blocks


def filter_blocks(blocks, block_size, marker_num, terminal):
    nb, out = len(blocks), []
    for r, b in enumerate(blocks):
        end_block = terminal > 0 and (r == 0 or r == nb - 1)
        keep = b[3] >= terminal if end_block else (b[1] - b[0] >= block_size and b[3] >= marker_num)
        if keep:
            out.append(list(b))
    return out


def get_breakpoints(filtered):
    bps = []
    if len(filtered) > 1:
        acc = filtered[0][3]
        for r in range(1, len(filtered)):
            if filtered[r][2] != filtered[r - 1][2]:
                bps.append([filtered[r - 1][1], filtered[r][0], filtered[r - 1][2], filtered[r][2], acc])
                acc = filtered[r][3]
            else:
                acc += filtered[r][3]
    return bps


def fine_breakpoints(bps, pos, smt):
    out = []
    for b in bps:
        idx = [i for i, p in enumerate(pos) if b[0] <= p <= b[1]]
        seq = [smt[i] for i in idx]
        L = len(seq)
        nb = list(b)
        if L >= 2:
            z = [0] * (L + 1)
            o = [0] * (L + 1)
            for k, v in enumerate(seq):
                z[k + 1] = z[k] + (v == 0)
                o[k + 1] = o[k] + (v == 1)
            max_score, max_pos = 0, 1
            for p in range(1, L):
                af0L, af1L = z[p] / p, o[p] / p
                af0R, af1R = (z[L] - z[p]) / (L - p), (o[L] - o[p]) / (L - p)
                score = max(af0L * af1R, af1L * af0R)
                if score > max_score:
                    max_score, max_pos = score, p
            nb[0], nb[1] = pos[idx[max_pos - 1]], pos[idx[max_pos]]
        out.append(nb)
    return out


def make_final_blocks(filtered, chr_size):
    final, cur, acc, last = [], [1, 1, None, None], 0, 0
    if len(filtered) > 1:
        for r in range(1, len(filtered)):
            if filtered[r][2] != filtered[r - 1][2]:
                cur[1], cur[2], cur[3] = filtered[r - 1][1], filtered[r - 1][2], acc
                final.append(list(cur))
                cur[0], cur[1] = filtered[r][0], filtered[r][1]
            else:
                acc += filtered[r][3]
        last = len(filtered) - 1
    if filtered:
        cur[1], cur[2], cur[3] = filtered[last][1], filtered[last][2], acc
        final.append(list(cur))
        if final[-1][1] < chr_size:
            final[-1][1] = chr_size
    return final


def fine_blocks(bps, final):
    for k in range(1, len(final)):
        final[k - 1][1] = bps[k - 1][0]
        final[k][0] = bps[k - 1][1]
    return final


def blocks_to_calls(pos, smt, chr_size, genotype, block_size, marker_num, terminal):
    """From one chromosome's smoothed calls to (breakpoints, final blocks); None if no block passes."""
    filtered = filter_blocks(genotype_blocks(pos, smt, genotype), block_size, marker_num, terminal)
    if not filtered:
        return None
    bps = get_breakpoints(filtered)
    bpu = fine_breakpoints(bps, pos, smt) if bps else []
    final = make_final_blocks(filtered, chr_size)
    if bps:
        final = fine_blocks(bpu, final)
    return bpu, final


def prepare(rows, chrom_map, cell_markers):
    """Informative rows per chromosome (file order) as (pos, raw REF share); None below the cell floor."""
    inf = [r for r in rows if r[2] > 0 or r[3] > 0]
    if len(inf) < cell_markers:
        return None
    per = {}
    for c, p, rc, ac in inf:
        per.setdefault(c, []).append((p, rc / (rc + ac)))
    return {c: per.get(c, []) for c, _ in chrom_map}


def call_cell(rows, chrom_map, prm, smoothed=None):
    """One cell at one setting. 'smoothed' may carry {chrom: (pos, smt)} from smooth_cell for this baseAF."""
    prep = prepare(rows, chrom_map, prm["cell_markers"]) if smoothed is None else None
    if smoothed is None and prep is None:
        return None
    if smoothed is None:
        smoothed = smooth_cell(prep, prm["base_af"], prm["window_af"])
    out = {}
    for c, size in chrom_map:
        if c not in smoothed:
            continue
        pos, smt = smoothed[c]
        r = blocks_to_calls(pos, smt, size, prm["genotype"], prm["block_size"], prm["marker_num"],
                            prm["terminal_marker_num"])
        if r is not None:
            out[c] = r
    return out


def smooth_cell(prep, base_af, window_af):
    sm = {}
    for c, x in prep.items():
        if len(x) < 2:
            continue                                     # hapCO: 'too few markers' -> nothing written
        pos = [p for p, _ in x]
        sm[c] = (pos, smooth([a for _, a in x], base_af, window_af))
    return sm


def _r_num(v):
    v = float(v)
    return str(int(v)) if v == int(v) else repr(v)


def write_like_hapco(res, outpath, prefix):
    """The two text files hapCO writes (numbers as plain integers; hapCO's R may print 1e+08 style)."""
    bp_out = os.path.join(outpath, prefix + "_co_pred.txt")
    bk_out = os.path.join(outpath, prefix + "_co_block_pred.txt")
    for p in (bp_out, bk_out):
        if os.path.exists(p):
            os.remove(p)
    if res is None:
        return
    with open(bp_out, "a") as fb, open(bk_out, "a") as fk:
        for c, (bps, final) in res.items():
            for b in bps:
                fb.write(" ".join([c] + [_r_num(x) for x in b]) + "\n")
            for b in final:
                fk.write(" ".join([c] + [_r_num(x) for x in b]) + "\n")
    for p in (bp_out, bk_out):
        if os.path.exists(p) and os.path.getsize(p) == 0:
            os.remove(p)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", "-i", required=True); ap.add_argument("--prefix", "-p", default="out")
    ap.add_argument("--chrom_map", "-g", required=True); ap.add_argument("--outpath", "-o", default=".")
    ap.add_argument("--cell_markers", "-c", type=int, default=200)
    ap.add_argument("--block_size", "-s", type=float, default=2000000)
    ap.add_argument("--marker_num", "-n", type=int, default=8)
    ap.add_argument("--baseAF", type=float, default=0.3); ap.add_argument("--windowAF", type=float, default=0.4)
    ap.add_argument("--genotype", type=float, default=0.2)
    ap.add_argument("--terminal_marker_num", type=int, default=0)
    A = ap.parse_args()
    os.makedirs(A.outpath, exist_ok=True)
    res = call_cell(read_cell(A.input), read_chrom_map(A.chrom_map),
                    dict(base_af=A.baseAF, window_af=A.windowAF, genotype=A.genotype, block_size=A.block_size,
                         marker_num=A.marker_num, terminal_marker_num=A.terminal_marker_num,
                         cell_markers=A.cell_markers))
    write_like_hapco(res, A.outpath, A.prefix)
