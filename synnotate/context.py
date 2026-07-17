"""Build the +/-W genomic-context windows Synnotate consumes — one masked window per gene,
exactly as the model was trained: target -> HYP_TARGET (function hidden), neighbours carry their
family token + strand + intergenic-distance bin; contig edges -> CONTIG_EDGE."""
from __future__ import annotations
import numpy as np
from synnotate.features import distance_bin, STRAND_IDX


def build_windows(genes, families, bundle):
    """genes: list[Gene] sorted by (contig,start); families: parallel list of family strings
    (curated neighbour annotation). Returns arrays (tokens, strand, dist, seg, target_pos)
    each of shape (n_genes, 2W+1) [target_pos is (n_genes,)] plus the per-gene family window used
    for MSA (centre = HYP)."""
    W = bundle.W; L = 2*W+1
    HYP, EDGE = bundle.HYP_TARGET, bundle.CONTIG_EDGE
    n = len(genes)
    toks = np.zeros((n, L), np.int64); strd = np.full((n, L), STRAND_IDX["."], np.int64)
    dist = np.zeros((n, L), np.int64); seg = np.tile(bundle.seg, (n, 1)).astype(np.int64)
    fam = np.zeros((n, L), np.int64)     # family window for MSA (same as toks but centre stays family-less)
    # contig boundaries
    contig = [g.contig for g in genes]
    lo = np.zeros(n, int); hi = np.zeros(n, int); s = 0
    for i in range(1, n+1):
        if i == n or contig[i] != contig[s]:
            lo[s:i] = s; hi[s:i] = i; s = i
    tok_of = [bundle.token_of(f) for f in families]
    starts = [g.start for g in genes]; ends = [g.end for g in genes]; strands = [g.strand for g in genes]
    for i in range(n):
        for off in range(-W, W+1):
            k = off + W; j = i + off
            if off == 0:
                toks[i, k] = HYP; fam[i, k] = HYP; strd[i, k] = STRAND_IDX.get(strands[i], 2)
                lg = (starts[i]-ends[i-1]-1) if i-1 >= lo[i] else None
                rg = (starts[i+1]-ends[i]-1) if i+1 < hi[i] else None
                g = [x for x in (lg, rg) if x is not None]
                dist[i, k] = distance_bin(min(g)) if g else distance_bin(None)
                continue
            if j < lo[i] or j >= hi[i]:
                toks[i, k] = EDGE; fam[i, k] = EDGE; dist[i, k] = distance_bin(None); continue
            toks[i, k] = tok_of[j]; fam[i, k] = tok_of[j]; strd[i, k] = STRAND_IDX.get(strands[j], 2)
            gap = (starts[j+1]-ends[j]-1) if off < 0 else (starts[j]-ends[j-1]-1)
            dist[i, k] = distance_bin(int(gap))
    target_pos = np.full(n, W, np.int64)
    return toks, strd, dist, seg, target_pos, fam
