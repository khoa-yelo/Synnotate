"""Interpretation: kNN retrieval + MSA synteny verification + exact additive attribution.

For a masked target gene's +/-W family window, retrieve the k most similar training contexts,
Needleman-Wunsch align the query family window against each, and score:
  mean_flank_ident       = avg fraction of matching flank families across the k loci
  target_slot_conserved  = fraction of loci whose CENTRE slot the query centre aligns to
  synteny_conf           = mean_flank_ident * target_slot_conserved
The MSA-adjusted confidence is geom-mean(softmax_conf, synteny_conf). The additive attribution
comes straight from the sparsemax read-out (contrib_j + bias == logit_pred, exact)."""
from __future__ import annotations
import numpy as np


def nw(Q, R, W, SPECIAL, hyp, edge, gap=-1.0):
    """Needleman-Wunsch over two token windows. Returns (flank_identity, aligned_centre_token).
    Substitution scoring matches the deployed scorer (scripts/fig/synteny_k5_bact.py): identical
    non-reserved tokens +2, the masked-target (hyp) or contig-edge (edge) positions 0 (neutral),
    every other pair -1; gaps -1. Aligning these exactly keeps the synteny score on the same scale
    the trusted-region contours were fit on."""
    n, m = len(Q), len(R)
    S = np.full((n, m), -1.0)
    for i in range(n):
        qi = Q[i]
        for j in range(m):
            rj = R[j]
            if qi == hyp or rj == hyp or qi == edge or rj == edge:
                S[i, j] = 0.0
            elif qi == rj and qi not in SPECIAL:
                S[i, j] = 2.0
    H = np.zeros((n+1, m+1)); P = np.zeros((n+1, m+1), np.int8)
    for i in range(1, n+1): H[i, 0] = i*gap; P[i, 0] = 1
    for j in range(1, m+1): H[0, j] = j*gap; P[0, j] = 2
    for i in range(1, n+1):
        Hi, Hp, Si, Pi = H[i], H[i-1], S[i-1], P[i]
        for j in range(1, m+1):
            d = Hp[j-1] + Si[j-1]; u = Hp[j] + gap; l = Hi[j-1] + gap
            if d >= u and d >= l: Hi[j] = d; Pi[j] = 0
            elif u >= l: Hi[j] = u; Pi[j] = 1
            else: Hi[j] = l; Pi[j] = 2
    i, j, ctr, aR, match, tot = n, m, W, None, 0, 0
    while i > 0 or j > 0:
        p = P[i, j]
        if p == 0:
            qi, rj = i-1, j-1
            if qi == ctr: aR = R[rj]
            elif Q[qi] not in SPECIAL and R[rj] not in SPECIAL: tot += 1; match += (Q[qi] == R[rj])
            i -= 1; j -= 1
        elif p == 1: i -= 1
        else: j -= 1
    return (match/tot if tot else 0.0), aR


def synteny(query_emb, query_fam, query_genome, bundle, k=5):
    """kNN + MSA for a single query. query_emb: (D,) normalised context embedding; query_fam:
    the +/-W curated/predicted family window (centre = HYP). Returns a dict of the components."""
    import faiss  # noqa: F401
    W, SPECIAL = bundle.W, {bundle.PAD, bundle.HYP_TARGET, bundle.CONTIG_EDGE, bundle.OTHER, bundle.UNKNOWN}
    meta = bundle.meta(); loci = bundle.loci()
    mg = meta.genome_id.values; ml = meta.label_id.values
    D, I = bundle.index().search(query_emb[None, :].astype(np.float32), k + 16)
    idents, slots, hit_rows = [], [], []
    seen = {query_genome}
    for jj in I[0]:
        if mg[jj] in seen:
            continue
        seen.add(mg[jj]); hit_rows.append(jj)
        R = loci[jj].astype(np.int64).copy()
        R[W] = bundle.vocab.get(bundle.ilabel.get(int(ml[jj])), bundle.OTHER)   # unmask retrieved centre
        # some retrieved loci are reverse-complemented (reversed gene order) vs the query -> try both
        # orientations and keep the better-aligning one (window is odd-length, so centre is invariant)
        ident, aR = nw(query_fam, R, W, SPECIAL, bundle.HYP_TARGET, bundle.CONTIG_EDGE)
        ident_r, aR_r = nw(query_fam, R[::-1].copy(), W, SPECIAL, bundle.HYP_TARGET, bundle.CONTIG_EDGE)
        if ident_r > ident:
            ident, aR = ident_r, aR_r
        idents.append(ident)
        slots.append(1.0 if (aR is not None and int(aR) == R[W]) else 0.0)
        if len(hit_rows) >= k:
            break
    mfi = float(np.mean(idents)) if idents else 0.0
    slot = float(np.mean(slots)) if slots else 0.0
    return {"n_hits": len(hit_rows), "mean_flank_ident": mfi, "target_slot_conserved": slot,
            "synteny_conf": mfi*slot,
            "support_genomes": [str(mg[r]) for r in hit_rows],
            "support_families": [bundle.ilabel.get(int(ml[r]), "?") for r in hit_rows]}


def adjusted_conf(softmax_conf, synteny_conf):
    """MSA-adjusted confidence = geometric mean of softmax and synteny confidence."""
    return float(np.sqrt(max(softmax_conf, 0.0) * max(synteny_conf, 0.0)))


def driving_neighbours(bundle, tokens, strand, dist, seg, target_pos, top=3, device="cpu"):
    """Exact per-neighbour attribution for the predicted label (sparsemax additive read-out).
    Returns the top driving neighbour positions (offset from target) with their contributions."""
    import torch
    m = bundle.model(device)
    T = torch.tensor(tokens[None, :], device=device); S = torch.tensor(strand[None, :], device=device)
    Dd = torch.tensor(dist[None, :], device=device); G = torch.tensor(seg[None, :], device=device)
    P = torch.tensor([target_pos], device=device)
    pred, alpha, contrib = m.attributions(T, S, Dd, G, P)
    c = contrib[0].cpu().numpy(); order = np.argsort(-np.abs(c))
    out = []
    for pos in order:
        if pos == target_pos or tokens[pos] == bundle.PAD:
            continue
        out.append({"offset": int(pos - target_pos), "family": bundle.ivocab.get(int(tokens[pos]), "?"),
                    "contribution": float(c[pos])})
        if len(out) >= top:
            break
    return int(pred[0]), out
