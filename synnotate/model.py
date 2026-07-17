"""Context transformer over a gene-neighbourhood window (vendored, standalone).
The deployment tracks are `ContextTransformerAdditive`: a sparse additive attention read-out
that gives an EXACT per-neighbour decomposition of the prediction (interpretable by design)."""
from __future__ import annotations
import torch
import torch.nn as nn
from synnotate.features import DIST_BINS


def sparsemax(z):
    """Sparsemax over the last dim (Martins & Astudillo 2016): most entries become exactly
    0, so the surviving attention weights ARE the explanation."""
    z = z - z.max(dim=-1, keepdim=True).values
    zs, _ = torch.sort(z, dim=-1, descending=True)
    rng = torch.arange(1, z.size(-1) + 1, device=z.device, dtype=z.dtype)
    cssv = zs.cumsum(-1) - 1
    cond = (zs - cssv / rng) > 0
    k = cond.sum(-1, keepdim=True).clamp(min=1)
    tau = cssv.gather(-1, (k - 1)) / k
    return torch.clamp(z - tau, min=0.0)


class ContextTransformerAdditive(nn.Module):
    """Same encoder as a small transformer, but the readout is an explicit sparse additive
    attention over neighbours. Because the head is linear,
        logit_c = sum_j alpha_j * (W v_j)_c + b_c
    is an exact per-neighbour decomposition — no occlusion/gradients needed."""
    def __init__(self, n_tokens, n_labels, max_len=21, dim=128, n_layers=2,
                 n_heads=4, ffn=512, dropout=0.1, n_strand=3, n_seg=3, use_sparsemax=True):
        super().__init__()
        self.use_sparsemax = use_sparsemax
        self.tok = nn.Embedding(n_tokens, dim, padding_idx=0)
        self.pos = nn.Embedding(max_len, dim)
        self.strand = nn.Embedding(n_strand, dim)
        self.dist = nn.Embedding(len(DIST_BINS), dim)
        self.seg = nn.Embedding(n_seg, dim)
        layer = nn.TransformerEncoderLayer(d_model=dim, nhead=n_heads, dim_feedforward=ffn,
                                           dropout=dropout, batch_first=True, activation="gelu")
        self.encoder = nn.TransformerEncoder(layer, num_layers=n_layers)
        self.norm = nn.LayerNorm(dim)
        self.q_proj = nn.Linear(dim, dim)
        self.k_proj = nn.Linear(dim, dim)
        self.v_proj = nn.Linear(dim, dim)
        self.scale = dim ** 0.5
        self.head = nn.Linear(dim, n_labels)
        self.dropout = nn.Dropout(dropout)

    def _encode(self, tokens, strand, dist, seg):
        B, L = tokens.shape
        pos_ids = torch.arange(L, device=tokens.device).unsqueeze(0).expand(B, L)
        x = (self.tok(tokens) + self.pos(pos_ids) + self.strand(strand)
             + self.dist(dist) + self.seg(seg))
        x = self.dropout(x)
        return self.norm(self.encoder(x, src_key_padding_mask=tokens.eq(0)))

    def _weights(self, h, tokens, target_pos):
        B, L, _ = h.shape
        idx = target_pos.view(B, 1, 1).expand(B, 1, h.size(-1))
        q = self.q_proj(h.gather(1, idx).squeeze(1))
        k = self.k_proj(h)
        scores = (q.unsqueeze(1) * k).sum(-1) / self.scale
        m = tokens.eq(0).clone()
        m.scatter_(1, target_pos.view(B, 1), True)
        scores = scores.masked_fill(m, -1e9)
        return sparsemax(scores) if self.use_sparsemax else torch.softmax(scores, -1)

    def forward(self, tokens, strand, dist, seg, target_pos, return_emb=False):
        h = self._encode(tokens, strand, dist, seg)
        alpha = self._weights(h, tokens, target_pos)
        v = self.v_proj(h)
        pooled = (alpha.unsqueeze(-1) * v).sum(1)
        if return_emb:
            return pooled
        return self.head(pooled)

    @torch.no_grad()
    def attributions(self, tokens, strand, dist, seg, target_pos):
        """Exact per-neighbour decomposition for the predicted label:
        contrib[:,j] = alpha_j*(W v_j)_pred, sum_j contrib_j + b_pred == logit_pred."""
        self.eval()
        h = self._encode(tokens, strand, dist, seg)
        alpha = self._weights(h, tokens, target_pos)
        v = self.v_proj(h)
        logits = self.head((alpha.unsqueeze(-1) * v).sum(1))
        pred = logits.argmax(-1)
        Wp = self.head.weight[pred]
        contrib = alpha * (v * Wp.unsqueeze(1)).sum(-1)
        return pred, alpha, contrib
