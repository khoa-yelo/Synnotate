"""Load a Synnotate interpretation bundle (models + vocab + config + FAISS index + loci).

Bundle layout (per organism type, produced by `synnotate setup`):
  config.json            DIM/layers/heads/W, calibration temp, organism, vocab special tokens
  curated_context.pt     ContextTransformerAdditive on curated neighbour tokens
  vocab.json             family -> context-token id
  label_vocab.json       family -> label id
  lab2tok.npy            label id -> context-token id
  seg.npy                segment vector
  label2cat.json         family -> functional category (PHROG / eggNOG)
  index/ctx_cur.faiss    curated-context retrieval index
  index/meta.parquet     row -> (genome_id, gene_id, label_id)
  index/loci_tokens.npy  curated +/-W family windows (row-aligned, for kNN + MSA)
  calibration.json       (optional) deploy temperature + isotonic map + trusted-region contours
"""
from __future__ import annotations
import json, os
import numpy as np


class Bundle:
    def __init__(self, path):
        self.path = path
        self.config = json.load(open(os.path.join(path, "config.json")))
        self.vocab = json.load(open(os.path.join(path, "vocab.json")))          # family -> token id
        self.label_vocab = json.load(open(os.path.join(path, "label_vocab.json")))
        self.ivocab = {v: k for k, v in self.vocab.items()}
        self.ilabel = {v: k for k, v in self.label_vocab.items()}
        self.label2cat = self._maybe_json("label2cat.json")
        self.W = int(self.config.get("W", 10))
        self.NT = len(self.vocab); self.NL = len(self.label_vocab)
        for k in ("PAD", "HYP_TARGET", "CONTIG_EDGE", "OTHER", "UNKNOWN"):
            setattr(self, k, self.vocab[k])
        self.seg = np.load(os.path.join(path, "seg.npy")) if os.path.exists(os.path.join(path, "seg.npy")) \
            else np.array([0]*self.W + [1] + [2]*self.W, dtype=np.int64)
        # ---- calibration (optional): temperature + isotonic map + trusted-region contours -------
        # Without it the tool degrades gracefully to raw softmax confidence and no trust flag.
        self.calib = self._maybe_json("calibration.json")
        self.deploy_temp = float(self.calib.get("deploy_temp", 1.0))
        self._iso_x = np.asarray(self.calib.get("isotonic_x", []), dtype=float)
        self._iso_y = np.asarray(self.calib.get("isotonic_y", []), dtype=float)
        self.trusted_region = self.calib.get("trusted_region", {})   # {"0.95":{a,p,q}, "0.99":{...}}
        # lazy heavy assets
        self._model = None; self._index = None; self._meta = None; self._loci = None

    def calibrate(self, conf):
        """Map a temperature-scaled max-probability to a calibrated expected accuracy via the
        isotonic map. Identity if the bundle ships no isotonic map."""
        if self._iso_x.size == 0:
            return float(conf)
        return float(np.interp(float(conf), self._iso_x, self._iso_y))

    def trust_level(self, conf, syn):
        """Strictest expected-accuracy threshold whose trusted-region contour the prediction meets,
        as a string ('0.99' > '0.95'); '' if it meets none / no contours are shipped. The contour is
        `syn >= a - p*conf - q*conf^2` (clipped to [0,1]), fit on calibrated confidence x synteny."""
        best = ""
        for thr in sorted(self.trusted_region):                       # '0.95' before '0.99'
            r = self.trusted_region[thr]
            bound = min(max(r["a"] - r["p"]*conf - r["q"]*conf*conf, 0.0), 1.0)
            if syn >= bound:
                best = thr
        return best

    def _maybe_json(self, name):
        p = os.path.join(self.path, name)
        return json.load(open(p)) if os.path.exists(p) else {}

    # ---- context model (curated track) ----
    def model(self, device="cpu"):
        if self._model is None:
            import torch
            from synnotate.model import ContextTransformerAdditive
            c = self.config
            m = ContextTransformerAdditive(self.NT, self.NL, max_len=2*self.W+1,
                                           dim=c.get("dim", 128), n_layers=c.get("n_layers", 2),
                                           n_heads=c.get("n_heads", 4), ffn=c.get("ffn", 512),
                                           use_sparsemax=True).to(device).eval()
            m.load_state_dict(torch.load(os.path.join(self.path, "curated_context.pt"), map_location=device))
            self._model = m
        return self._model

    # ---- retrieval index + row-aligned loci windows (for kNN + MSA) ----
    def index(self):
        if self._index is None:
            import faiss
            self._index = faiss.read_index(os.path.join(self.path, "index", "ctx_cur.faiss"))
        return self._index

    def meta(self):
        if self._meta is None:
            import pandas as pd
            self._meta = pd.read_parquet(os.path.join(self.path, "index", "meta.parquet"))
        return self._meta

    def loci(self):
        if self._loci is None:
            self._loci = np.load(os.path.join(self.path, "index", "loci_tokens.npy"))
        return self._loci

    def token_of(self, family):
        """Map an annotation family string to its context-token id (OTHER if unseen, UNKNOWN if dark)."""
        if family in (None, "", "UNKNOWN", "hypothetical protein"):
            return self.UNKNOWN
        return self.vocab.get(family, self.OTHER)
