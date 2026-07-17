"""End-to-end Synnotate pipeline: gene-call -> annotate context -> context model -> (interpret)
-> per-gene table. The neighbour annotation vocabulary and the retrieval index match: eggNOG-mapper
(prokaryote) / Pharokka (phage) families feed the curated context model and the ctx_cur index."""
from __future__ import annotations
import os, json, subprocess, numpy as np   # `json` was MISSING: every json.load here
# (family_map on BOTH the phage and bacteria paths) raised NameError -- these branches had
# never been executed end-to-end.
from synnotate import genecall, context as ctx, interpret as itp


# ---------- annotation backends (fill the per-gene family used to build context) ----------
def annotate_from_tsv(genes, path):
    """Read a precomputed gene_id<TAB>family TSV (skip a header if present)."""
    fam = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "gene_id":
            fam[p[0]] = p[1]
    return [fam.get(g.gene_id, "UNKNOWN") for g in genes]


def _product_matcher(bundle):
    """Build a f(product_string) -> vocabulary label / OTHER / None matcher, the SAME rule the corpus
    used for RefSeq/PGAP product strings: exact `refseq_map.json`, then `canonical` token-set match,
    else OTHER (named but off-vocabulary), else None (uninformative -> dark). Shared by the GFF-product
    and eggNOG backends so both honour one matching rule."""
    from . import naming
    rpath = os.path.join(bundle.path, "refseq_map.json")
    rmap = json.load(open(rpath)) if os.path.exists(rpath) else {}     # product_raw -> label
    labels = set(bundle.label_vocab)
    canon = naming.build_canonical_index(labels)

    def match(product):
        prod = (product or "").strip()
        if not naming.is_informative(prod):
            return None                                  # no usable name -> dark
        hit = rmap.get(prod) or canon.get(naming.canonical(prod))
        return hit if (hit and hit in labels) else "OTHER"
    return match


def has_products(genes, frac=0.3):
    """True if at least `frac` of genes carry an informative GFF product string (so the input GFF is
    already annotated and can be string-matched without running eggNOG-mapper)."""
    from . import naming
    if not genes:
        return False
    n = sum(1 for g in genes if naming.is_informative((getattr(g, "product", "") or "").strip()))
    return n >= frac * len(genes)


def annotate_from_gff_products(genes, bundle):
    """Match each gene's GFF `product=` string onto the vocabulary WITHOUT running eggNOG-mapper —
    for genomes already annotated by PGAP / RefSeq / Prokka / Bakta. exact refseq_map.json -> canonical
    token-set -> OTHER (named, off-vocab) -> UNKNOWN (no informative product)."""
    match = _product_matcher(bundle)
    fams, n_rs, n_other, n_dark = [], 0, 0, 0
    for g in genes:
        lab = match(getattr(g, "product", ""))
        if lab is None:
            n_dark += 1; fams.append("UNKNOWN")
        elif lab == "OTHER":
            n_other += 1; fams.append("OTHER")
        else:
            n_rs += 1; fams.append(lab)
    print(f"[annotate] GFF-product match: named->vocab {n_rs:,} | OTHER(named, off-vocab) {n_other:,} | "
          f"dark {n_dark:,} ({n_dark/max(len(genes),1)*100:.1f}%)", flush=True)
    return fams


def annotate_pharokka(genome_fasta, workdir, db=None, threads=8, pharokka="pharokka"):
    """Run Pharokka v1.10 (phage) on the GENOME FASTA — Pharokka does its own gene-calling
    (PHANOTATE) — and map each of ITS genes to a PHROG family via `family_map.json`
    (phrog→product_norm). Returns (called_genes, families). `pharokka` may be a wrapper such as
    "micromamba run -n pharokka pharokka"."""
    out = os.path.join(workdir, "pharokka"); os.makedirs(out, exist_ok=True)
    cmd = pharokka.split() + ["run", "-i", genome_fasta, "-o", out, "-t", str(threads), "--fast", "-f"]
    if db: cmd += ["-d", db]
    subprocess.run(cmd, check=True)
    tsv = os.path.join(out, "pharokka_cds_final_merged_output.tsv")
    import csv
    from synnotate.genecall import Gene
    called, phrogs = [], []
    with open(tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gid = row.get("gene") or row.get("ID") or row.get("cds_id", "")
            called.append(Gene(gene_id=gid, contig=row.get("contig", ""), start=int(row.get("start", 0)),
                               end=int(row.get("stop", row.get("end", 0))), strand=row.get("frame", "+")[:1] or "+",
                               product=row.get("annot") or row.get("product", "")))
            phrogs.append(row.get("phrog", ""))
    return called, phrogs


# emapper column indices -- MUST match scripts/18_build_eggnog_gene_table.py:25, which built the
# corpus. QUERY, OGS, COGCAT, PNAME, KO, PFAM.
_QUERY, _OGS, _COGCAT, _PNAME, _KO, _PFAM = 0, 4, 6, 8, 11, 20


def _group_key(ko, pfam, cog):
    """KO -> 1st Pfam -> COG:<cat> -> NO_ID. Verbatim from 18_build_eggnog_gene_table.py:28."""
    if ko and ko not in ("-", ""):
        return "KO:" + ko.split(",")[0].replace("ko:", "")
    if pfam and pfam not in ("-", ""):
        return "PF:" + pfam.split(",")[0]
    if cog and cog not in ("-", ""):
        return "COG:" + cog[0]
    return "NO_ID"


def annotate_eggnog(genes, proteins_fasta, workdir, bundle, emapper="emapper.py", db=None, threads=16):
    """Label each gene the SAME way the corpus was labelled, or the context tokens are meaningless.

    Precedence (identical to publication_v2/scripts/s00_vocab_v2.py):
      1. eggNOG  -> group_key (KO->Pfam->COG->NO_ID) -> family_map.json -> label
      2. RefSeq  -> if eggNOG is silent and the GFF product names a function, match it onto the
                    vocabulary (exact via refseq_map.json, then `canonical` token-set)
      3. OTHER   -> RefSeq named it but the vocabulary has no entry (named, NOT dark)
      4. UNKNOWN -> genuinely dark: eggNOG silent AND RefSeq silent

    Three bugs this replaces, all silent:
      * it returned a BARE "K01234" while the corpus key is "KO:K01234";
      * it never applied family_map.json (that lookup existed only on the phage path), so
        token_of() fell through to OTHER for EVERY eggNOG-annotated bacterial gene -- the model
        received no context at all;
      * it never consulted the RefSeq product, so a gene RefSeq names (e.g. "cell division protein
        FtsL") was labelled UNKNOWN -- exactly the v1 corpus bug that publication_v2 exists to fix,
        re-committed on every genome a user annotates.
    """
    out = os.path.join(workdir, "eggnog"); os.makedirs(out, exist_ok=True)
    cmd = [emapper, "-i", proteins_fasta, "-o", "emap", "--output_dir", out,
           "--cpu", str(threads), "-m", "diamond", "--override"]
    if db: cmd += ["--data_dir", db]
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        raise SystemExit(
            f"[synnotate] eggNOG-mapper ('{emapper}') not found on PATH.\n"
            "  Install it (conda: `mamba install -c bioconda eggnog-mapper`) and download its DB\n"
            "  (`download_eggnog_data.py`), or — if your GFF already has product= strings — re-run\n"
            "  with `--backend gff` to string-match those products instead (no eggNOG-mapper needed).")
    ann = os.path.join(out, "emap.emapper.annotations")

    fmap = json.load(open(os.path.join(bundle.path, "family_map.json")))     # group_key -> label
    match = _product_matcher(bundle)                                         # RefSeq product -> label

    gk = {}
    for line in open(ann):
        if line.startswith("#") or line.startswith("query"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) <= _PFAM: continue
        gk[f[_QUERY]] = _group_key(f[_KO], f[_PFAM], f[_COGCAT])

    fams, n_egg, n_rs, n_other, n_dark = [], 0, 0, 0, 0
    for g in genes:
        lab = fmap.get(gk.get(g.gene_id, ""), "")
        if lab:
            n_egg += 1; fams.append(lab); continue
        hit = match(getattr(g, "product", ""))                               # eggNOG silent -> GFF product
        if hit is None:
            n_dark += 1; fams.append("UNKNOWN")
        elif hit == "OTHER":
            n_other += 1; fams.append("OTHER")   # NAMED but not in the vocabulary -> not dark
        else:
            n_rs += 1; fams.append(hit)
    print(f"[annotate] eggNOG {n_egg:,} | RefSeq {n_rs:,} | OTHER(named, off-vocab) {n_other:,} | "
          f"dark {n_dark:,} ({n_dark/max(len(genes),1)*100:.1f}%)", flush=True)
    return fams


# ---------- orchestration ----------
def run(fasta, bundle, gff=None, annotations=None, interpret=False, k=5, device="cpu",
        workdir="synnotate_work", threads=16, phage=False, backend="auto"):
    import torch
    os.makedirs(workdir, exist_ok=True)
    genes = genecall.read_gff(gff) if gff else genecall.call_genes(fasta, phage=phage)
    prot = genecall.write_proteins(genes, os.path.join(workdir, "proteins.faa"))
    # --- context families: pick the backend that fills each neighbour's family token ---
    if annotations:                                       # explicit gene_id -> family TSV
        families = annotate_from_tsv(genes, annotations)
    elif phage:
        # Pharokka does its own gene-calling (PHANOTATE) AND PHROG annotation -> use its genes,
        # map phrog -> product_norm via the bundle's family_map.json.
        genes, phrogs = annotate_pharokka(fasta, workdir, db=bundle.config.get("pharokka_db"),
                                          threads=threads, pharokka=bundle.config.get("pharokka", "pharokka"))
        fmap = json.load(open(os.path.join(bundle.path, "family_map.json")))
        families = [fmap.get(p, "UNKNOWN") for p in phrogs]
    else:
        # prokaryote: 'gff' = string-match the GFF product= column; 'eggnog' = run eggNOG-mapper;
        # 'auto' = use the GFF products when the input GFF is already annotated, else eggNOG-mapper.
        b = backend
        if b == "auto":
            b = "gff" if (gff and has_products(genes)) else "eggnog"
        print(f"[synnotate] annotation backend: {b}", flush=True)
        if b == "gff":
            families = annotate_from_gff_products(genes, bundle)
        else:
            families = annotate_eggnog(genes, prot, workdir, bundle,
                                       db=bundle.config.get("eggnog_db"), threads=threads)
    toks, strd, dist, seg, tpos, fam = ctx.build_windows(genes, families, bundle)
    # --- context model inference ---
    m = bundle.model(device); rows = []
    B = 4096
    for s in range(0, len(genes), B):
        sl = slice(s, s+B)
        T = torch.tensor(toks[sl], device=device); S = torch.tensor(strd[sl], device=device)
        D = torch.tensor(dist[sl], device=device); G = torch.tensor(seg[sl], device=device)
        P = torch.tensor(tpos[sl], device=device)
        with torch.no_grad():
            logits = m(T, S, D, G, P); prob = torch.softmax(logits, 1)
            conf, idx = prob.max(1)
            emb = m(T, S, D, G, P, return_emb=True)
            emb = torch.nn.functional.normalize(emb, dim=1).cpu().numpy()
        idx = idx.cpu().numpy(); conf = conf.cpu().numpy()
        top5 = prob.topk(5, 1).indices.cpu().numpy()
        for a, gi in enumerate(range(s, min(s+B, len(genes)))):
            g = genes[gi]
            pred = bundle.ilabel.get(int(idx[a]), "?")
            rec = {"gene_id": g.gene_id, "contig": g.contig, "start": g.start, "end": g.end,
                   "strand": g.strand, "annotation": families[gi], "prediction": pred,
                   "confidence": round(float(conf[a]), 4),
                   "top5": ";".join(bundle.ilabel.get(int(t), "?") for t in top5[a]),
                   "category": bundle.label2cat.get(pred, "")}
            if interpret:
                syn = itp.synteny(emb[a], fam[gi], g.contig, bundle, k=k)
                _, drv = itp.driving_neighbours(bundle, toks[gi], strd[gi], dist[gi], seg[gi], int(tpos[gi]), device=device)
                rec["adjusted_confidence"] = round(itp.adjusted_conf(float(conf[a]), syn["synteny_conf"]), 4)
                rec["synteny_support"] = round(syn["synteny_conf"], 4)
                rec["mean_flank_ident"] = round(syn["mean_flank_ident"], 4)
                rec["target_slot_conserved"] = round(syn["target_slot_conserved"], 4)
                rec["driving_neighbours"] = ";".join(f"{d['offset']:+d}:{d['family']}({d['contribution']:.2f})" for d in drv)
            rows.append(rec)
    return genes, rows
