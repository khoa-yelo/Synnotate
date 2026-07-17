# Synnotate

**Synteny-transformer gene-function annotation from genomic context.**

Given a genome (FASTA, or a GFF of CDS calls), Synnotate builds each gene's ±10-gene
neighbourhood as a sequence of protein-family tokens and runs an interpretable transformer to
predict the gene's function — **including the dark proteome** (hypothetical proteins with no
sequence homolog), which it annotates from *where the gene sits* rather than *what it looks like*.

Every prediction comes with:

1. a **softmax confidence**,
2. an optional **kNN + MSA synteny-adjusted confidence** — does the same neighbourhood recur in
   other genomes, with the predicted family in the same slot?, and
3. an optional **exact per-neighbour attribution** — which neighbouring genes drove the call. The
   additive sparsemax read-out makes this decomposition faithful by construction, not a post-hoc
   approximation.

---

## Installation

```bash
pip install git+https://github.com/khoa-yelo/Synnotate.git
# or, from a local clone:
git clone https://github.com/khoa-yelo/Synnotate.git && pip install ./Synnotate
```

Python ≥ 3.10. Core dependencies (numpy, torch, faiss-cpu, pandas, pyarrow, biopython) install
automatically. Gene calling uses `pyrodigal` / `pyrodigal-gv` in-process.

### Download the interpretation bundle

The model, vocabulary, and retrieval index are downloaded separately (too large for a Python wheel).
Pick a directory and download the prokaryote bundle into it:

```bash
synnotate setup --type prokaryote --dir ./synnotate_db
```

This fetches a versioned tarball, verifies its checksum, and extracts it into `./synnotate_db`. Pass
that directory to every `annotate` call with `--bundle ./synnotate_db`. Bundles are published on
Zenodo (URL/checksum in the release notes).

Other options:

```bash
synnotate setup --type prokaryote --dir ./synnotate_db --from ./bundle   # install a local dir or .tar.gz
synnotate setup --type prokaryote --dir ./synnotate_db --check           # verify the download (sha256)
```

### Try it on the example data

The repository ships a small example under [`examples/`](examples/): a division/cell-wall (`dcw`)
gene cluster where `ftsL` has been blanked to `hypothetical protein`. From the neighbouring gene
families alone, Synnotate recovers it:

```bash
synnotate annotate examples/demo.fna --type prokaryote --bundle ./synnotate_db \
    --gff examples/demo.gff --backend gff --interpret --out demo
```

```
gene_id  annotation  prediction                   confidence  synteny_support
dcw_03   UNKNOWN     cell division protein FtsL    0.55        0.63
```

(driven by its neighbours: `-1 RsmH`, `+1 FtsI`, `+2 MurF`.) See [`examples/README.md`](examples/README.md).

---

## Neighbour annotation (an external step)

Synnotate reads the **family labels of a gene's neighbours**, so a new genome must first have those
neighbours named. Pick the cheapest option that fits your input:

| Your input | Use | External tool needed |
|---|---|---|
| A GFF already annotated (PGAP / RefSeq / Prokka / Bakta) | `--backend gff` | **none** — product strings are string-matched onto the vocabulary |
| Your own `gene_id → family` table | `--annotations fams.tsv` | **none** |
| An un-annotated genome | `--backend eggnog` (prokaryote) | [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) + its DB |
| A phage genome | `--type phage` | [Pharokka](https://github.com/gbouras13/pharokka) |

Product-string matching uses the **same rule the training corpus used**: exact match, then a
`canonical` token-set match (word-order-insensitive, e.g.
`N-acetyltransferase family GNAT` → `GNAT family N-acetyltransferase`), else `OTHER`
(named but off-vocabulary), else dark.

To install eggNOG-mapper (only for `--backend eggnog`):
```bash
mamba install -c bioconda eggnog-mapper
download_eggnog_data.py            # downloads its ~50 GB database
```

---

## Quick start

```bash
# 1. an already-annotated genome — no external annotation tool, add synteny + attributions
synnotate annotate genome.fna --type prokaryote --bundle ./synnotate_db --gff pgap.gff --interpret --out result

# 2. an un-annotated genome — run eggNOG-mapper to name neighbours, then predict
synnotate annotate genome.fna --type prokaryote --bundle ./synnotate_db --backend eggnog --out result

# 3. a phage genome
synnotate annotate phage.fna --type phage --bundle ./synnotate_db --interpret --out result

# 4. bring your own neighbour families
synnotate annotate genome.fna --type prokaryote --bundle ./synnotate_db --annotations fams.tsv --out result
```

Outputs `result.synnotate.tsv` (and, with `--gff-out`, an annotated GFF).

### Key options

| Option | Meaning |
|---|---|
| `--type {prokaryote,phage}` | selects the bundle and default backend |
| `--backend {auto,eggnog,gff}` | prokaryote neighbour source; `auto` uses GFF products if present, else eggNOG-mapper |
| `--gff FILE` | use these CDS coordinates (and `product=` strings) instead of calling genes |
| `--annotations TSV` | `gene_id⇥family` table; skip the backend entirely |
| `--interpret` | add kNN(k=5)+MSA synteny-adjusted confidence and per-neighbour attributions |
| `--k N` | neighbours retrieved for synteny (default 5) |
| `--device {auto,cpu,cuda}` | inference device (default auto) |
| `--out PREFIX` | output prefix (default `synnotate_out`) |
| `--gff-out FILE` | also write an annotated GFF |

---

## Output

One row per gene in `<prefix>.synnotate.tsv`:

| column | meaning |
|---|---|
| `gene_id, contig, start, end, strand` | gene coordinates |
| `annotation` | the neighbour-annotation family fed to the model (`UNKNOWN` = dark) |
| `prediction` | Synnotate's predicted family |
| `confidence` | softmax confidence of the prediction |
| `category` | functional category of the prediction |
| `top5` | top-5 predicted families |

With `--interpret`, five more columns:

| column | meaning |
|---|---|
| `adjusted_confidence` | confidence combined with synteny support |
| `synteny_support` | flank-identity × target-slot-conservation over the k retrieved loci |
| `mean_flank_ident` | mean neighbourhood-token identity across retrieved loci |
| `target_slot_conserved` | fraction of retrieved loci placing the predicted family in the central slot |
| `driving_neighbours` | per-neighbour attribution, e.g. `+1:lipoprotein LipL21(4.05);-2:...` |

A dark gene is well-supported when it is **both** context-confident **and** syntenically supported —
a two-dimensional criterion, not a single recalibrated probability.

---

## How it works

- **Input.** A ±10-gene window (21 positions) of family tokens centred on the target; the target's
  own token is masked. Each position also carries strand, intergenic-distance, and segment embeddings.
- **Model.** A 2-layer transformer encoder feeding an **additive sparsemax attention head**: it pools
  neighbours with one explicit, sparse attention distribution and a single linear layer, so each
  neighbour's contribution to the logit is exactly recoverable.
- **Dark proteome.** Because the prediction reads *context*, not sequence, it annotates genes with no
  detectable homolog — the ones homology-based tools leave as "hypothetical protein".
- **Synteny corroboration.** A FAISS index of context embeddings retrieves the k nearest
  neighbourhoods from other genomes; a Needleman–Wunsch alignment of the family-token windows (both
  orientations) scores how conserved the neighbourhood is and whether it consistently hosts the
  predicted family.

---

## Reproducing / building a bundle

Bundles are produced by the training pipeline (not in this repo). A bundle directory contains
`config.json`, `curated_context.pt`, the vocabulary + `family_map.json` / `refseq_map.json`,
`index/{ctx_cur.faiss, meta.parquet, loci_tokens.npy, loci_strand.npy}`, a `MODEL_CARD.md`, and a
sha256 `MANIFEST.json`. Point Synnotate at any bundle with `--bundle DIR`.

## Citation

If you use Synnotate, please cite the accompanying paper (in preparation). See `MODEL_CARD.md` in the
bundle for the vocabulary-construction rule and dataset provenance.

## License

MIT — see [LICENSE](LICENSE).
