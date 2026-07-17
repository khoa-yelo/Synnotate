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

The model, vocabulary, and retrieval index ship as a separate **bundle** (they are too large for a
Python wheel). `synnotate setup` downloads and verifies it, exactly like Pharokka's
`install_databases.py`:

```bash
synnotate setup --type prokaryote                 # download to ~/.synnotate/bundle/prokaryote
synnotate setup --type prokaryote --dir /data/db  # ... or a location you choose
synnotate setup --type prokaryote --from ./bundle # ... or install a local copy / .tar.gz
synnotate setup --type prokaryote --check         # verify an installed bundle (sha256 MANIFEST)
```

`setup` fetches a versioned tarball, md5-checks it, extracts, and verifies a per-file sha256
manifest. The bundle URL/md5 live in `synnotate/cli.py::BUNDLE_REGISTRY`, or override with
`$SYNNOTATE_BUNDLE_URL_PROKARYOTE`. Bundles are published on Zenodo (DOI in the release notes).

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
synnotate annotate genome.fna --type prokaryote --gff pgap.gff --interpret --out result

# 2. an un-annotated genome — run eggNOG-mapper to name neighbours, then predict
synnotate annotate genome.fna --type prokaryote --backend eggnog --out result

# 3. a phage genome (Pharokka backend)
synnotate annotate phage.fna --type phage --interpret --out result

# 4. bring your own neighbour families
synnotate annotate genome.fna --type prokaryote --annotations fams.tsv --out result
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
