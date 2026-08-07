# Synnotate

**Genomic-context annotation of the prokaryotic and viral dark proteome.**

Synnotate is an interpretable framework for annotating uncharacterized microbial and viral proteins
from genomic context, beyond sequence and structural homology. It represents a gene's neighbourhood
as standardized functional-annotation tokens and infers the target gene's function with a
Transformer, reporting a calibrated confidence, exact per-neighbour attribution, and retrieved
reference neighbourhoods.

![Synnotate overview](https://raw.githubusercontent.com/khoa-yelo/Synnotate/main/docs/synnotate_overview.png)

> **(A) Training** — a masked-token Transformer learns to predict the central gene of a 21-gene window
> from its neighbours' standardized annotations plus strand, intergenic distance, and relative position.
> **(B) Inference** — for a hypothetical protein it returns a prediction with a calibrated confidence,
> per-neighbour attribution, and the *k* nearest reference neighbourhoods from a vector database.
> **(C) Post-hoc alignment** — a synteny score (target alignment × flank identity) corroborates each call.

## Install

```bash
pip install synnotate                                         # from PyPI
# or the latest development version from GitHub:
pip install git+https://github.com/khoa-yelo/Synnotate.git
```

**With conda (recommended on clusters, or if `pip` tries to build numpy/torch from source).**
Synnotate depends on PyTorch and FAISS; installing those as prebuilt conda binaries avoids
compiler-toolchain issues. Create a ready-to-use environment from the provided file:

```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate synnotate
```

or into an existing environment:

```bash
conda install -c conda-forge -c bioconda python=3.11 numpy pandas pyarrow pytorch faiss-cpu biopython pyrodigal
pip install --no-deps synnotate
```

Synnotate targets **Python 3.10–3.12**.

Then download the data Synnotate needs into a folder of your choice (do this once):

```bash
synnotate setup --type prokaryote --dir ./synnotate_db
```

You point to that folder with `--bundle ./synnotate_db` whenever you annotate.

## Try it

The repository includes a tiny example: a cluster of ribosomal-protein genes with one gene hidden as
"hypothetical protein." Synnotate works out what the hidden gene is, from its neighbours alone:

```bash
synnotate annotate examples/demo.fna --type prokaryote --bundle ./synnotate_db \
    --gff examples/demo.gff --interpret --out demo
```

```
gene_id  annotation  prediction                   confidence  synteny_support  trusted
dcw_03   UNKNOWN     cell division protein FtsL    0.77        0.65             0.95
```

`confidence` is a **calibrated expected accuracy** (a 0.95 means ~95% of such calls are right), and
`trusted` flags the strictest accuracy tier — `0.99`, `0.95`, or blank — that the call meets from
confidence and synteny together.

## Annotate your genome

Pick the line that matches what you have:

```bash
# You already have an annotated genome (a GFF with product names — e.g. from NCBI, Prokka, or Bakta)
synnotate annotate genome.fna --type prokaryote --bundle ./synnotate_db --gff annotations.gff --interpret --out result

# You have only a genome sequence, not annotated yet
synnotate annotate genome.fna --type prokaryote --bundle ./synnotate_db --backend eggnog --out result

# A phage genome
synnotate annotate phage.fna --type phage --bundle ./synnotate_db --interpret --out result
```

Add `--interpret` to include the supporting evidence (synteny and neighbour breakdown). Leave it off
for a faster, prediction-only run.

### What Synnotate needs from you

Synnotate reads the functions of a gene's neighbours, so those neighbours have to be named first. How
that happens depends on your input:

- **An already-annotated genome** (a GFF with `product=` names) — Synnotate reads the names directly.
  Nothing else to install.
- **A genome sequence only** — Synnotate first names the genes with
  [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) (install it separately, below).
- **A phage** — Synnotate uses [Pharokka](https://github.com/gbouras13/pharokka).
- **Your own names** — give Synnotate a simple `gene_id<TAB>function` table with `--annotations`.

To use the sequence-only path, install eggNOG-mapper once:

```bash
mamba install -c bioconda eggnog-mapper
download_eggnog_data.py
```

## Your results

Synnotate writes `result.synnotate.tsv`, one row per gene:

- **prediction** — the predicted function
- **confidence** — a **calibrated expected accuracy** (isotonic regression), so 0.95 means ~95% of
  such predictions are correct
- **confidence_raw** — the uncalibrated softmax score, for reference
- **top5** — the five most likely functions
- **category** — a broad functional category

With `--interpret`, you also get:

- **synteny_support** — how strongly other genomes back up the call, from 0 to 1
- **trusted** — the strictest expected-accuracy tier the call clears from calibrated confidence ×
  synteny (`0.99`, `0.95`, or blank); the same trusted-region gate used in the paper
- **driving_neighbours** — which neighbouring genes led to the prediction

Add `--gff-out result.gff` to also get an annotated GFF you can load into a genome browser.

## Options at a glance

| Option | What it does |
|---|---|
| `--type {prokaryote,phage}` | which organism |
| `--gff FILE` | use the gene coordinates and names from this GFF |
| `--annotations FILE` | supply your own `gene_id<TAB>function` table |
| `--interpret` | include synteny support and the neighbour breakdown |
| `--out NAME` | name for the output files |
| `--gff-out FILE` | also write an annotated GFF |
| `--device {auto,cpu,cuda}` | run on CPU or GPU (auto by default) |

## Citation & license

Paper in preparation. Released under the MIT license — see [LICENSE](LICENSE).
