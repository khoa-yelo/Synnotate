# Synnotate

**Context-based annotation of bacterial and phage genomes"

Synnotate annotates genes remains unannotated after sequence and structural similarity search ("hypothetical/uncharacterized protein") based on neighbor genes annnotation

## Install

```bash
pip install git+https://github.com/khoa-yelo/Synnotate.git
```

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
gene_id  prediction                   confidence  synteny_support
rps_09   50S ribosomal protein L16    1.00        0.99
```

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
- **confidence** — how sure Synnotate is, from 0 to 1
- **top5** — the five most likely functions
- **category** — a broad functional category

With `--interpret`, you also get:

- **synteny_support** — how strongly other genomes back up the call, from 0 to 1
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
