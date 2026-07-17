# Example data

A tiny worked example that runs with **no external annotation tool** — it uses the
GFF-product-matching backend (`--backend gff`).

## Files

- **`demo.gff`** — a bacterial **division/cell-wall (`dcw`) gene cluster**: 13 CDS in operon order
  (MraZ, RsmH, FtsL, FtsI, MurF, MraY, FtsW, MurG, MurC, Ddl, FtsQ, FtsA, FtsZ), each with a
  `product=` string. The third gene, **`ftsL`**, has been deliberately blanked to
  `product=hypothetical protein` — it is *dark* to homology.
- **`demo.fna`** — a matching contig (synthetic sequence; the `--gff` path reads coordinates, not
  bases).

## Run

```bash
synnotate setup    --type prokaryote --dir ./synnotate_db          # once
synnotate annotate examples/demo.fna --type prokaryote --bundle ./synnotate_db \
    --gff examples/demo.gff --backend gff --interpret --out demo
```

## Expected result

Synnotate recovers the blanked gene from its neighbourhood alone:

```
gene_id  annotation  prediction                   confidence  synteny_support
dcw_03   UNKNOWN     cell division protein FtsL    0.55        0.63
```

`driving_neighbours` shows why: `-1 RsmH`, `+1 FtsI`, `+2 MurF` — the conserved `dcw` context. This
is the core capability in miniature: a gene with **no usable product name** is annotated from
**where it sits**, corroborated by synteny, with an explicit per-neighbour attribution.

The other 12 genes carry real product strings, which Synnotate matches onto its vocabulary to build
the context (`named->vocab 12`). This example is synthetic (one contig, placeholder sequence); on
real genomes give Synnotate a full assembly plus its GFF, or let it call genes itself.
