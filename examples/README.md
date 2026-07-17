# Example data

A tiny worked example that runs with **no external annotation tool** — it uses the
GFF-product-matching backend (`--backend gff`).

## Files

- **`demo.gff`** — the conserved **S10–spc ribosomal-protein operon**: 19 CDS in operon order, each
  with a `product=` string. The central gene, **`rps_09` (ribosomal protein L16)**, has been
  deliberately blanked to `product=hypothetical protein` — it is *dark* to homology.
- **`demo.fna`** — a matching contig (synthetic sequence; the `--gff` path reads coordinates, not
  bases).

## Run

```bash
synnotate setup    --type prokaryote --dir ./synnotate_db          # once
synnotate annotate examples/demo.fna --type prokaryote --bundle ./synnotate_db \
    --gff examples/demo.gff --backend gff --interpret --out demo
```

## Expected result

Synnotate recovers the blanked gene from its neighbourhood alone, with near-certain confidence and
strong synteny support:

```
gene_id  annotation  prediction                    confidence  adjusted_confidence  synteny_support
rps_09   UNKNOWN     50S ribosomal protein L16      1.00        0.99                 0.99
```

`driving_neighbours` shows why: `-1 S3`, `+1 L29`, `+2 S17` — the conserved ribosomal context. This
is the core capability in miniature: a gene with **no usable product name** is annotated from
**where it sits**, corroborated by synteny, with an explicit per-neighbour attribution.

The other 18 genes carry real product strings, which Synnotate matches onto its vocabulary to build
the context (`named->vocab 18`). This example is synthetic (one contig, placeholder sequence); on
real genomes give Synnotate a full assembly plus its GFF, or let it call genes itself.
