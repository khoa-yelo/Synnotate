"""Gene calling + FASTA/GFF parsing (standalone). Uses pyrodigal (prokaryote) / pyrodigal-gv
(phage) when no GFF is supplied. Returns a per-contig-ordered list of genes with coordinates,
strand, and translated protein sequence — the input Synnotate needs to build genomic context."""
from __future__ import annotations
import gzip, os
from dataclasses import dataclass, field


@dataclass
class Gene:
    gene_id: str
    contig: str
    start: int          # 1-based inclusive
    end: int
    strand: str         # '+' / '-'
    protein: str = ""   # amino-acid sequence (may be empty if from a GFF without translation)
    product: str = ""   # annotation family/product (filled by the annotation step)
    extra: dict = field(default_factory=dict)


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    """Yield (name, sequence) records from a (optionally gzipped) FASTA."""
    name, seq = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]; seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def call_genes(fasta, phage=False):
    """Call genes with pyrodigal(-gv). Returns list[Gene] in (contig, start) order with proteins."""
    try:
        if phage:
            import pyrodigal_gv as _pg
            orf = _pg.ViralGeneFinder(meta=True)
        else:
            import pyrodigal as _pg
            orf = _pg.GeneFinder(meta=True)
    except ImportError as e:  # pragma: no cover
        raise SystemExit(f"gene-calling needs pyrodigal{'-gv' if phage else ''}: {e}")
    genes = []
    for contig, seq in read_fasta(fasta):
        for i, g in enumerate(orf.find_genes(seq.encode()), 1):
            genes.append(Gene(gene_id=f"{contig}_{i:05d}", contig=contig,
                              start=g.begin, end=g.end, strand="+" if g.strand == 1 else "-",
                              protein=g.translate().rstrip("*")))
    genes.sort(key=lambda x: (x.contig, x.start, x.end))
    return genes


def read_gff(gff, fasta=None):
    """Parse CDS features from a GFF3 (coordinates + strand + ID). Proteins are left empty unless a
    matching FASTA of translations is supplied separately by the caller."""
    genes = []
    with _open(gff) as fh:
        for line in fh:
            if not line.strip() or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "gene"):
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            gid = attrs.get("ID") or attrs.get("locus_tag") or f"{f[0]}_{f[3]}"
            genes.append(Gene(gene_id=gid, contig=f[0], start=int(f[3]), end=int(f[4]),
                              strand=f[6] if f[6] in "+-" else "+",
                              product=attrs.get("product", "")))
    genes.sort(key=lambda x: (x.contig, x.start, x.end))
    return genes


def write_proteins(genes, path):
    """Write called proteins to a FASTA (for eggNOG-mapper / Pharokka input)."""
    with open(path, "w") as fh:
        for g in genes:
            if g.protein:
                fh.write(f">{g.gene_id}\n{g.protein}\n")
    return path
