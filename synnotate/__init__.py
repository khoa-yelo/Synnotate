"""Synnotate — synteny-transformer gene-function annotation from genomic context.

Standalone package (does not import the training repo). Given a genome FASTA (and optional
GFF), it calls genes, builds a curated-vocabulary neighbourhood context, runs the interpretable
additive context transformer, and reports per-gene function predictions with softmax confidence,
a k=5 kNN+MSA synteny-adjusted confidence, and exact per-neighbour attributions."""
__version__ = "0.1.0"
from synnotate.model import ContextTransformerAdditive, sparsemax  # noqa: F401
