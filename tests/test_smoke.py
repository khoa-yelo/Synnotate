"""Smoke tests: the package imports standalone (no training repo) and gene-calling + feature
binning work offline. The model-inference / interpretation paths need a bundle (tested separately)."""
from synnotate import genecall, features, __version__


def test_version():
    assert __version__


def test_genecall_and_fasta(tmp_path):
    fa = tmp_path / "g.fa"
    seq = "".join("ATG" + "GCA" * 200 + "TAA" for _ in range(3))
    fa.write_text(">c1\n" + seq + "\n")
    genes = genecall.call_genes(str(fa), phage=False)
    assert len(genes) >= 1
    assert genes[0].strand in "+-"
    assert genes[0].protein


def test_features_bins():
    assert features.distance_bin(120) == features.DIST_BIN_IDX["51_200"]
    assert features.distance_bin(None) == features.DIST_BIN_IDX["edge"]
    assert features.STRAND_IDX["+"] == 0


def test_standalone_imports():
    # importing the package must not require the training repo (cohyp)
    import sys
    import synnotate.model, synnotate.context, synnotate.interpret, synnotate.bundle, synnotate.pipeline  # noqa
    assert "cohyp" not in sys.modules
