"""`synnotate` command-line interface.

  synnotate annotate genome.fasta --type {phage|prokaryote} [--gff G] \
                     [--annotations A.tsv] [--interpret] [--k 5] [--out prefix]
  synnotate setup   --type {phage|prokaryote} [--dir DIR]     # fetch the interpretation bundle
  synnotate --version
"""
from __future__ import annotations
import argparse, os, sys, csv, json


def _default_bundle_dir(organism):
    root = os.environ.get("SYNNOTATE_HOME", os.path.expanduser("~/.synnotate"))
    return os.path.join(root, "bundle", organism)


def cmd_annotate(a):
    import torch
    from synnotate.bundle import Bundle
    from synnotate import pipeline
    organism = "phage" if a.type == "phage" else "prokaryote"
    bdir = a.bundle or _default_bundle_dir(organism)
    if not os.path.isdir(bdir):
        sys.exit(f"bundle not found at {bdir} — run `synnotate setup --type {a.type}` first")
    device = "cuda" if (a.device == "auto" and torch.cuda.is_available()) or a.device == "cuda" else "cpu"
    bundle = Bundle(bdir)
    genes, rows = pipeline.run(a.fasta, bundle, gff=a.gff, annotations=a.annotations,
                               interpret=a.interpret, k=a.k, device=device,
                               workdir=a.out + "_work", threads=a.threads,
                               phage=(organism == "phage"), backend=a.backend)
    tsv = a.out + ".synnotate.tsv"
    cols = ["gene_id", "contig", "start", "end", "strand", "annotation", "prediction", "confidence",
            "category", "top5"]
    if a.interpret:
        cols += ["adjusted_confidence", "synteny_support", "mean_flank_ident",
                 "target_slot_conserved", "driving_neighbours"]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", delimiter="\t")
        w.writeheader(); w.writerows(rows)
    n_dark = sum(1 for r in rows if r["annotation"] in ("UNKNOWN", "hypothetical protein", ""))
    print(f"[synnotate] {len(rows)} genes | {n_dark} hypothetical | -> {tsv}")
    if a.gff_out:
        _write_gff(rows, a.gff_out); print(f"[synnotate] annotated GFF -> {a.gff_out}")


def _write_gff(rows, path):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for r in rows:
            attrs = f"ID={r['gene_id']};product={r['prediction']};synnotate_confidence={r['confidence']}"
            if "adjusted_confidence" in r:
                attrs += f";synnotate_adjusted={r['adjusted_confidence']}"
            fh.write(f"{r['contig']}\tsynnotate\tCDS\t{r['start']}\t{r['end']}\t.\t{r['strand']}\t.\t{attrs}\n")


# Bundle registry — the download URL + md5 of each published bundle tarball. Fill `url`/`md5` at
# release time (e.g. a Zenodo record) or override per-run with $SYNNOTATE_BUNDLE_URL_<ORGANISM>.
BUNDLE_REGISTRY = {
    "prokaryote": {"version": "v3",
                   "url": os.environ.get("SYNNOTATE_BUNDLE_URL_PROKARYOTE", ""),
                   "md5": os.environ.get("SYNNOTATE_BUNDLE_MD5_PROKARYOTE", "")},
    "phage":      {"version": "v3",
                   "url": os.environ.get("SYNNOTATE_BUNDLE_URL_PHAGE", ""),
                   "md5": os.environ.get("SYNNOTATE_BUNDLE_MD5_PHAGE", "")},
}


def _md5(path, buf=1 << 20):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(buf), b""):
            h.update(c)
    return h.hexdigest()


def _sha256(path, buf=1 << 20):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(buf), b""):
            h.update(c)
    return h.hexdigest()


# a handful of required files -> a CHEAP presence check for the idempotent skip: don't re-hash a
# multi-GB bundle on every run; full sha256 verification is reserved for `setup --check`.
_REQUIRED = ["config.json", "curated_context.pt", "vocab.json", "label_vocab.json",
             "index/ctx_cur.faiss", "index/meta.parquet"]


def _installed(bdir):
    return all(os.path.exists(os.path.join(bdir, f)) for f in _REQUIRED)


def _download(url, dest_path, ua):
    """Stream a URL to disk with a progress line and a User-Agent header (Zenodo prefers one)."""
    import urllib.request, sys as _sys
    req = urllib.request.Request(url, headers={"User-Agent": ua})
    with urllib.request.urlopen(req) as resp, open(dest_path, "wb") as out:
        total = int(resp.headers.get("Content-Length", 0)); done = 0
        while True:
            buf = resp.read(1 << 20)
            if not buf:
                break
            out.write(buf); done += len(buf)
            if total:
                _sys.stderr.write(f"\r[synnotate] downloading {done/1e6:7.0f}/{total/1e6:.0f} MB "
                                  f"({100*done/total:5.1f}%)")
            else:
                _sys.stderr.write(f"\r[synnotate] downloading {done/1e6:7.0f} MB")
            _sys.stderr.flush()
        _sys.stderr.write("\n")


def _verify_manifest(bdir):
    """Check every file listed in the bundle MANIFEST.json against its recorded sha256."""
    mp = os.path.join(bdir, "MANIFEST.json")
    if not os.path.exists(mp):
        print("[synnotate] no MANIFEST.json — skipping integrity check"); return True
    man = json.load(open(mp))
    bad = []
    for name, ent in man.items():
        want = ent.get("sha256") if isinstance(ent, dict) else None
        if not want:
            continue
        f = os.path.join(bdir, name)
        if not os.path.exists(f) or _sha256(f) != want:
            bad.append(name)
    if bad:
        print(f"[synnotate] MANIFEST verify FAILED for {len(bad)} file(s): {bad[:5]}"); return False
    print(f"[synnotate] MANIFEST verify OK ({len(man)} files)"); return True


def _extract_into(tar_path, dest):
    """Extract a bundle tarball into `dest`, flattening a single top-level wrapper dir if present."""
    import tarfile, tempfile, shutil
    os.makedirs(dest, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        with tarfile.open(tar_path) as tf:
            try:
                tf.extractall(td, filter="data")     # py>=3.12: reject unsafe members
            except TypeError:
                tf.extractall(td)                    # py<3.12: no filter kwarg
        entries = os.listdir(td)
        src = os.path.join(td, entries[0]) if (len(entries) == 1 and
              os.path.isdir(os.path.join(td, entries[0]))) else td
        for name in os.listdir(src):
            s, d = os.path.join(src, name), os.path.join(dest, name)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)


def cmd_setup(a):
    organism = "phage" if a.type == "phage" else "prokaryote"
    dest = a.dir or _default_bundle_dir(organism)
    reg = BUNDLE_REGISTRY.get(organism, {})

    # --check: just validate an already-installed bundle
    if a.check:
        if not os.path.isdir(dest):
            sys.exit(f"[synnotate] no bundle at {dest}")
        sys.exit(0 if _verify_manifest(dest) else 1)

    # already installed (cheap presence check), unless --force. Full sha256 is `setup --check`.
    if _installed(dest) and not a.force:
        print(f"[synnotate] {organism} bundle already at {dest} (use --force to reinstall, "
              f"--check to verify sha256)")
        return

    os.makedirs(dest, exist_ok=True)

    # (1) install from a local tarball or directory
    if a.from_path:
        src = a.from_path
        if os.path.isdir(src):
            import shutil
            for name in os.listdir(src):
                s, d = os.path.join(src, name), os.path.join(dest, name)
                shutil.copytree(s, d, dirs_exist_ok=True) if os.path.isdir(s) else shutil.copy2(s, d)
        else:
            _extract_into(src, dest)
        ok = _installed(dest)
        print(f"[synnotate] installed {organism} bundle from {src} -> {dest}"
              f"{'' if ok else '  [WARN: required files missing]'}")
        return

    # (2) download from the registry / an explicit --url
    from synnotate import __version__
    url = a.url or reg.get("url", "")
    md5 = reg.get("md5", "")
    if not url:
        sys.exit(
            f"[synnotate] no download URL for the {organism} bundle yet.\n"
            f"  Provide one with `--url <tarball>` or $SYNNOTATE_BUNDLE_URL_{organism.upper()},\n"
            f"  or install a local copy with `synnotate setup --type {a.type} --from <dir-or-tarball>`.")
    tar = os.path.join(dest, f"_{organism}_bundle.tar.gz")
    print(f"[synnotate] downloading {organism} bundle ({reg.get('version','?')}) from {url}")
    _download(url, tar, ua=f"synnotate/{__version__}")
    if md5:
        got = _md5(tar)
        if got != md5:
            os.remove(tar); sys.exit(f"[synnotate] md5 mismatch: got {got}, expected {md5}")
        print(f"[synnotate] md5 OK ({md5})")
    else:
        print("[synnotate] WARNING: no md5 in registry — skipping tarball integrity check")
    _extract_into(tar, dest)
    os.remove(tar)
    ok = _installed(dest)
    print(f"[synnotate] installed {organism} bundle -> {dest}"
          f"{'' if ok else '  [WARN: required files missing]'}")


def main(argv=None):
    from synnotate import __version__
    ap = argparse.ArgumentParser(prog="synnotate", description="Synteny-transformer gene annotation")
    ap.add_argument("--version", action="version", version=f"synnotate {__version__}")
    sub = ap.add_subparsers(dest="cmd")

    an = sub.add_parser("annotate", help="annotate a genome FASTA")
    an.add_argument("fasta")
    an.add_argument("--type", required=True, choices=["phage", "prokaryote"])
    an.add_argument("--backend", choices=["auto", "eggnog", "gff"], default="auto",
                    help="prokaryote neighbour-annotation source: 'gff' string-matches the GFF "
                         "product= column, 'eggnog' runs eggNOG-mapper, 'auto' (default) uses the "
                         "GFF products when the input GFF is already annotated else eggNOG-mapper")
    an.add_argument("--gff", help="use CDS coordinates (and product= strings) from this GFF instead of gene-calling")
    an.add_argument("--annotations", help="precomputed gene_id<TAB>family TSV (skip the backend)")
    an.add_argument("--interpret", action="store_true", help="add kNN+MSA synteny + attributions")
    an.add_argument("--k", type=int, default=5)
    an.add_argument("--out", default="synnotate_out")
    an.add_argument("--gff-out", help="also write an annotated GFF here")
    an.add_argument("--bundle", help="bundle directory (default ~/.synnotate/bundle/<organism>)")
    an.add_argument("--threads", type=int, default=16)
    an.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    an.set_defaults(func=cmd_annotate)

    st = sub.add_parser("setup", help="download/install the interpretation bundle (models+index)")
    st.add_argument("--type", required=True, choices=["phage", "prokaryote"])
    st.add_argument("--dir", help="install location (default ~/.synnotate/bundle/<organism>)")
    st.add_argument("--from", dest="from_path", help="install from a local bundle dir or .tar.gz instead of downloading")
    st.add_argument("--url", help="download the bundle tarball from this URL")
    st.add_argument("--force", action="store_true", help="reinstall even if a bundle is already present")
    st.add_argument("--check", action="store_true", help="only verify an installed bundle against its MANIFEST")
    st.set_defaults(func=cmd_setup)

    a = ap.parse_args(argv)
    if not getattr(a, "cmd", None):
        ap.print_help(); return
    a.func(a)


if __name__ == "__main__":
    main()
