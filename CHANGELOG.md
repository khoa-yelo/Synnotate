# Changelog

## 0.2.0 (2026-08-02)

Calibrated, trusted-region release on the v3 dim-384 models.

- **Calibrated confidence.** `confidence` is now an isotonic-calibrated expected accuracy (a 0.95
  means ~95% of such calls are correct) instead of a raw softmax score. The uncalibrated value is
  still reported as `confidence_raw`. Calibration parameters ship in the bundle's `calibration.json`.
- **Trusted-region gate.** With `--interpret`, a `trusted` column reports the strictest
  expected-accuracy tier (`0.99`, `0.95`, or blank) the call meets from calibrated confidence ×
  synteny support — the same gate used in the paper. Replaces the previous geometric-mean
  `adjusted_confidence`.
- **Bundles rebuilt at dim-384** (v3) for both organisms, including the first **phage** bundle
  (no-unknown-function vocabulary, 582 labels). Each bundle is self-contained with a MANIFEST.json
  (per-file sha256) and passes `synnotate setup --check`.
- Synteny NW scorer aligned to the deployed scorer (mismatch −1, mask/edge 0) so the synteny score
  matches the trusted-region contours.
- Added `LICENSE` (MIT) and this changelog.

## 0.1.0

Initial standalone package: gene-calling, curated-context annotation, additive context transformer,
kNN+MSA synteny, and per-neighbour attribution.
