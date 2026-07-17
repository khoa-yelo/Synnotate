"""Feature helpers (vendored, standalone): distance/length/GC bins, strand + segment
encoders, and context-window construction. Independent of the training repo (no cohyp import)."""
from __future__ import annotations

# 6 distance bins. Index 5 reserved for contig-edge / unknown.
DIST_BINS = ["overlap", "0_50", "51_200", "201_1000", "gt_1000", "edge"]
DIST_BIN_IDX = {b: i for i, b in enumerate(DIST_BINS)}


def distance_bin(gap: int | None) -> int:
    """Map an intergenic gap (bp) to a bin index. None -> edge/unknown."""
    if gap is None:
        return DIST_BIN_IDX["edge"]
    if gap < 0:
        return DIST_BIN_IDX["overlap"]
    if gap <= 50:
        return DIST_BIN_IDX["0_50"]
    if gap <= 200:
        return DIST_BIN_IDX["51_200"]
    if gap <= 1000:
        return DIST_BIN_IDX["201_1000"]
    return DIST_BIN_IDX["gt_1000"]


SEG_LEFT, SEG_TARGET, SEG_RIGHT = 0, 1, 2
STRAND_IDX = {"+": 0, "-": 1, ".": 2}
STRAND_MASK = STRAND_IDX["."]  # target's own strand never leaks (context only)

FEATURE_TYPES = ["PAD", "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                 "riboswitch", "pseudogene", "EDGE", "MASK"]
FEATURE_TYPE_IDX = {t: i for i, t in enumerate(FEATURE_TYPES)}


def feature_type_idx(ft: str | None) -> int:
    return FEATURE_TYPE_IDX.get(ft, FEATURE_TYPE_IDX["CDS"])


LEN_BINS = ["PAD", "EDGE", "MASK", "lt150", "150_300", "300_600",
            "600_900", "900_1500", "1500_3000", "gt3000"]
LEN_BIN_IDX = {b: i for i, b in enumerate(LEN_BINS)}


def length_bin(bp: int | None) -> int:
    if bp is None:
        return LEN_BIN_IDX["EDGE"]
    if bp < 150:
        return LEN_BIN_IDX["lt150"]
    if bp < 300:
        return LEN_BIN_IDX["150_300"]
    if bp < 600:
        return LEN_BIN_IDX["300_600"]
    if bp < 900:
        return LEN_BIN_IDX["600_900"]
    if bp < 1500:
        return LEN_BIN_IDX["900_1500"]
    if bp < 3000:
        return LEN_BIN_IDX["1500_3000"]
    return LEN_BIN_IDX["gt3000"]


GC_BINS = ["PAD", "EDGE", "MASK", "lt30", "30_40", "40_50",
           "50_60", "60_70", "gt70", "unknown"]
GC_BIN_IDX = {b: i for i, b in enumerate(GC_BINS)}


def gc_bin(gc: float | None) -> int:
    if gc is None or gc < 0:
        return GC_BIN_IDX["unknown"]
    if gc < 0.30:
        return GC_BIN_IDX["lt30"]
    if gc < 0.40:
        return GC_BIN_IDX["30_40"]
    if gc < 0.50:
        return GC_BIN_IDX["40_50"]
    if gc < 0.60:
        return GC_BIN_IDX["50_60"]
    if gc < 0.70:
        return GC_BIN_IDX["60_70"]
    return GC_BIN_IDX["gt70"]


FEATURE_MASK = FEATURE_TYPE_IDX["MASK"]; FEATURE_EDGE = FEATURE_TYPE_IDX["EDGE"]
LEN_MASK = LEN_BIN_IDX["MASK"]; LEN_EDGE = LEN_BIN_IDX["EDGE"]
GC_MASK = GC_BIN_IDX["MASK"]; GC_EDGE = GC_BIN_IDX["EDGE"]
