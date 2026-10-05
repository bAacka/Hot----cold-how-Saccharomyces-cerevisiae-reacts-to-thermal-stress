#!/usr/bin/env python3

from pathlib import Path
import hashlib
import re

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Circle

from scipy.ndimage import gaussian_filter
from scipy.stats import spearmanr, pearsonr


                                                              
       
                                                              

ROOT = Path("/")

BRIDGE = ROOT / "44_phenotype_bridge"

FORMAL = (
    BRIDGE
    / "focus_metric_validation"
    / "190321_WT"
    / "formal_manual_comparison"
)

HELDOUT = (
    BRIDGE
    / "focus_metric_validation"
    / "190321_WT"
    / "frozen_v3"
)

PHENO = (
    BRIDGE
    / "46_phenotype_state"
    / "46n_Guk1_recovery_phenotype"
)

MOVIE = (
    BRIDGE
    / "45_published_validation"
    / "published_movie_frames_v2"
)

TIME_ALIGNED = (
    FORMAL
    / "TIME_ALIGNED_COMPARISON.tsv"
)

CONCORDANCE = (
    FORMAL
    / "TRAJECTORY_CONCORDANCE.tsv"
)

ENDPOINT = (
    FORMAL
    / "ENDPOINT_DISTRIBUTION_SUMMARY.tsv"
)

CELL_LONG = (
    HELDOUT
    / "HELDOUT_FOCUS_LONG_COMPLETE_T1_T12.tsv"
)

PHENO_CELL = (
    PHENO
    / "46n_frozen_cell_Guk1_recovery_phenotype.tsv"
)

PROJECTION_CACHE = (
    BRIDGE
    / "projection_cache_190321_WT"
)

GFP_STACK = (
    PROJECTION_CACHE
    / "egfp_max_projection.npy"
)

DIC_STACK = (
    PROJECTION_CACHE
    / "dic_best_gfp_z.npy"
)

OUTROOT = ROOT / "54_figures"

CURATED = (
    OUTROOT
    / "02_curated_inputs"
    / "Fig4_5"
)

PANELS = (
    OUTROOT
    / "04_panels"
    / "Fig4_5_microscopy"
)

FINAL = OUTROOT / "05_final"
LEGENDS = OUTROOT / "06_legends"

for d in [
    CURATED,
    PANELS,
    FINAL,
    LEGENDS,
]:
    d.mkdir(
        parents=True,
        exist_ok=True,
    )


                                                              
                         
                                                              

EXPECTED_N = 301
EXPECTED_LOWER = 296
EXPECTED_LOWER_FRAC = 0.9834

EXPECTED_AUTO_ENDPOINT = 0.491
EXPECTED_MANUAL_ENDPOINT = 0.579

EXPECTED_RAW_SPEARMAN = 0.972
EXPECTED_RAW_PEARSON = 0.726

EXPECTED_NORM_SPEARMAN = 0.986
EXPECTED_NORM_PEARSON = 0.799


                                                              
         
                                                              

def fail(msg):
    raise RuntimeError(msg)


def sha256(path):

    h = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def norm_name(x):

    return re.sub(
        r"[^a-z0-9]+",
        "_",
        str(x).lower(),
    ).strip("_")


def numeric(df, col):

    return pd.to_numeric(
        df[col],
        errors="coerce",
    )


def choose_col(
    df,
    required=(),
    forbidden=(),
):

    hits = []

    for c in df.columns:

        n = norm_name(c)

        if not all(
            x in n
            for x in required
        ):
            continue

        if any(
            x in n
            for x in forbidden
        ):
            continue

        x = numeric(
            df,
            c,
        )

        if x.notna().sum() == 0:
            continue

        hits.append(c)

    if len(hits) == 1:
        return hits[0]

    return None


def parse_time_value(x):

    if pd.isna(x):
        return np.nan

    if isinstance(
        x,
        (int, float, np.integer, np.floating),
    ):
        return float(x)

    s = str(x)

    m = re.search(
        r"(\d+(?:\.\d+)?)",
        s,
    )

    if not m:
        return np.nan

    return float(
        m.group(1)
    )


def image_gray(img):

    arr = np.asarray(img)

    if arr.ndim == 2:
        return arr.astype(float)

    if arr.shape[-1] >= 3:

        rgb = arr[..., :3].astype(float)

        return (
            0.2126 * rgb[..., 0]
            + 0.7152 * rgb[..., 1]
            + 0.0722 * rgb[..., 2]
        )

    return arr[..., 0].astype(float)


                                                              
                
                                                              

required_files = [
    TIME_ALIGNED,
    CONCORDANCE,
    ENDPOINT,
    CELL_LONG,
    PHENO_CELL,
    GFP_STACK,
    DIC_STACK,
]

missing = [
    str(p)
    for p in required_files
    if not p.exists()
]

if missing:
    fail(
        "Missing required Figure 4.5 inputs:\n"
        + "\n".join(missing)
    )


                                                              
                                                   
                                                              

cmp = pd.read_csv(
    TIME_ALIGNED,
    sep="\t",
    low_memory=False,
)

print(
    "TIME_ALIGNED columns:",
    list(cmp.columns),
)


                              
             
                              

time_col = None

for candidate in [
    "time_min",
    "minutes",
    "minute",
    "time",
    "timepoint",
]:

    for c in cmp.columns:

        if norm_name(c) == candidate:
            time_col = c
            break

    if time_col:
        break

if time_col is None:

    for c in cmp.columns:

        n = norm_name(c)

        if (
            "time" in n
            or "minute" in n
        ):
            time_col = c
            break

if time_col is None:
    fail(
        "Could not identify time column in "
        "TIME_ALIGNED_COMPARISON.tsv"
    )

cmp["_time"] = (
    cmp[time_col]
    .map(parse_time_value)
)


                              
                                                    
                              

                                                
                                                          
                                                            
                                 

manual_raw = "manual_median_ratio"
auto_raw = "automated_median_ratio"

manual_norm = "manual_median_relative_T1"
auto_norm = "automated_median_relative_T1"

required_trajectory_cols = [
    manual_raw,
    auto_raw,
    manual_norm,
    auto_norm,
]

missing = [
    c for c in required_trajectory_cols
    if c not in cmp.columns
]

if missing:
    fail(
        "TIME_ALIGNED_COMPARISON.tsv missing required columns: "
        + ", ".join(missing)
    )


plot_cmp = (
    cmp[
        [
            "_time",
            manual_raw,
            auto_raw,
        ]
    ]
    .copy()
)

plot_cmp[manual_raw] = pd.to_numeric(
    plot_cmp[manual_raw],
    errors="coerce",
)

plot_cmp[auto_raw] = pd.to_numeric(
    plot_cmp[auto_raw],
    errors="coerce",
)

plot_cmp = (
    plot_cmp
    .dropna()
    .sort_values("_time")
)

plot_cmp.columns = [
    "time",
    "manual_raw",
    "automated_raw",
]


norm_cmp = (
    cmp[
        [
            "_time",
            manual_norm,
            auto_norm,
        ]
    ]
    .copy()
)

norm_cmp[manual_norm] = pd.to_numeric(
    norm_cmp[manual_norm],
    errors="coerce",
)

norm_cmp[auto_norm] = pd.to_numeric(
    norm_cmp[auto_norm],
    errors="coerce",
)

norm_cmp = (
    norm_cmp
    .dropna()
    .sort_values("_time")
)

norm_cmp.columns = [
    "time",
    "manual_norm",
    "automated_norm",
]


                              
                       
                              

raw_rho = spearmanr(
    plot_cmp["manual_raw"],
    plot_cmp["automated_raw"],
).statistic

raw_r = pearsonr(
    plot_cmp["manual_raw"],
    plot_cmp["automated_raw"],
).statistic

norm_rho = spearmanr(
    norm_cmp["manual_norm"],
    norm_cmp["automated_norm"],
).statistic

norm_r = pearsonr(
    norm_cmp["manual_norm"],
    norm_cmp["automated_norm"],
).statistic


for label, obs, exp in [
    (
        "raw Spearman",
        raw_rho,
        EXPECTED_RAW_SPEARMAN,
    ),
    (
        "raw Pearson",
        raw_r,
        EXPECTED_RAW_PEARSON,
    ),
    (
        "normalized Spearman",
        norm_rho,
        EXPECTED_NORM_SPEARMAN,
    ),
    (
        "normalized Pearson",
        norm_r,
        EXPECTED_NORM_PEARSON,
    ),
]:

    if abs(obs - exp) > 0.035:
        fail(
            f"{label} mismatch: "
            f"observed={obs:.4f}, expected≈{exp:.4f}"
        )


plot_cmp.to_csv(
    CURATED
    / "Fig4_5_manual_automated_raw.tsv",
    sep="\t",
    index=False,
)

norm_cmp.to_csv(
    CURATED
    / "Fig4_5_manual_automated_normalized.tsv",
    sep="\t",
    index=False,
)


                                                              
                                  
                                                              

cells = pd.read_csv(
    CELL_LONG,
    sep="\t",
    low_memory=False,
)

print(
    "CELL_LONG columns:",
    list(cells.columns),
)


                              
               
                              

id_cols = []

for wanted in [
    "acquisition",
    "acquisition_id",
    "scene",
    "scene_id",
    "track_id",
    "cell_id",
    "cell",
]:

    for c in cells.columns:

        if norm_name(c) == wanted:
            if c not in id_cols:
                id_cols.append(c)

                                                                   
strong_ids = [
    c for c in id_cols
    if (
        "track" in norm_name(c)
        or "cell_id" in norm_name(c)
    )
]

if strong_ids:

    context_ids = [
        c for c in id_cols
        if (
            "scene" in norm_name(c)
            or "acquisition" in norm_name(c)
        )
    ]

    id_cols = (
        context_ids
        + strong_ids[:1]
    )

if not id_cols:
    fail(
        "Could not identify cell identity columns in "
        "HELDOUT_FOCUS_LONG_COMPLETE_T1_T12.tsv"
    )

cells["_cell"] = (
    cells[id_cols]
    .astype(str)
    .agg(
        "|".join,
        axis=1,
    )
)


                              
      
                              

cell_time_col = None

for c in cells.columns:

    n = norm_name(c)

    if n in {
        "time_min",
        "minutes",
        "minute",
        "time",
        "timepoint",
        "t",
    }:
        cell_time_col = c
        break

if cell_time_col is None:

    for c in cells.columns:

        if "time" in norm_name(c):
            cell_time_col = c
            break

if cell_time_col is None:
    fail(
        "Could not identify time column in cell trajectory table."
    )

cells["_time"] = (
    cells[cell_time_col]
    .map(parse_time_value)
)


                              
                     
                              

                                         
                                                            
                                                              
ratio_col = "focus_to_cytosol_ratio"

if ratio_col not in cells.columns:
    fail(
        f"Missing required trajectory column: {ratio_col}"
    )

cells["_ratio"] = numeric(
    cells,
    ratio_col,
)

cells = (
    cells[
        [
            "_cell",
            "_time",
            "_ratio",
        ]
    ]
    .dropna()
    .sort_values(
        [
            "_cell",
            "_time",
        ]
    )
)


                              
                               
                              

counts = (
    cells.groupby("_cell")
    .size()
)

n_cells = int(
    len(counts)
)

if n_cells != EXPECTED_N:
    fail(
        f"Expected {EXPECTED_N} complete cells; "
        f"found {n_cells}"
    )


first = (
    cells.groupby("_cell")
    .first()["_ratio"]
)

last = (
    cells.groupby("_cell")
    .last()["_ratio"]
)

lower = int(
    (last < first).sum()
)

if lower != EXPECTED_LOWER:
    fail(
        f"Expected {EXPECTED_LOWER}/{EXPECTED_N} "
        f"cells lower at endpoint; found {lower}/{n_cells}"
    )


                              
                   
                              

baseline = (
    cells.groupby("_cell")["_ratio"]
    .transform("first")
)

cells["_Y"] = np.log2(
    cells["_ratio"]
    / baseline
)

cells.to_csv(
    CURATED
    / "Fig4_5_301_cell_normalized_trajectories.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


                                                              
                                 
                                                              

ph = pd.read_csv(
    PHENO_CELL,
    sep="\t",
    low_memory=False,
)

print(
    "PHENO columns:",
    list(ph.columns),
)

auc_candidates = []

for c in ph.columns:

    n = norm_name(c)

    if "auc" in n:

        x = numeric(
            ph,
            c,
        )

        if x.notna().sum() > 0:

            score = 0

            if "signed" in n:
                score += 3

            if (
                "t1" in n
                and "t8" in n
            ):
                score += 2

            if (
                "guk1" in n
                or "recovery" in n
                or "phenotype" in n
            ):
                score += 1

            auc_candidates.append(
                (
                    score,
                    c,
                )
            )

if not auc_candidates:
    fail(
        "Could not identify signed-AUC phenotype column."
    )

auc_candidates.sort(
    reverse=True
)

auc_col = auc_candidates[0][1]

ph["_Pi"] = numeric(
    ph,
    auc_col,
)

phenotype = (
    ph["_Pi"]
    .dropna()
    .to_numpy(float)
)

if len(phenotype) == 0:
    fail(
        "Frozen phenotype column contains no numeric values."
    )

pd.DataFrame(
    {
        "P_i": phenotype,
    }
).to_csv(
    CURATED
    / "Fig4_5_frozen_signed_AUC.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                                
                                                              

cells_A = pd.read_csv(
    CELL_LONG,
    sep="\t",
    low_memory=False,
)

required_A_cols = [
    "scene",
    "track_id",
    "time_index",
    "biological_min",
    "centroid_y",
    "centroid_x",
    "focus_y",
    "focus_x",
    "cytosol_y",
    "cytosol_x",
    "focus_to_cytosol_ratio",
    "initial_focus_prominence_z",
]

missing_A = [
    c for c in required_A_cols
    if c not in cells_A.columns
]

if missing_A:
    fail(
        "Panel A held-out table missing columns: "
        + ", ".join(missing_A)
    )

for c in [
    "scene",
    "track_id",
    "time_index",
    "biological_min",
    "centroid_y",
    "centroid_x",
    "focus_y",
    "focus_x",
    "cytosol_y",
    "cytosol_x",
    "focus_to_cytosol_ratio",
    "initial_focus_prominence_z",
]:
    cells_A[c] = pd.to_numeric(
        cells_A[c],
        errors="coerce",
    )

cells_A = cells_A.dropna(
    subset=required_A_cols
)

cells_A["cell_key"] = (
    cells_A["scene"].astype(int).astype(str)
    + "|"
    + cells_A["track_id"].astype(int).astype(str)
)


                                                              
                                
                                                   
                                                              

gfp_raw = np.load(
    GFP_STACK,
    mmap_mode="r",
)

dic_raw = np.load(
    DIC_STACK,
    mmap_mode="r",
)


def canonicalize_projection(arr, name, n_scenes, max_t):

    arr = np.asarray(arr)

                                            
    arr = np.squeeze(arr)

    if arr.ndim != 4:
        fail(
            f"{name}: expected 4D projection stack after squeeze; "
            f"got shape {arr.shape}"
        )

    shape = arr.shape

                                                                  
    spatial_axes = sorted(
        range(4),
        key=lambda i: shape[i],
        reverse=True,
    )[:2]

    nonspatial = [
        i for i in range(4)
        if i not in spatial_axes
    ]

                                          
    spatial_axes = sorted(spatial_axes)

    a, b = nonspatial

    candidates = []

                     
    if (
        shape[a] >= n_scenes
        and shape[b] > max_t
    ):
        score = (
            abs(shape[a] - n_scenes)
            + abs(shape[b] - (max_t + 1))
        )

        candidates.append(
            (score, a, b)
        )

                     
    if (
        shape[b] >= n_scenes
        and shape[a] > max_t
    ):
        score = (
            abs(shape[b] - n_scenes)
            + abs(shape[a] - (max_t + 1))
        )

        candidates.append(
            (score, b, a)
        )

    if not candidates:
        fail(
            f"{name}: cannot infer scene/time axes from shape "
            f"{shape}; n_scenes={n_scenes}, max_t={max_t}"
        )

    candidates.sort()
    _, scene_axis, time_axis = candidates[0]

    order = [
        scene_axis,
        time_axis,
        *spatial_axes,
    ]

    out = np.transpose(
        arr,
        axes=order,
    )

    print(
        f"{name}: raw shape={shape} -> "
        f"canonical scene,time,y,x={out.shape}"
    )

    return out


n_scenes_A = int(
    cells_A["scene"].max()
) + 1

max_t_A = int(
    cells_A["time_index"].max()
)

gfp_stack = canonicalize_projection(
    gfp_raw,
    "GFP",
    n_scenes_A,
    max_t_A,
)

dic_stack = canonicalize_projection(
    dic_raw,
    "DIC",
    n_scenes_A,
    max_t_A,
)

if gfp_stack.shape != dic_stack.shape:
    fail(
        "GFP/DIC canonical stack shapes differ: "
        f"{gfp_stack.shape} vs {dic_stack.shape}"
    )

_, _, IMG_H, IMG_W = gfp_stack.shape


                                                              
                                                     
 
               
                                    
                             
                       
                                             
 
                                                              
                          
                                                              

DISPLAY_T = [1, 6, 12]

summary_rows = []

for key, g in cells_A.groupby(
    "cell_key",
    sort=True,
):

    gg = (
        g[g["time_index"].isin(DISPLAY_T)]
        .sort_values("time_index")
    )

    if list(
        gg["time_index"].astype(int)
    ) != DISPLAY_T:
        continue

    first = gg.iloc[0]
    last = gg.iloc[-1]

    initial_ratio = float(
        first["focus_to_cytosol_ratio"]
    )

    endpoint_relative = float(
        last["focus_to_cytosol_ratio"]
        / initial_ratio
    )

    z0 = float(
        first["initial_focus_prominence_z"]
    )

                                                      
    half = 32

    boundary_ok = True

    for _, row in gg.iterrows():

        cx0 = float(row["centroid_x"])
        cy0 = float(row["centroid_y"])

        if not (
            half <= cx0 < IMG_W - half
            and half <= cy0 < IMG_H - half
        ):
            boundary_ok = False
            break

    if not boundary_ok:
        continue

    summary_rows.append(
        {
            "cell_key": key,
            "scene": int(first["scene"]),
            "track_id": int(first["track_id"]),
            "initial_ratio": initial_ratio,
            "endpoint_relative": endpoint_relative,
            "initial_z": z0,
        }
    )

rep = pd.DataFrame(
    summary_rows
)

if len(rep) == 0:
    fail(
        "No eligible representative cell for Panel A."
    )

target_endpoint = float(
    rep["endpoint_relative"].median()
)

target_z = float(
    rep["initial_z"].quantile(0.75)
)

endpoint_scale = float(
    rep["endpoint_relative"].quantile(0.75)
    - rep["endpoint_relative"].quantile(0.25)
)

z_scale = float(
    rep["initial_z"].quantile(0.75)
    - rep["initial_z"].quantile(0.25)
)

endpoint_scale = max(
    endpoint_scale,
    1e-6,
)

z_scale = max(
    z_scale,
    1e-6,
)

rep["_score"] = (
    abs(
        rep["endpoint_relative"]
        - target_endpoint
    )
    / endpoint_scale
    +
    0.35
    * abs(
        rep["initial_z"]
        - target_z
    )
    / z_scale
)

rep = rep.sort_values(
    [
        "_score",
        "scene",
        "track_id",
    ]
)

chosen = rep.iloc[0]

REP_SCENE = int(
    chosen["scene"]
)

REP_TRACK = int(
    chosen["track_id"]
)

rep_cell = (
    cells_A[
        (cells_A["scene"].astype(int) == REP_SCENE)
        &
        (cells_A["track_id"].astype(int) == REP_TRACK)
        &
        (cells_A["time_index"].astype(int).isin(DISPLAY_T))
    ]
    .copy()
    .sort_values("time_index")
)

if len(rep_cell) != 3:
    fail(
        "Representative cell does not have all three display timepoints."
    )


                                                              
                                                                
                                                                
                                                              

CROP_HALF = 32

panel_A_frames = []

for _, row in rep_cell.iterrows():

    scene = int(
        row["scene"]
    )

    t = int(
        row["time_index"]
    )

    cx0 = int(
        round(
            float(row["centroid_x"])
        )
    )

    cy0 = int(
        round(
            float(row["centroid_y"])
        )
    )

    x0 = cx0 - CROP_HALF
    x1 = cx0 + CROP_HALF
    y0 = cy0 - CROP_HALF
    y1 = cy0 + CROP_HALF

    gfp = np.asarray(
        gfp_stack[
            scene,
            t,
            y0:y1,
            x0:x1,
        ],
        dtype=float,
    )

    dic = np.asarray(
        dic_stack[
            scene,
            t,
            y0:y1,
            x0:x1,
        ],
        dtype=float,
    )

    if (
        gfp.shape
        != (
            2 * CROP_HALF,
            2 * CROP_HALF,
        )
        or dic.shape != gfp.shape
    ):
        fail(
            f"Unexpected crop shape at scene={scene}, T={t}: "
            f"GFP={gfp.shape}, DIC={dic.shape}"
        )

    panel_A_frames.append(
        {
            "time_index": t,
            "biological_min": int(
                round(
                    float(
                        row["biological_min"]
                    )
                )
            ),
            "gfp": gfp,
            "dic": dic,
            "focus_x": float(
                row["focus_x"]
            ) - x0,
            "focus_y": float(
                row["focus_y"]
            ) - y0,
            "cytosol_x": float(
                row["cytosol_x"]
            ) - x0,
            "cytosol_y": float(
                row["cytosol_y"]
            ) - y0,
            "ratio": float(
                row["focus_to_cytosol_ratio"]
            ),
        }
    )


                                             
all_gfp = np.concatenate(
    [
        x["gfp"].ravel()
        for x in panel_A_frames
    ]
)

all_dic = np.concatenate(
    [
        x["dic"].ravel()
        for x in panel_A_frames
    ]
)

GFP_VMIN = float(
    np.quantile(
        all_gfp,
        0.05,
    )
)

GFP_VMAX = float(
    np.quantile(
        all_gfp,
        0.995,
    )
)

DIC_VMIN = float(
    np.quantile(
        all_dic,
        0.01,
    )
)

DIC_VMAX = float(
    np.quantile(
        all_dic,
        0.99,
    )
)

if GFP_VMAX <= GFP_VMIN:
    fail("Invalid GFP display range.")

if DIC_VMAX <= DIC_VMIN:
    fail("Invalid DIC display range.")


                                                            
rep_cell[
    [
        "scene",
        "track_id",
        "time_index",
        "biological_min",
        "centroid_y",
        "centroid_x",
        "focus_y",
        "focus_x",
        "cytosol_y",
        "cytosol_x",
        "focus_to_cytosol_ratio",
        "initial_focus_prominence_z",
    ]
].to_csv(
    CURATED
    / "Fig4_5_panelA_representative_real_cell.tsv",
    sep="\t",
    index=False,
)


                                                              
            
                                                              

manifest = []

for role, path in [
    (
        "manual_automated_time_aligned",
        TIME_ALIGNED,
    ),
    (
        "trajectory_concordance",
        CONCORDANCE,
    ),
    (
        "endpoint_distribution_summary",
        ENDPOINT,
    ),
    (
        "heldout_301_cell_trajectories",
        CELL_LONG,
    ),
    (
        "frozen_Guk1_recovery_phenotype",
        PHENO_CELL,
    ),
]:

    manifest.append(
        {
            "role": role,
            "path": str(path),
            "sha256": sha256(path),
        }
    )

for role, path in [
    (
        "190321_WT_frozen_GFP_projection",
        GFP_STACK,
    ),
    (
        "190321_WT_frozen_DIC_projection",
        DIC_STACK,
    ),
]:

    manifest.append(
        {
            "role": role,
            "path": str(path),
            "sha256": sha256(path),
        }
    )

pd.DataFrame(
    manifest
).to_csv(
    CURATED
    / "Fig4_5_source_manifest.tsv",
    sep="\t",
    index=False,
)


                                                              
              
                                                              

plt.rcParams.update(
    {
        "font.family":
            "DejaVu Sans",
        "font.size": 9.0,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "legend.fontsize": 7.8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


fig = plt.figure(
    figsize=(18.0, 11.6),
)

outer = fig.add_gridspec(
    2,
    12,
    height_ratios=[
        0.95,
        1.05,
    ],
    hspace=0.34,
    wspace=0.55,
)


                                                              
                                               
                                                              

Aouter = outer[
    0,
    0:6
].subgridspec(
    3,
    1,
    height_ratios=[
        0.11,
        0.82,
        0.07,
    ],
    hspace=0.03,
)


                                                              
                                  
                                                              

axAkey = fig.add_subplot(
    Aouter[0, 0]
)

axAkey.set_xlim(0, 1)
axAkey.set_ylim(0, 1)
axAkey.axis("off")

axAkey.text(
    0.00,
    0.90,
    "A",
    transform=axAkey.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)

               
focus_key = Circle(
    (0.20, 0.50),
    radius=0.025,
    transform=axAkey.transAxes,
    fill=False,
    edgecolor="#D14C9A",
    linewidth=2.0,
)

axAkey.add_patch(focus_key)

axAkey.plot(
    [0.20],
    [0.50],
    marker="o",
    markersize=2.8,
    color="#D14C9A",
    transform=axAkey.transAxes,
)

axAkey.text(
    0.235,
    0.50,
    "Focus measurement ROI",
    transform=axAkey.transAxes,
    ha="left",
    va="center",
    fontsize=8.2,
    color="#D14C9A",
)

                 
cyto_key = Circle(
    (0.61, 0.50),
    radius=0.025,
    transform=axAkey.transAxes,
    fill=False,
    edgecolor="#23A6B8",
    linewidth=1.8,
    linestyle="--",
)

axAkey.add_patch(cyto_key)

axAkey.plot(
    [0.61],
    [0.50],
    marker="o",
    markersize=2.6,
    color="#23A6B8",
    transform=axAkey.transAxes,
)

axAkey.text(
    0.645,
    0.50,
    "Cytosol reference ROI",
    transform=axAkey.transAxes,
    ha="left",
    va="center",
    fontsize=8.2,
    color="#23A6B8",
)


                                                              
            
                                                              

Agrid = Aouter[
    1,
    0
].subgridspec(
    1,
    3,
    wspace=0.045,
)

for i, frame in enumerate(
    panel_A_frames
):

    ax = fig.add_subplot(
        Agrid[0, i]
    )

                    
    ax.imshow(
        frame["dic"],
        cmap="gray",
        vmin=DIC_VMIN,
        vmax=DIC_VMAX,
        interpolation="nearest",
    )

                 
    g = np.clip(
        (
            frame["gfp"]
            - GFP_VMIN
        )
        /
        (
            GFP_VMAX
            - GFP_VMIN
        ),
        0,
        1,
    )

    alpha = np.clip(
        0.08
        + 0.78 * g,
        0.08,
        0.86,
    )

    ax.imshow(
        frame["gfp"],
        cmap="Greens",
        vmin=GFP_VMIN,
        vmax=GFP_VMAX,
        alpha=alpha,
        interpolation="nearest",
    )

                                               
                              
    focus_roi = Circle(
        (
            frame["focus_x"],
            frame["focus_y"],
        ),
        radius=4,
        fill=False,
        edgecolor="#D14C9A",
        linewidth=2.0,
    )

    cyto_roi = Circle(
        (
            frame["cytosol_x"],
            frame["cytosol_y"],
        ),
        radius=4,
        fill=False,
        edgecolor="#23A6B8",
        linewidth=1.8,
        linestyle="--",
    )

    ax.add_patch(focus_roi)
    ax.add_patch(cyto_roi)

                           
    ax.plot(
        frame["focus_x"],
        frame["focus_y"],
        marker="o",
        markersize=2.8,
        color="#D14C9A",
    )

    ax.plot(
        frame["cytosol_x"],
        frame["cytosol_y"],
        marker="o",
        markersize=2.6,
        color="#23A6B8",
    )

    ax.set_xlim(
        0,
        2 * CROP_HALF,
    )

    ax.set_ylim(
        2 * CROP_HALF,
        0,
    )

    ax.axis("off")

    ax.text(
        0.50,
        -0.045,
        f'{frame["biological_min"]} min',
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9.0,
    )


                                                              
                                         
                                                              

axAeq = fig.add_subplot(
    Aouter[2, 0]
)

axAeq.axis("off")

axAeq.text(
    0.50,
    0.50,
    r"$R_i(t)=I_{\mathrm{focus}}/I_{\mathrm{cytosol}}$",
    transform=axAeq.transAxes,
    ha="center",
    va="center",
    fontsize=9.3,
)

                                                              
                                        
                                                              

axB = fig.add_subplot(
    outer[
        0,
        6:9
    ]
)

axB.plot(
    plot_cmp["time"],
    plot_cmp["manual_raw"],
    marker="o",
    markersize=4.0,
    linewidth=1.7,
    label="Ръчно",
)

axB.plot(
    plot_cmp["time"],
    plot_cmp["automated_raw"],
    marker="o",
    markersize=4.0,
    linewidth=1.7,
    label="Автоматизирано",
)

axB.set_xlabel(
    "Време на възстановяване (min)"
)

axB.set_ylabel(
    "Фокус / цитозол"
)

axB.legend(
    frameon=False,
    loc="best",
)

axB.spines[
    "top"
].set_visible(False)

axB.spines[
    "right"
].set_visible(False)

axB.text(
    -0.20,
    1.04,
    "B",
    transform=axB.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                      
                                                              

axC = fig.add_subplot(
    outer[
        0,
        9:12
    ]
)

axC.axhline(
    1.0,
    linewidth=0.7,
    linestyle="--",
    color="#AAAAAA",
)

axC.plot(
    norm_cmp["time"],
    norm_cmp["manual_norm"],
    marker="o",
    markersize=4.0,
    linewidth=1.7,
    label="Ръчно",
)

axC.plot(
    norm_cmp["time"],
    norm_cmp["automated_norm"],
    marker="o",
    markersize=4.0,
    linewidth=1.7,
    label="Автоматизирано",
)

axC.set_xlabel(
    "Време на възстановяване (min)"
)

axC.set_ylabel(
    "Нормализирано към T1"
)

axC.legend(
    frameon=False,
    loc="best",
)

axC.spines[
    "top"
].set_visible(False)

axC.spines[
    "right"
].set_visible(False)

axC.text(
    -0.20,
    1.04,
    "C",
    transform=axC.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                        
                                                              

axD = fig.add_subplot(
    outer[
        1,
        0:7
    ]
)


for cell_id, g in cells.groupby(
    "_cell",
    sort=False,
):

    axD.plot(
        g["_time"],
        g["_Y"],
        linewidth=0.45,
        alpha=0.055,
        color="#555555",
    )


                      
med = (
    cells.groupby("_time")["_Y"]
    .median()
    .sort_index()
)

q25 = (
    cells.groupby("_time")["_Y"]
    .quantile(0.25)
    .sort_index()
)

q75 = (
    cells.groupby("_time")["_Y"]
    .quantile(0.75)
    .sort_index()
)


axD.fill_between(
    med.index.to_numpy(float),
    q25.to_numpy(float),
    q75.to_numpy(float),
    alpha=0.16,
)

axD.plot(
    med.index,
    med.values,
    linewidth=2.2,
    color="#222222",
)

axD.axhline(
    0,
    linewidth=0.7,
    linestyle="--",
    color="#AAAAAA",
)

axD.set_xlabel(
    "Време на възстановяване"
)

axD.set_ylabel(
    r"$Y_i(t)=\log_2[R_i(t)/R_i(T1)]$"
)

axD.spines[
    "top"
].set_visible(False)

axD.spines[
    "right"
].set_visible(False)

axD.text(
    -0.09,
    1.035,
    "D",
    transform=axD.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                          
                                                              

axE = fig.add_subplot(
    outer[
        1,
        7:12
    ]
)

bins = max(
    20,
    min(
        45,
        int(
            np.sqrt(
                len(phenotype)
            )
            * 2.2
        ),
    ),
)

axE.hist(
    phenotype,
    bins=bins,
    edgecolor="white",
    linewidth=0.5,
)

axE.axvline(
    0,
    linewidth=1.0,
    linestyle="--",
    color="#666666",
)

axE.set_xlabel(
    r"$P_i$: signed AUC, T1–T8"
)

axE.set_ylabel(
    "Брой клетки"
)

axE.spines[
    "top"
].set_visible(False)

axE.spines[
    "right"
].set_visible(False)

axE.text(
    -0.12,
    1.035,
    "E",
    transform=axE.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
        
                                                              

fig.subplots_adjust(
    left=0.065,
    right=0.985,
    top=0.975,
    bottom=0.085,
)

OUTBASE = (
    FINAL
    / "Fig4_5_Guk1_recovery_phenotype"
)

fig.savefig(
    str(OUTBASE) + ".png",
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE) + ".pdf",
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE) + ".svg",
    bbox_inches="tight",
)

plt.close(fig)


                                                              
                   
                                                              

legend = f"""Фигура 4.5. Автоматизираният Guk1-7-GFP показател възпроизвежда кинетиката на постстресовото възстановяване.

(A) Една реална проследена WT клетка от frozen 190321_WT held-out анализа при 10, 60 и 120 min от периода на възстановяване. DIC изображението показва клетъчната морфология, а зеленият overlay — Guk1-7-GFP сигнала. Магентовият и прекъснатият cyan кръг представляват измервателни ROI, а не граници на биологични структури: първата ROI е инициализирана върху Guk1 фокуса при T1, а втората служи като цитозолна референция. След T1 двете позиции се пренасят единствено според проследеното движение на клетката, без повторно търсене на най-ярката флуоресцентна област.

(B) Сравнение между публикуваното ръчно и автоматизираното ненормализирано съотношение фокус/цитозол. Въпреки систематичната разлика в амплитудата времевото подреждане е силно съгласувано (Spearman ρ={raw_rho:.3f}; Pearson r={raw_r:.3f}).

(C) Същите времеви профили след нормализация спрямо първата времева точка. Съответствието се увеличава до Spearman ρ={norm_rho:.3f} и Pearson r={norm_r:.3f}, което показва, че автоматизираният анализ надеждно възпроизвежда относителната кинетика на възстановяването.

(D) Нормализирани едноклетъчни траектории за {n_cells} WT клетки с пълно проследяване. Тънките линии представят отделните клетки, плътната линия — медианата, а затъмнената област — интерквартилния диапазон. При {lower}/{n_cells} клетки ({100*lower/n_cells:.2f}%) съотношението фокус/цитозол е по-ниско в крайната спрямо началната времева точка.

(E) Разпределение на frozen количествения фенотип P_i, дефиниран като знакова площ под нормализираната log₂ траектория за T1–T8. Отрицателните стойности означават нетно намаляване на Guk1 фокусното обогатяване; по-отрицателните стойности отразяват по-силно или по-продължително изчистване.

Валидацията оценява измервателния метод, а не представлява независима биологична репликация, тъй като публикуваният ръчен и автоматизираният анализ на валидационния етап произлизат от една и съща биологична серия.
"""

legend_path = (
    LEGENDS
    / "Fig4_5_legend_bg.txt"
)

legend_path.write_text(
    legend,
    encoding="utf-8",
)


                                                              
                 
                                                              

print()
print("===== FIGURE 4.5 COMPLETE =====")
print(
    f"raw:  rho={raw_rho:.4f}, r={raw_r:.4f}"
)
print(
    f"norm: rho={norm_rho:.4f}, r={norm_r:.4f}"
)
print(
    f"cells: {n_cells}"
)
print(
    f"endpoint lower: {lower}/{n_cells} "
    f"({100*lower/n_cells:.2f}%)"
)
print(
    f"phenotype column: {auc_col}"
)
print(
    f"phenotype n: {len(phenotype)}"
)
print(
    "Panel A representative real cell: "
    f"scene={REP_SCENE}, track_id={REP_TRACK}"
)

print(
    "Panel A source: frozen 190321_WT "
    "egfp_max_projection.npy + dic_best_gfp_z.npy"
)
print()
print(str(OUTBASE) + ".pdf")
