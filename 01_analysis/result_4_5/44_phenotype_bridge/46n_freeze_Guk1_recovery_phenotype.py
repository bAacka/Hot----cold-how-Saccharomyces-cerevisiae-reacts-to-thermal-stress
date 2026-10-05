#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress


ROOT = Path("/")
OUT = ROOT / "44_phenotype_bridge"
STATE = OUT / "46_phenotype_state"

DEST = STATE / "46n_Guk1_recovery_phenotype"
DEST.mkdir(parents=True, exist_ok=True)

AUDIT = (
    STATE /
    "46m_replicate_comparison_T1_T8"
)

DATASETS = {
    "190319":
        STATE /
        "190319_tracked_focus_T0_T8",

    "190326":
        STATE /
        "190326_tracked_focus_T0_T8",
}

PRIMARY_RADIUS = 4

                                                          
PRIMARY_TIMES = list(range(1, 9))


def require_columns(x, cols, label):
    missing = set(cols) - set(x.columns)

    if missing:
        raise RuntimeError(
            f"{label}: missing columns "
            f"{sorted(missing)}"
        )


                                                              
                              
                                                              

rows = []

for acquisition, base in DATASETS.items():

    path = (
        base /
        "44g_focus_cytosol_measurements.tsv.gz"
    )

    if not path.is_file():
        raise FileNotFoundError(path)

    x = pd.read_csv(
        path,
        sep="\t",
    )

    require_columns(
        x,
        [
            "scene",
            "baseline_label",
            "time_index",
            "radius",
            "ratio",
        ],
        acquisition,
    )

    x = x[
        x["radius"].eq(PRIMARY_RADIUS)
        &
        x["time_index"].isin(PRIMARY_TIMES)
    ].copy()

    cells = (
        x[
            [
                "scene",
                "baseline_label",
            ]
        ]
        .drop_duplicates()
    )

    if len(cells) != 20:
        raise RuntimeError(
            f"{acquisition}: expected 20 selected cells, "
            f"observed {len(cells)}"
        )

    for (
        scene,
        baseline_label
    ), g in x.groupby(
        [
            "scene",
            "baseline_label",
        ]
    ):

        g = g.sort_values(
            "time_index"
        )

        times = (
            g["time_index"]
            .astype(int)
            .tolist()
        )

        if times != PRIMARY_TIMES:
            raise RuntimeError(
                f"{acquisition} "
                f"S{scene} cell={baseline_label}: "
                f"expected {PRIMARY_TIMES}, got {times}"
            )

        t = (
            g["time_index"]
            .to_numpy(float)
        )

        r = (
            g["ratio"]
            .to_numpy(float)
        )

        if (
            (~np.isfinite(r)).any()
            or
            (r <= 0).any()
        ):
            raise RuntimeError(
                f"{acquisition} "
                f"S{scene} cell={baseline_label}: "
                "non-finite/non-positive ratio"
            )

        r1 = float(r[0])

        logrel = np.log2(
            r / r1
        )

                            
        auc = float(
            np.trapezoid(
                logrel,
                t,
            )
        )

                                     
        fit = linregress(
            t,
            logrel,
        )

        endpoint = float(
            logrel[-1]
        )

        jmin = int(
            np.argmin(logrel)
        )

        rows.append({
            "acquisition":
                acquisition,

            "scene":
                int(scene),

            "baseline_label":
                int(baseline_label),

                                            
                     
                                            
            "primary_Guk1_recovery_AUC_log2_rel_T1":
                auc,

                                             
                                 
                                       
                                      
            "primary_AUC_mean_equivalent":
                auc / 7.0,

                                            
                              
            "derived_clearance_score":
                -auc,

                                            
                       
                                            
            "secondary_log2_T8_over_T1":
                endpoint,

            "secondary_T8_over_T1":
                float(r[-1] / r1),

            "secondary_linear_slope":
                float(fit.slope),

            "secondary_linear_slope_r2":
                float(fit.rvalue ** 2),

            "secondary_minimum_log2_rel_T1":
                float(logrel[jmin]),

            "secondary_minimum_time_index":
                int(t[jmin]),

                        
            "T1_ratio":
                r1,

            "T8_ratio":
                float(r[-1]),

            "n_primary_timepoints":
                len(t),

            "primary_radius_px":
                PRIMARY_RADIUS,
        })


cell = pd.DataFrame(
    rows
)

cell = cell.sort_values(
    [
        "acquisition",
        "scene",
        "baseline_label",
    ]
).reset_index(drop=True)

if len(cell) != 40:
    raise RuntimeError(
        f"Expected 40 cells total, observed {len(cell)}"
    )


                                                              
                                                    
 
                                                          
                                                      
                                                              

audit_path = (
    AUDIT /
    "46m_candidate_cell_phenotypes.tsv"
)

if not audit_path.is_file():
    raise FileNotFoundError(
        audit_path
    )

old = pd.read_csv(
    audit_path,
    sep="\t",
    dtype={
        "acquisition": str,
    },
)

                                                            
cell["acquisition"] = (
    cell["acquisition"]
    .astype(str)
)

old["acquisition"] = (
    old["acquisition"]
    .astype(str)
)

require_columns(
    old,
    [
        "acquisition",
        "scene",
        "baseline_label",
        "auc_log2_relative_T1_T1_T8",
        "log2_T8_over_T1",
        "linear_slope_log2_relative_T1",
    ],
    "46m audit",
)

check = cell.merge(
    old[
        [
            "acquisition",
            "scene",
            "baseline_label",
            "auc_log2_relative_T1_T1_T8",
            "log2_T8_over_T1",
            "linear_slope_log2_relative_T1",
        ]
    ],
    on=[
        "acquisition",
        "scene",
        "baseline_label",
    ],
    how="outer",
    validate="one_to_one",
    indicator=True,
)

if not check["_merge"].eq("both").all():
    raise RuntimeError(
        "46m/46n cell identity mismatch"
    )

check[
    "abs_delta_primary_AUC"
] = np.abs(
    check[
        "primary_Guk1_recovery_AUC_log2_rel_T1"
    ]
    -
    check[
        "auc_log2_relative_T1_T1_T8"
    ]
)

check[
    "abs_delta_endpoint"
] = np.abs(
    check[
        "secondary_log2_T8_over_T1"
    ]
    -
    check[
        "log2_T8_over_T1"
    ]
)

check[
    "abs_delta_slope"
] = np.abs(
    check[
        "secondary_linear_slope"
    ]
    -
    check[
        "linear_slope_log2_relative_T1"
    ]
)

tol = 1e-10

for c in [
    "abs_delta_primary_AUC",
    "abs_delta_endpoint",
    "abs_delta_slope",
]:
    if check[c].max() > tol:
        raise RuntimeError(
            f"Freeze does not reproduce 46m: "
            f"{c} max={check[c].max()}"
        )


                                                              
                        
 
                                                     
                                                       
                        
                                                              

PRIMARY = (
    "primary_Guk1_recovery_AUC_log2_rel_T1"
)


def summarize(df, groups):

    z = (
        df.groupby(
            groups,
            as_index=False,
        )
        .agg(
            n_cells=(
                PRIMARY,
                "size",
            ),

            mean_primary_AUC=(
                PRIMARY,
                "mean",
            ),

            median_primary_AUC=(
                PRIMARY,
                "median",
            ),

            sd_primary_AUC=(
                PRIMARY,
                "std",
            ),

            q25_primary_AUC=(
                PRIMARY,
                lambda q:
                    q.quantile(0.25),
            ),

            q75_primary_AUC=(
                PRIMARY,
                lambda q:
                    q.quantile(0.75),
            ),

            mean_log2_T8_over_T1=(
                "secondary_log2_T8_over_T1",
                "mean",
            ),

            median_log2_T8_over_T1=(
                "secondary_log2_T8_over_T1",
                "median",
            ),

            mean_linear_slope=(
                "secondary_linear_slope",
                "mean",
            ),
        )
    )

    return z


scene = summarize(
    cell,
    [
        "acquisition",
        "scene",
    ],
)

acquisition = summarize(
    cell,
    [
        "acquisition",
    ],
)


                                                              
               
                                                              

cell.to_csv(
    DEST /
    "46n_frozen_cell_Guk1_recovery_phenotype.tsv",
    sep="\t",
    index=False,
)

scene.to_csv(
    DEST /
    "46n_scene_summary.tsv",
    sep="\t",
    index=False,
)

acquisition.to_csv(
    DEST /
    "46n_acquisition_summary.tsv",
    sep="\t",
    index=False,
)

check.to_csv(
    DEST /
    "46n_validation_against_46m.tsv",
    sep="\t",
    index=False,
)


                                                              
                      
                                                              

spec = """\
GUK1 SINGLE-CELL RECOVERY PHENOTYPE — FROZEN SPECIFICATION
===========================================================

STATUS
------
Frozen after measurement validation in two independent ssa1-2DD
acquisitions (190319 and 190326), and before testing any molecular
perturbation -> phenotype association.

This is a phenotype-definition freeze, not a preregistered inferential
endpoint.

MEASUREMENT
-----------
Input:
tracked single-cell Guk1 focus/cytosol ratio.

Primary ROI:
radius = 4 px.

Sensitivity radii 3 and 5 remain technical sensitivity analyses only.

CELL SELECTION
--------------
Five baseline-focus cells are selected per scene from complete tracks
using T0 focus/cytosol ratio.

T0 is used for track initialization and baseline-focus selection only.

T0 is NOT part of the primary phenotype.

PRIMARY WINDOW
--------------
T1 through T8 inclusive.

This corresponds approximately to 10-80 min recovery.

The window was established by the phenotype-blind technical-QC process
documented in 46e and 46j.

PRIMARY CELL-LEVEL PHENOTYPE
----------------------------
For cell i:

    R_i(t) = radius-4 Guk1 focus/cytosol ratio at time t

    Y_i(t) = log2( R_i(t) / R_i(T1) )

Primary phenotype:

    P_i = trapezoidal AUC of Y_i(t), t = T1,...,T8

Time is represented by the equally ordered acquisition frame index.

Interpretation:

    P_i < 0:
        net loss of focus enrichment relative to T1.

    More negative P_i:
        greater / more sustained Guk1 focus clearance over T1-T8.

    P_i near 0:
        persistence near the T1 state.

The raw signed AUC is the scientific primary phenotype.

A sign-flipped value (-P_i) may be displayed as a convenient
"clearance score", but it is not a different phenotype.

WHY AUC
-------
AUC was an explicitly considered phenotype before this freeze.

It is selected because:

1. it uses the complete prespecified T1-T8 trajectory;
2. it accommodates the observed non-monotonic late rebound;
3. it is less dependent on one terminal frame than T8/T1;
4. it does not impose a globally linear recovery model;
5. it showed close replication across the two independent acquisitions;
6. no molecular predictor or perturbation result was used to choose it.

SECONDARY DESCRIPTORS
---------------------
Not primary:

- log2(T8/T1)
- T8/T1
- linear slope of log2-relative trajectory
- linear-slope R-squared
- minimum log2-relative value
- time of minimum

These are retained for interpretation and sensitivity analysis.

HIERARCHY / INFERENCE
---------------------
Cells are nested within scenes.
Scenes are nested within acquisitions.

Cells must NOT be treated as independent biological replicates.

The two current acquisitions validate the measurement/phenotype
coordinate; they do not provide a large-N biological replication study.

Future perturbational analyses must preserve the experimental hierarchy
and use acquisition/experiment as the biological replication level.

PROHIBITIONS
------------
Do not:

- reintroduce T0 into the primary AUC;
- change the primary radius from 4 px based on downstream results;
- alter T1-T8 based on molecular association strength;
- replace AUC with endpoint or slope because another metric gives a
  stronger perturbation effect;
- treat individual cells as independent biological replicates.

Any future change to this phenotype definition must be versioned as a
new phenotype rather than silently modifying this freeze.
"""

(
    DEST /
    "46n_GUK1_RECOVERY_PHENOTYPE_FREEZE.txt"
).write_text(
    spec
)


                                                              
                       
                                                              

print()
print("=" * 100)
print("FROZEN PRIMARY PHENOTYPE")
print("=" * 100)

print(
    "AUC_T1_T8[ log2(R_t / R_T1) ], "
    "primary radius = 4 px"
)

print()
print("=" * 100)
print("ACQUISITION SUMMARIES — CELL-LEVEL PRIMARY PHENOTYPE")
print("=" * 100)

print(
    acquisition.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6g}",
    )
)

print()
print("=" * 100)
print("SCENE SUMMARIES")
print("=" * 100)

print(
    scene.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6g}",
    )
)

print()
print("=" * 100)
print("46m -> 46n REPRODUCTION CHECK")
print("=" * 100)

print(
    "max |delta primary AUC| =",
    check[
        "abs_delta_primary_AUC"
    ].max(),
)

print(
    "max |delta endpoint| =",
    check[
        "abs_delta_endpoint"
    ].max(),
)

print(
    "max |delta slope| =",
    check[
        "abs_delta_slope"
    ].max(),
)

print()
print("OUTPUT:", DEST)
