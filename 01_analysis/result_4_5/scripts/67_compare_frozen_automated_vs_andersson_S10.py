#!/usr/bin/env python3

from pathlib import Path
import hashlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import (
    pearsonr,
    spearmanr,
)


ROOT = Path(
    "/"
)

BASE = (
    ROOT
    / "44_phenotype_bridge"
)

AUTO = (
    BASE
    / "focus_metric_validation"
    / "190321_WT"
    / "frozen_v3"
)

REF = (
    BASE
    / "reference_validation"
    / "andersson_2021"
    / "parsed_WT"
)

OUT = (
    BASE
    / "focus_metric_validation"
    / "190321_WT"
    / "formal_manual_comparison"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
        
                                                              

MANUAL_RAW = (
    REF
    / "ANDERSSON_S10_WT_REFERENCE_TRAJECTORY.tsv"
)

MANUAL_NORM = (
    REF
    / "ANDERSSON_S10_WT_REFERENCE_TRAJECTORY_NORMALIZED.tsv"
)

MANUAL_ENDPOINT = (
    REF
    / "ANDERSSON_S10_WT_ENDPOINT_BY_CELL.tsv"
)

MANUAL_HASHES = (
    REF
    / "SHA256SUMS.txt"
)

AUTO_POOLED = (
    AUTO
    / "HELDOUT_POOLED_TRAJECTORY.tsv"
)

AUTO_SCENE = (
    AUTO
    / "HELDOUT_PER_SCENE_TRAJECTORY.tsv"
)

AUTO_LONG = (
    AUTO
    / "HELDOUT_FOCUS_LONG_COMPLETE_T1_T12.tsv"
)

AUTO_ENDPOINT = (
    AUTO
    / "HELDOUT_ENDPOINT_BY_CELL.tsv"
)

AUTO_HASHES = (
    AUTO
    / "HELDOUT_AUTOMATED_PRECOMPARISON_SHA256SUMS.txt"
)


                                                              
                   
                                                              

def file_sha256(
    path,
):
    return hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


def parse_manifest_by_basename(
    manifest,
):

    out = {}

    for raw in manifest.read_text().splitlines():

        raw = raw.strip()

        if not raw:
            continue

        parts = raw.split(
            None,
            1,
        )

        if len(parts) != 2:
            raise RuntimeError(
                f"Cannot parse hash line:\n{raw}"
            )

        digest = parts[0].strip()
        filename = parts[1].strip()

        if filename.startswith("*"):
            filename = filename[1:]

        basename = Path(
            filename
        ).name

        if basename in out:
            raise RuntimeError(
                f"Duplicate basename in hash manifest: {basename}"
            )

        out[
            basename
        ] = digest

    return out


def verify_against_manifest(
    manifest,
    paths,
    label,
):

    expected = parse_manifest_by_basename(
        manifest
    )

    print()
    print("=" * 120)
    print(
        f"{label} HASH VERIFICATION"
    )
    print("=" * 120)

    for p in paths:

        if not p.exists():
            raise FileNotFoundError(
                p
            )

        name = p.name

        if name not in expected:
            raise RuntimeError(
                f"{name} is absent from {manifest}"
            )

        observed = file_sha256(
            p
        )

        wanted = expected[
            name
        ]

        ok = (
            observed
            == wanted
        )

        print(
            f"{name}:",
            "VERIFIED"
            if ok
            else "MISMATCH",
        )

        if not ok:
            raise RuntimeError(
                f"Hash mismatch for {p}\n"
                f"Expected: {wanted}\n"
                f"Observed: {observed}"
            )


verify_against_manifest(
    MANUAL_HASHES,
    [
        MANUAL_RAW,
        MANUAL_NORM,
        MANUAL_ENDPOINT,
    ],
    "MANUAL REFERENCE",
)

verify_against_manifest(
    AUTO_HASHES,
    [
        AUTO_POOLED,
        AUTO_SCENE,
        AUTO_LONG,
        AUTO_ENDPOINT,
    ],
    "AUTOMATED HELD-OUT",
)


                                                              
                    
                                                              

manual_raw = pd.read_csv(
    MANUAL_RAW,
    sep="\t",
)

manual_norm = pd.read_csv(
    MANUAL_NORM,
    sep="\t",
)

manual_endpoint = pd.read_csv(
    MANUAL_ENDPOINT,
    sep="\t",
)

auto = pd.read_csv(
    AUTO_POOLED,
    sep="\t",
)

auto_scene = pd.read_csv(
    AUTO_SCENE,
    sep="\t",
)

auto_long = pd.read_csv(
    AUTO_LONG,
    sep="\t",
)

auto_endpoint = pd.read_csv(
    AUTO_ENDPOINT,
    sep="\t",
)


                                                              
                      
                                                              

required_manual_raw = {
    "time_min",
    "n",
    "mean",
    "median",
    "q25",
    "q75",
}

required_manual_norm = {
    "time_min",
    "n",
    "mean",
    "median",
    "q25",
    "q75",
}

required_manual_endpoint = {
    "cell_id",
    "ratio_10",
    "ratio_120",
    "relative_120_over_10",
    "decreased",
}

required_auto = {
    "biological_min",
    "n_cells",
    "median_ratio",
    "mean_ratio",
    "q25_ratio",
    "q75_ratio",
    "median_ratio_relative_T1",
}

required_auto_scene = {
    "scene",
    "biological_min",
    "median_ratio",
    "median_ratio_relative_T1",
}

required_auto_long = {
    "scene",
    "cell_id",
    "biological_min",
    "ratio_relative_to_T1",
}

required_auto_endpoint = {
    "cell_id",
    "scene",
    "focus_to_cytosol_ratio_T1",
    "focus_to_cytosol_ratio_T12",
    "ratio_relative_to_T1_T12",
    "T12_lower_than_T1",
}


checks = [
    (
        manual_raw,
        required_manual_raw,
        "manual raw trajectory",
    ),
    (
        manual_norm,
        required_manual_norm,
        "manual normalized trajectory",
    ),
    (
        manual_endpoint,
        required_manual_endpoint,
        "manual endpoint",
    ),
    (
        auto,
        required_auto,
        "automated pooled trajectory",
    ),
    (
        auto_scene,
        required_auto_scene,
        "automated per-scene trajectory",
    ),
    (
        auto_long,
        required_auto_long,
        "automated long table",
    ),
    (
        auto_endpoint,
        required_auto_endpoint,
        "automated endpoint",
    ),
]


for df, required, name in checks:

    missing = (
        required
        - set(
            df.columns
        )
    )

    if missing:
        raise RuntimeError(
            f"{name}: missing columns {sorted(missing)}"
        )


                                                              
                      
                                                              

expected_times = list(
    range(
        10,
        121,
        10,
    )
)


for df, col, name in [
    (
        manual_raw,
        "time_min",
        "manual raw",
    ),
    (
        manual_norm,
        "time_min",
        "manual normalized",
    ),
    (
        auto,
        "biological_min",
        "automated",
    ),
]:

    observed = sorted(
        int(x)
        for x in df[
            col
        ].unique()
    )

    if observed != expected_times:
        raise RuntimeError(
            f"{name} time axis differs from 10-120 min:\n"
            f"{observed}"
        )


mraw = (
    manual_raw[
        [
            "time_min",
            "n",
            "median",
            "mean",
            "q25",
            "q75",
        ]
    ]
    .rename(
        columns={
            "n":
                "manual_n",

            "median":
                "manual_median_ratio",

            "mean":
                "manual_mean_ratio",

            "q25":
                "manual_q25_ratio",

            "q75":
                "manual_q75_ratio",
        }
    )
)


mnorm = (
    manual_norm[
        [
            "time_min",
            "median",
            "mean",
            "q25",
            "q75",
        ]
    ]
    .rename(
        columns={
            "median":
                "manual_median_relative_T1",

            "mean":
                "manual_mean_relative_T1",

            "q25":
                "manual_q25_relative_T1",

            "q75":
                "manual_q75_relative_T1",
        }
    )
)


a = (
    auto[
        [
            "biological_min",
            "n_cells",
            "median_ratio",
            "mean_ratio",
            "q25_ratio",
            "q75_ratio",
            "median_ratio_relative_T1",
        ]
    ]
    .rename(
        columns={
            "biological_min":
                "time_min",

            "n_cells":
                "automated_n",

            "median_ratio":
                "automated_median_ratio",

            "mean_ratio":
                "automated_mean_ratio",

            "q25_ratio":
                "automated_q25_ratio",

            "q75_ratio":
                "automated_q75_ratio",

            "median_ratio_relative_T1":
                "automated_median_relative_T1",
        }
    )
)


aligned = (
    mraw
    .merge(
        mnorm,
        on="time_min",
        how="inner",
        validate="1:1",
    )
    .merge(
        a,
        on="time_min",
        how="inner",
        validate="1:1",
    )
    .sort_values(
        "time_min"
    )
    .reset_index(
        drop=True
    )
)


if len(
    aligned
) != 12:

    raise RuntimeError(
        f"Expected 12 aligned timepoints; got {len(aligned)}"
    )


aligned[
    "raw_automated_minus_manual"
] = (
    aligned[
        "automated_median_ratio"
    ]
    - aligned[
        "manual_median_ratio"
    ]
)

aligned[
    "normalized_automated_minus_manual"
] = (
    aligned[
        "automated_median_relative_T1"
    ]
    - aligned[
        "manual_median_relative_T1"
    ]
)


aligned.to_csv(
    OUT
    / "TIME_ALIGNED_COMPARISON.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                         
                                                              

auto_norm_dist = (
    auto_long.groupby(
        "biological_min"
    )[
        "ratio_relative_to_T1"
    ]
    .agg(
        automated_norm_n="count",
        automated_norm_median="median",
        automated_norm_q25=lambda x:
            x.quantile(
                0.25
            ),
        automated_norm_q75=lambda x:
            x.quantile(
                0.75
            ),
    )
    .reset_index()
    .rename(
        columns={
            "biological_min":
                "time_min"
        }
    )
)


                                                              
                    
 
                                                 
                                                                 
                          
                                                              

def trajectory_metrics(
    comparison,
    manual_col,
    automated_col,
    scale,
):

    x = comparison[
        manual_col
    ].to_numpy(
        dtype=float
    )

    y = comparison[
        automated_col
    ].to_numpy(
        dtype=float
    )

    t = comparison[
        "time_min"
    ].to_numpy(
        dtype=float
    )


    pr = pearsonr(
        x,
        y,
    )

    sr = spearmanr(
        x,
        y,
    )


    error = (
        y - x
    )


    duration = float(
        t.max()
        - t.min()
    )


    manual_auc_mean = float(
        np.trapezoid(
            x,
            t,
        )
        / duration
    )

    automated_auc_mean = float(
        np.trapezoid(
            y,
            t,
        )
        / duration
    )


    return {
        "scale":
            scale,

        "n_timepoints":
            len(
                x
            ),

        "pearson_r":
            float(
                pr.statistic
            ),

        "pearson_p_descriptive_only":
            float(
                pr.pvalue
            ),

        "spearman_rho":
            float(
                sr.statistic
            ),

        "spearman_p_descriptive_only":
            float(
                sr.pvalue
            ),

        "MAE":
            float(
                np.mean(
                    np.abs(
                        error
                    )
                )
            ),

        "RMSE":
            float(
                np.sqrt(
                    np.mean(
                        error ** 2
                    )
                )
            ),

        "mean_bias_automated_minus_manual":
            float(
                np.mean(
                    error
                )
            ),

        "manual_time_averaged_AUC":
            manual_auc_mean,

        "automated_time_averaged_AUC":
            automated_auc_mean,

        "AUC_difference_automated_minus_manual":
            (
                automated_auc_mean
                - manual_auc_mean
            ),

        "manual_start":
            float(
                x[
                    0
                ]
            ),

        "automated_start":
            float(
                y[
                    0
                ]
            ),

        "manual_end":
            float(
                x[
                    -1
                ]
            ),

        "automated_end":
            float(
                y[
                    -1
                ]
            ),

        "endpoint_difference_automated_minus_manual":
            float(
                y[
                    -1
                ]
                - x[
                    -1
                ]
            ),
    }


metrics = pd.DataFrame(
    [
        trajectory_metrics(
            aligned,
            "manual_median_ratio",
            "automated_median_ratio",
            "raw_focus_to_cytosol_ratio",
        ),

        trajectory_metrics(
            aligned,
            "manual_median_relative_T1",
            "automated_median_relative_T1",
            "within_cell_T1_normalized_ratio",
        ),
    ]
)


metrics.to_csv(
    OUT
    / "TRAJECTORY_CONCORDANCE.tsv",
    sep="\t",
    index=False,
)


                                                              
                       
                                                              

scene_rows = []


for scene, sg in auto_scene.groupby(
    "scene"
):

    sg = (
        sg[
            [
                "biological_min",
                "median_ratio",
                "median_ratio_relative_T1",
            ]
        ]
        .rename(
            columns={
                "biological_min":
                    "time_min",

                "median_ratio":
                    "scene_median_ratio",

                "median_ratio_relative_T1":
                    "scene_median_relative_T1",
            }
        )
    )


    z = (
        aligned[
            [
                "time_min",
                "manual_median_ratio",
                "manual_median_relative_T1",
            ]
        ]
        .merge(
            sg,
            on="time_min",
            how="inner",
            validate="1:1",
        )
    )


    raw_r = pearsonr(
        z[
            "manual_median_ratio"
        ],
        z[
            "scene_median_ratio"
        ],
    )

    raw_s = spearmanr(
        z[
            "manual_median_ratio"
        ],
        z[
            "scene_median_ratio"
        ],
    )

    norm_r = pearsonr(
        z[
            "manual_median_relative_T1"
        ],
        z[
            "scene_median_relative_T1"
        ],
    )

    norm_s = spearmanr(
        z[
            "manual_median_relative_T1"
        ],
        z[
            "scene_median_relative_T1"
        ],
    )


    raw_error = (
        z[
            "scene_median_ratio"
        ].to_numpy()
        - z[
            "manual_median_ratio"
        ].to_numpy()
    )

    norm_error = (
        z[
            "scene_median_relative_T1"
        ].to_numpy()
        - z[
            "manual_median_relative_T1"
        ].to_numpy()
    )


    scene_rows.append({
        "scene":
            int(
                scene
            ),

        "raw_pearson_r":
            float(
                raw_r.statistic
            ),

        "raw_spearman_rho":
            float(
                raw_s.statistic
            ),

        "raw_MAE":
            float(
                np.mean(
                    np.abs(
                        raw_error
                    )
                )
            ),

        "normalized_pearson_r":
            float(
                norm_r.statistic
            ),

        "normalized_spearman_rho":
            float(
                norm_s.statistic
            ),

        "normalized_MAE":
            float(
                np.mean(
                    np.abs(
                        norm_error
                    )
                )
            ),

        "scene_T120_ratio":
            float(
                z.loc[
                    z[
                        "time_min"
                    ] == 120,
                    "scene_median_ratio",
                ].iloc[
                    0
                ]
            ),

        "scene_T120_relative_T1":
            float(
                z.loc[
                    z[
                        "time_min"
                    ] == 120,
                    "scene_median_relative_T1",
                ].iloc[
                    0
                ]
            ),
    })


scene_metrics = pd.DataFrame(
    scene_rows
)

scene_metrics.to_csv(
    OUT
    / "PER_SCENE_MANUAL_CONCORDANCE.tsv",
    sep="\t",
    index=False,
)


                                                              
                        
 
                   
 
                                                         
                                                            
                                                              

manual_rel = manual_endpoint[
    "relative_120_over_10"
].astype(
    float
)

auto_rel = auto_endpoint[
    "ratio_relative_to_T1_T12"
].astype(
    float
)


endpoint_summary = pd.DataFrame(
    [
        {
            "source":
                "Andersson_S10_manual",

            "n_cells":
                len(
                    manual_rel
                ),

            "median_T120_over_T10":
                float(
                    manual_rel.median()
                ),

            "q25_T120_over_T10":
                float(
                    manual_rel.quantile(
                        0.25
                    )
                ),

            "q75_T120_over_T10":
                float(
                    manual_rel.quantile(
                        0.75
                    )
                ),

            "fraction_decreased":
                float(
                    manual_endpoint[
                        "decreased"
                    ].astype(
                        bool
                    ).mean()
                ),
        },

        {
            "source":
                "frozen_automated_heldout",

            "n_cells":
                len(
                    auto_rel
                ),

            "median_T120_over_T10":
                float(
                    auto_rel.median()
                ),

            "q25_T120_over_T10":
                float(
                    auto_rel.quantile(
                        0.25
                    )
                ),

            "q75_T120_over_T10":
                float(
                    auto_rel.quantile(
                        0.75
                    )
                ),

            "fraction_decreased":
                float(
                    auto_endpoint[
                        "T12_lower_than_T1"
                    ].astype(
                        bool
                    ).mean()
                ),
        },
    ]
)


endpoint_summary.to_csv(
    OUT
    / "ENDPOINT_DISTRIBUTION_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
                        
                                                              

fig, ax = plt.subplots(
    figsize=(
        8,
        5,
    )
)


ax.plot(
    aligned[
        "time_min"
    ],
    aligned[
        "manual_median_ratio"
    ],
    marker="o",
    label="Andersson S10 manual",
)

ax.fill_between(
    aligned[
        "time_min"
    ],
    aligned[
        "manual_q25_ratio"
    ],
    aligned[
        "manual_q75_ratio"
    ],
    alpha=0.15,
)


ax.plot(
    aligned[
        "time_min"
    ],
    aligned[
        "automated_median_ratio"
    ],
    marker="o",
    label="Frozen automated",
)

ax.fill_between(
    aligned[
        "time_min"
    ],
    aligned[
        "automated_q25_ratio"
    ],
    aligned[
        "automated_q75_ratio"
    ],
    alpha=0.15,
)


ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xlabel(
    "Recovery time (min)"
)

ax.set_ylabel(
    "Focus / cytosol ratio"
)

ax.set_title(
    "190321 WT — manual vs frozen automated recovery phenotype"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUT
    / "RAW_TRAJECTORY_MANUAL_VS_AUTOMATED.png",
    dpi=200,
)

plt.close(
    fig
)


                                                              
                               
                                                              

norm_plot = (
    aligned[
        [
            "time_min",
            "manual_median_relative_T1",
            "manual_q25_relative_T1",
            "manual_q75_relative_T1",
            "automated_median_relative_T1",
        ]
    ]
    .merge(
        auto_norm_dist,
        on="time_min",
        how="inner",
        validate="1:1",
    )
)


fig, ax = plt.subplots(
    figsize=(
        8,
        5,
    )
)


ax.plot(
    norm_plot[
        "time_min"
    ],
    norm_plot[
        "manual_median_relative_T1"
    ],
    marker="o",
    label="Andersson S10 manual",
)

ax.fill_between(
    norm_plot[
        "time_min"
    ],
    norm_plot[
        "manual_q25_relative_T1"
    ],
    norm_plot[
        "manual_q75_relative_T1"
    ],
    alpha=0.15,
)


ax.plot(
    norm_plot[
        "time_min"
    ],
    norm_plot[
        "automated_median_relative_T1"
    ],
    marker="o",
    label="Frozen automated",
)

ax.fill_between(
    norm_plot[
        "time_min"
    ],
    norm_plot[
        "automated_norm_q25"
    ],
    norm_plot[
        "automated_norm_q75"
    ],
    alpha=0.15,
)


ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xlabel(
    "Recovery time (min)"
)

ax.set_ylabel(
    "Within-cell ratio relative to 10 min"
)

ax.set_title(
    "190321 WT — normalized manual vs automated recovery"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUT
    / "NORMALIZED_TRAJECTORY_MANUAL_VS_AUTOMATED.png",
    dpi=200,
)

plt.close(
    fig
)


                                                              
                              
                                                              

raw = metrics.loc[
    metrics[
        "scale"
    ]
    == "raw_focus_to_cytosol_ratio"
].iloc[
    0
]

norm = metrics.loc[
    metrics[
        "scale"
    ]
    == "within_cell_T1_normalized_ratio"
].iloc[
    0
]


manual_ep = endpoint_summary.loc[
    endpoint_summary[
        "source"
    ]
    == "Andersson_S10_manual"
].iloc[
    0
]

auto_ep = endpoint_summary.loc[
    endpoint_summary[
        "source"
    ]
    == "frozen_automated_heldout"
].iloc[
    0
]


text = f"""FORMAL MANUAL VS AUTOMATED PHENOTYPE COMPARISON

DATASET
190321 WT.

COMPARISON
Frozen automated phenotype versus published Andersson S10 manual
focus/cytosol measurements.

TIME ALIGNMENT
Exact 10-120 min alignment, 12 timepoints.

IMPORTANT DESIGN LIMITATION
The automated method was developed on 190320 and frozen before
application to 190321.

However, the Andersson S10 manual measurements and the automated
validation derive from the same 190321 biological acquisition/date.
Therefore this is a held-out METHOD validation, not an independent
biological replication.

The 20 manually scored cells are not identified one-to-one with
the 301 automatically scored cells. Population trajectories and
distributions are compared; no paired-cell claim is made.

RAW MEDIAN TRAJECTORY
Pearson r = {raw['pearson_r']:.6f}
Spearman rho = {raw['spearman_rho']:.6f}
MAE = {raw['MAE']:.6f}
RMSE = {raw['RMSE']:.6f}
mean automated-manual bias = {raw['mean_bias_automated_minus_manual']:.6f}

10-min manual median = {raw['manual_start']:.6f}
10-min automated median = {raw['automated_start']:.6f}

120-min manual median = {raw['manual_end']:.6f}
120-min automated median = {raw['automated_end']:.6f}

NORMALIZED WITHIN-CELL TRAJECTORY
Pearson r = {norm['pearson_r']:.6f}
Spearman rho = {norm['spearman_rho']:.6f}
MAE = {norm['MAE']:.6f}
RMSE = {norm['RMSE']:.6f}
mean automated-manual bias = {norm['mean_bias_automated_minus_manual']:.6f}

Manual normalized time-averaged AUC = {norm['manual_time_averaged_AUC']:.6f}
Automated normalized time-averaged AUC = {norm['automated_time_averaged_AUC']:.6f}

ENDPOINT DISTRIBUTIONS
Manual n = {int(manual_ep['n_cells'])}
Manual median T120/T10 = {manual_ep['median_T120_over_T10']:.6f}
Manual fraction decreased = {manual_ep['fraction_decreased']:.6f}

Automated n = {int(auto_ep['n_cells'])}
Automated median T120/T10 = {auto_ep['median_T120_over_T10']:.6f}
Automated fraction decreased = {auto_ep['fraction_decreased']:.6f}

Endpoint median difference
automated minus manual =
{auto_ep['median_T120_over_T10'] - manual_ep['median_T120_over_T10']:.6f}

STATISTICAL INTERPRETATION
Pearson and Spearman p-values in TRAJECTORY_CONCORDANCE.tsv are
descriptive only. The 12 serial timepoints are not independent
replicates and their p-values must not be interpreted as formal
inferential evidence.

Likewise, no cell-level hypothesis test is performed between the
manual and automated endpoint distributions because the cohorts
are not paired and individual cells are nested within the same
biological acquisition.

NO PARAMETERS WERE RETUNED AFTER HELD-OUT APPLICATION.
"""


(
    OUT
    / "FORMAL_VALIDATION_RESULTS.txt"
).write_text(
    text
)


                                                              
                         
                                                              

hash_targets = [
    OUT
    / "TIME_ALIGNED_COMPARISON.tsv",

    OUT
    / "TRAJECTORY_CONCORDANCE.tsv",

    OUT
    / "PER_SCENE_MANUAL_CONCORDANCE.tsv",

    OUT
    / "ENDPOINT_DISTRIBUTION_SUMMARY.tsv",

    OUT
    / "FORMAL_VALIDATION_RESULTS.txt",

    ROOT
    / "scripts"
    / "67_compare_frozen_automated_vs_andersson_S10.py",
]


with (
    OUT
    / "FORMAL_VALIDATION_SHA256SUMS.txt"
).open(
    "w"
) as fh:

    for p in hash_targets:

        fh.write(
            f"{file_sha256(p)}  {p}\n"
        )


                                                              
         
                                                              

print()
print("=" * 120)
print("TIME-ALIGNED COMPARISON")
print("=" * 120)

print(
    aligned[
        [
            "time_min",
            "manual_median_ratio",
            "automated_median_ratio",
            "raw_automated_minus_manual",
            "manual_median_relative_T1",
            "automated_median_relative_T1",
            "normalized_automated_minus_manual",
        ]
    ].to_string(
        index=False
    )
)


print()
print("=" * 120)
print("TRAJECTORY CONCORDANCE")
print("=" * 120)

print(
    metrics.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("PER-SCENE SENSITIVITY")
print("=" * 120)

print(
    scene_metrics.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("ENDPOINT DISTRIBUTIONS")
print("=" * 120)

print(
    endpoint_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("RESULT RECORD")
print("=" * 120)

print(
    text
)


print(
    "Output:",
    OUT
)
