#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, linregress

ROOT = Path("/")
OUT = ROOT / "44_phenotype_bridge"
STATE = OUT / "46_phenotype_state"

DEST = STATE / "46m_replicate_comparison_T1_T8"
DEST.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "190319": STATE / "190319_tracked_focus_T0_T8",
    "190326": STATE / "190326_tracked_focus_T0_T8",
}

PRIMARY_RADIUS = 4
PRIMARY_TIMES = list(range(1, 9))


def read_required(path, cols):
    if not path.is_file():
        raise FileNotFoundError(path)

    x = pd.read_csv(path, sep="\t")

    missing = set(cols) - set(x.columns)

    if missing:
        raise RuntimeError(
            f"{path}\nmissing columns: {sorted(missing)}\n"
            f"observed: {list(x.columns)}"
        )

    return x


                                                              
                                        
                                                              

traj_rows = []

for acquisition, base in DATASETS.items():

    x = read_required(
        base / "44g_trajectory.tsv",
        [
            "radius",
            "time_index",
            "n",
            "mean_ratio",
            "median_ratio",
            "sd_ratio",
        ],
    )

    x = x[
        x["radius"].eq(PRIMARY_RADIUS)
        &
        x["time_index"].isin(PRIMARY_TIMES)
    ].copy()

    x = x.sort_values("time_index")

    if x["time_index"].tolist() != PRIMARY_TIMES:
        raise RuntimeError(
            f"{acquisition}: incomplete T1-T8 trajectory"
        )

    r1 = float(
        x.loc[
            x["time_index"].eq(1),
            "mean_ratio",
        ].iloc[0]
    )

    x["acquisition"] = acquisition
    x["ratio_relative_T1"] = x["mean_ratio"] / r1
    x["log2_ratio_relative_T1"] = np.log2(
        x["ratio_relative_T1"]
    )

    traj_rows.append(x)

traj = pd.concat(
    traj_rows,
    ignore_index=True,
)

traj.to_csv(
    DEST / "46m_acquisition_trajectories_T1_T8.tsv",
    sep="\t",
    index=False,
)


                                                              
                               
                                                              

wide_raw = traj.pivot(
    index="time_index",
    columns="acquisition",
    values="mean_ratio",
)

wide_log = traj.pivot(
    index="time_index",
    columns="acquisition",
    values="log2_ratio_relative_T1",
)

a = wide_raw["190319"].to_numpy(float)
b = wide_raw["190326"].to_numpy(float)

la = wide_log["190319"].to_numpy(float)
lb = wide_log["190326"].to_numpy(float)

rp_raw, pp_raw = pearsonr(a, b)
rs_raw, ps_raw = spearmanr(a, b)

rp_log, pp_log = pearsonr(la, lb)
rs_log, ps_log = spearmanr(la, lb)

agreement = pd.DataFrame([
    {
        "comparison": "raw_mean_ratio_T1_T8",
        "pearson_r": rp_raw,
        "pearson_p": pp_raw,
        "spearman_rho": rs_raw,
        "spearman_p": ps_raw,
        "rmse": float(
            np.sqrt(
                np.mean(
                    (a - b) ** 2
                )
            )
        ),
    },
    {
        "comparison": "log2_relative_to_T1_T1_T8",
        "pearson_r": rp_log,
        "pearson_p": pp_log,
        "spearman_rho": rs_log,
        "spearman_p": ps_log,
        "rmse": float(
            np.sqrt(
                np.mean(
                    (la - lb) ** 2
                )
            )
        ),
    },
])

agreement.to_csv(
    DEST / "46m_replicate_agreement.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                    
 
                   
                                            
                                                              

candidate_rows = []

for acquisition, g in traj.groupby("acquisition"):

    g = g.sort_values("time_index")

    t = g["time_index"].to_numpy(float)
    raw = g["mean_ratio"].to_numpy(float)
    logrel = g["log2_ratio_relative_T1"].to_numpy(float)

    fit = linregress(
        t,
        logrel,
    )

    t1 = float(
        g.loc[
            g["time_index"].eq(1),
            "mean_ratio",
        ].iloc[0]
    )

    t3 = float(
        g.loc[
            g["time_index"].eq(3),
            "mean_ratio",
        ].iloc[0]
    )

    t4 = float(
        g.loc[
            g["time_index"].eq(4),
            "mean_ratio",
        ].iloc[0]
    )

    t8 = float(
        g.loc[
            g["time_index"].eq(8),
            "mean_ratio",
        ].iloc[0]
    )

    candidate_rows.append({
        "acquisition": acquisition,

        "T1_mean_ratio":
        t1,

        "T8_mean_ratio":
        t8,

        "T8_over_T1":
        t8 / t1,

        "log2_T8_over_T1":
        np.log2(t8 / t1),

        "auc_log2_relative_T1_T1_T8":
        float(
            np.trapezoid(
                logrel,
                t,
            )
        ),

        "linear_slope_log2_relative_T1":
        float(fit.slope),

        "linear_slope_r_squared":
        float(fit.rvalue ** 2),

        "early_log2_T3_over_T1":
        float(
            np.log2(
                t3 / t1
            )
        ),

        "mid_log2_T4_over_T1":
        float(
            np.log2(
                t4 / t1
            )
        ),

        "late_log2_T8_over_T4":
        float(
            np.log2(
                t8 / t4
            )
        ),

        "minimum_log2_relative_T1":
        float(
            np.min(logrel)
        ),

        "minimum_time_index":
        int(
            g.iloc[
                int(
                    np.argmin(logrel)
                )
            ]["time_index"]
        ),
    })

candidates = pd.DataFrame(
    candidate_rows
)

candidates.to_csv(
    DEST / "46m_candidate_acquisition_phenotypes.tsv",
    sep="\t",
    index=False,
)


                                                              
                            
 
                                                             
                                                                
                                                     
                                                              

scene_candidate_rows = []
cell_candidate_rows = []

for acquisition, base in DATASETS.items():

    m = read_required(
        base / "44g_focus_cytosol_measurements.tsv.gz",
        [
            "scene",
            "baseline_label",
            "time_index",
            "radius",
            "ratio",
        ],
    )

    m = m[
        m["radius"].eq(PRIMARY_RADIUS)
        &
        m["time_index"].isin(PRIMARY_TIMES)
    ].copy()

                                                              
                                       
                                                           
                                                              

    for (
        scene,
        baseline_label
    ), g in m.groupby([
        "scene",
        "baseline_label",
    ]):

        g = g.sort_values("time_index")

        if g["time_index"].tolist() != PRIMARY_TIMES:
            raise RuntimeError(
                f"{acquisition} scene={scene} "
                f"cell={baseline_label}: incomplete T1-T8"
            )

        t = g["time_index"].to_numpy(float)
        raw = g["ratio"].to_numpy(float)

        r1 = raw[0]

        if not np.isfinite(r1) or r1 <= 0:
            raise RuntimeError(
                f"{acquisition} scene={scene} "
                f"cell={baseline_label}: invalid T1 ratio"
            )

        logrel = np.log2(
            raw / r1
        )

        fit = linregress(
            t,
            logrel,
        )

        cell_candidate_rows.append({
            "acquisition":
            acquisition,

            "scene":
            int(scene),

            "baseline_label":
            int(baseline_label),

            "T1_ratio":
            float(raw[0]),

            "T8_ratio":
            float(raw[-1]),

            "T8_over_T1":
            float(
                raw[-1] /
                raw[0]
            ),

            "log2_T8_over_T1":
            float(logrel[-1]),

            "auc_log2_relative_T1_T1_T8":
            float(
                np.trapezoid(
                    logrel,
                    t,
                )
            ),

            "linear_slope_log2_relative_T1":
            float(fit.slope),

            "minimum_log2_relative_T1":
            float(
                np.min(logrel)
            ),

            "minimum_time_index":
            int(
                t[
                    np.argmin(logrel)
                ]
            ),
        })

                                                              
                                                           
                                                              

    st = (
        m.groupby(
            [
                "scene",
                "time_index",
            ],
            as_index=False,
        )
        ["ratio"]
        .mean()
    )

    for scene, g in st.groupby("scene"):

        g = g.sort_values(
            "time_index"
        )

        if g["time_index"].tolist() != PRIMARY_TIMES:
            raise RuntimeError(
                f"{acquisition} S{scene}: "
                "incomplete scene T1-T8"
            )

        t = g["time_index"].to_numpy(float)
        raw = g["ratio"].to_numpy(float)

        logrel = np.log2(
            raw / raw[0]
        )

        fit = linregress(
            t,
            logrel,
        )

        scene_candidate_rows.append({
            "acquisition":
            acquisition,

            "scene":
            int(scene),

            "T1_mean_ratio":
            float(raw[0]),

            "T8_mean_ratio":
            float(raw[-1]),

            "T8_over_T1":
            float(
                raw[-1] /
                raw[0]
            ),

            "log2_T8_over_T1":
            float(logrel[-1]),

            "auc_log2_relative_T1_T1_T8":
            float(
                np.trapezoid(
                    logrel,
                    t,
                )
            ),

            "linear_slope_log2_relative_T1":
            float(fit.slope),

            "linear_slope_r_squared":
            float(fit.rvalue ** 2),

            "minimum_log2_relative_T1":
            float(
                np.min(logrel)
            ),

            "minimum_time_index":
            int(
                t[
                    np.argmin(logrel)
                ]
            ),
        })


cells = pd.DataFrame(
    cell_candidate_rows
)

scenes = pd.DataFrame(
    scene_candidate_rows
)

cells.to_csv(
    DEST / "46m_candidate_cell_phenotypes.tsv",
    sep="\t",
    index=False,
)

scenes.to_csv(
    DEST / "46m_candidate_scene_phenotypes.tsv",
    sep="\t",
    index=False,
)


                                                              
                                       
 
                                                          
                                                        
                                                              

metrics = [
    "T8_over_T1",
    "log2_T8_over_T1",
    "auc_log2_relative_T1_T1_T8",
    "linear_slope_log2_relative_T1",
]

summary_rows = []

for level, df in [
    ("cell", cells),
    ("scene", scenes),
]:

    for acquisition, g in df.groupby(
        "acquisition"
    ):

        for metric in metrics:

            z = pd.to_numeric(
                g[metric],
                errors="coerce",
            ).dropna()

            summary_rows.append({
                "level":
                level,

                "acquisition":
                acquisition,

                "metric":
                metric,

                "n":
                len(z),

                "mean":
                float(z.mean()),

                "median":
                float(z.median()),

                "q25":
                float(z.quantile(0.25)),

                "q75":
                float(z.quantile(0.75)),

                "sd":
                float(z.std(ddof=1))
                if len(z) > 1
                else np.nan,
            })

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    DEST / "46m_candidate_metric_distribution_summary.tsv",
    sep="\t",
    index=False,
)


                                                              
                         
                                                              

print()
print("=" * 100)
print("T1-T8 REPLICATE AGREEMENT")
print("=" * 100)

print(
    agreement.to_string(
        index=False,
        float_format=lambda z: f"{z:.6g}",
    )
)

print()
print("=" * 100)
print("CANDIDATE ACQUISITION-LEVEL PHENOTYPES")
print("DESCRIPTIVE ONLY — NOT FROZEN")
print("=" * 100)

print(
    candidates.to_string(
        index=False,
        float_format=lambda z: f"{z:.6g}",
    )
)

print()
print("=" * 100)
print("SCENE-LEVEL CANDIDATE PHENOTYPES")
print("=" * 100)

print(
    scenes.to_string(
        index=False,
        float_format=lambda z: f"{z:.6g}",
    )
)

print()
print("=" * 100)
print("CELL-LEVEL DISTRIBUTION SUMMARY")
print("CELLS ARE NESTED — NO INDEPENDENCE CLAIM")
print("=" * 100)

print(
    summary[
        summary["level"].eq("cell")
    ].to_string(
        index=False,
        float_format=lambda z: f"{z:.6g}",
    )
)

(
    DEST /
    "46m_STATUS.txt"
).write_text(
    "Replicate comparison complete.\n"
    "Primary analysis window: T1-T8.\n"
    "Primary ROI radius: 4 px.\n"
    "T0 excluded from phenotype construction.\n"
    "Candidate scalar phenotypes are descriptive only and not yet frozen.\n"
    "Cells are nested within scenes/acquisitions and are not treated as independent biological replicates.\n"
)

print()
print("OUTPUT:", DEST)
