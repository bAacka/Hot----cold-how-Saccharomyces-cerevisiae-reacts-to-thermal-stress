#!/usr/bin/env python3

from pathlib import Path
import ast
import sys

import numpy as np
import pandas as pd

from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    gaussian_filter,
    convolve,
)
from scipy.optimize import linear_sum_assignment

from skimage.measure import regionprops
from skimage.morphology import disk
from skimage.registration import phase_cross_correlation


ROOT = Path(
    "/"
)

BRIDGE = ROOT / "44_phenotype_bridge"

STATE = (
    BRIDGE /
    "46_phenotype_state"
)

SOURCE = (
    BRIDGE /
    "46s2_track_focus_190418_WT_hsp104D_T0_T8.py"
)

TRK = (
    STATE /
    "190418_WT_hsp104D_tracked_focus_T0_T8"
)

CACHE = (
    BRIDGE /
    "projection_cache_190418_WT_hsp104D"
)

SEG = (
    STATE /
    "190418_WT_hsp104D_dual_segmentation" /
    "yeast_BF_cp3"
)

OUT = (
    TRK /
    "46t_final_hsp104_validation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

MEAS0 = (
    TRK /
    "44g_focus_cytosol_measurements.tsv.gz"
)

                                                              
                           
                                                              

SCENE_TO_GENOTYPE = {
    0: "hsp104D",
    3: "WT",
}

PRIMARY_RADIUS = 4

                                
EARLY = list(range(1, 9))

                                       
                                                         
LATE = list(range(8, 16))

                                    
MAX_DISTANCE_PX = 30.0
MIN_AREA_RATIO = 0.40
MAX_AREA_RATIO = 2.50


                                                              
                                          
                                                              

src_text = SOURCE.read_text()

required_literals = [
    "TIMES = list(range(9))",
    "RADII = [3, 4, 5]",
    "PRIMARY_RADIUS = 4",
    "MAX_DISTANCE_PX = 30.0",
    "MIN_AREA_RATIO = 0.40",
    "MAX_AREA_RATIO = 2.50",
    "N_PER_SCENE = 5",
]

for q in required_literals:
    if q not in src_text:
        raise RuntimeError(
            f"Frozen-source check failed: {q}"
        )

print("=" * 80)
print("FROZEN SOURCE VERIFIED")
print("=" * 80)

for q in required_literals:
    print(q)


                                                              
                                                    
 
                                                          
                            
                                                              

wanted = {
    "mask_path",
    "get_background",
    "register_to_t0",
    "prop_table",
    "track_scene",
    "matched_roi_measure",
}

tree = ast.parse(src_text)

func_nodes = [
    node
    for node in tree.body
    if (
        isinstance(node, ast.FunctionDef)
        and node.name in wanted
    )
]

found = {
    node.name
    for node in func_nodes
}

missing = wanted - found

if missing:
    raise RuntimeError(
        "Could not recover frozen functions: "
        + ", ".join(sorted(missing))
    )

module = ast.Module(
    body=func_nodes,
    type_ignores=[],
)

ast.fix_missing_locations(module)

ns = {
    "np": np,
    "pd": pd,
    "SEG": SEG,

    "TIMES": list(range(16)),

    "MAX_DISTANCE_PX":
        MAX_DISTANCE_PX,

    "MIN_AREA_RATIO":
        MIN_AREA_RATIO,

    "MAX_AREA_RATIO":
        MAX_AREA_RATIO,

    "binary_dilation":
        binary_dilation,

    "binary_erosion":
        binary_erosion,

    "gaussian_filter":
        gaussian_filter,

    "convolve":
        convolve,

    "linear_sum_assignment":
        linear_sum_assignment,

    "regionprops":
        regionprops,

    "disk":
        disk,

    "phase_cross_correlation":
        phase_cross_correlation,
}

exec(
    compile(
        module,
        str(SOURCE),
        "exec",
    ),
    ns,
)

mask_path = ns["mask_path"]
get_background = ns["get_background"]
track_scene = ns["track_scene"]
matched_roi_measure = ns["matched_roi_measure"]


                                                              
        
                                                              

gfp = np.load(
    CACHE /
    "egfp_max_projection.npy",
    mmap_mode="r",
)

dic = np.load(
    CACHE /
    "dic_midplane.npy",
    mmap_mode="r",
)

if gfp.shape != dic.shape:
    raise RuntimeError(
        f"GFP/DIC mismatch: {gfp.shape} vs {dic.shape}"
    )

if gfp.shape[:2] != (4, 16):
    raise RuntimeError(
        f"Expected 4 x 16 acquisition; got {gfp.shape}"
    )

frozen = pd.read_csv(
    MEAS0,
    sep="\t",
    compression="gzip",
)

print()
print("=" * 80)
print("INPUT")
print("=" * 80)
print("cache shape:", gfp.shape)
print("frozen rows:", len(frozen))


                                                              
                                                       
            
 
                                 
                                                              

sel = (
    frozen[
        frozen["scene"].isin(
            SCENE_TO_GENOTYPE
        )
        &
        frozen["radius"].eq(
            PRIMARY_RADIUS
        )
    ][
        [
            "scene",
            "baseline_label",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "scene",
            "baseline_label",
        ]
    )
)

counts = (
    sel.groupby("scene")
    .size()
    .to_dict()
)

if counts != {0: 5, 3: 5}:
    raise RuntimeError(
        "Expected exactly five frozen cells "
        f"in S0 and S3; got {counts}"
    )

print()
print("Frozen selected cells:")
print(
    sel.to_string(
        index=False
    )
)


                                                              
                                                    
 
            
                           
                     
                             
                       
                          
                   
                                                              

all_tracks = {}

for scene in sorted(
    SCENE_TO_GENOTYPE
):

    print()
    print("=" * 80)
    print(
        "EXTEND TRACKING:",
        SCENE_TO_GENOTYPE[scene],
        f"S{scene}",
    )
    print("=" * 80)

    tracks = track_scene(
        dic,
        "190418_WT_hsp104D",
        scene,
    )

    all_tracks[scene] = tracks

    n0 = len(tracks)

    n16 = sum(
        set(range(16)).issubset(
            tr.keys()
        )
        for tr in tracks.values()
    )

    print(
        "baseline tracks:",
        n0,
    )

    print(
        "complete T0:T15:",
        n16,
    )


                                                              
                                                              
                                                              

rows = []

for r in sel.itertuples(
    index=False
):

    scene = int(r.scene)
    base_lab = int(r.baseline_label)

    genotype = (
        SCENE_TO_GENOTYPE[
            scene
        ]
    )

    tracks = all_tracks[
        scene
    ]

    if base_lab not in tracks:
        raise RuntimeError(
            f"Frozen selected cell vanished "
            f"from track table: S{scene} "
            f"label={base_lab}"
        )

    tr = tracks[
        base_lab
    ]

    for t in range(16):

        if t not in tr:

            rows.append({
                "genotype":
                    genotype,

                "scene":
                    scene,

                "baseline_label":
                    base_lab,

                "time_index":
                    t,

                "current_label":
                    np.nan,

                "tracked":
                    False,

                "ratio":
                    np.nan,
            })

            continue

        mask = np.load(
            mask_path(
                "190418_WT_hsp104D",
                scene,
                t,
            )
        )

        img = np.asarray(
            gfp[
                scene,
                t,
            ]
        )

        bg = get_background(
            img,
            mask,
        )

        current_label = int(
            tr[t]["label"]
        )

        z = matched_roi_measure(
            img=img,
            mask=mask,
            label=current_label,
            background=bg,
            radius=PRIMARY_RADIUS,
        )

        if z is None:

            rows.append({
                "genotype":
                    genotype,

                "scene":
                    scene,

                "baseline_label":
                    base_lab,

                "time_index":
                    t,

                "current_label":
                    current_label,

                "tracked":
                    True,

                "ratio":
                    np.nan,
            })

            continue

        row = {
            "genotype":
                genotype,

            "scene":
                scene,

            "baseline_label":
                base_lab,

            "time_index":
                t,

            "current_label":
                current_label,

            "tracked":
                True,

            "tracking_distance":
                float(
                    tr[t]["distance"]
                ),

            "tracking_area_ratio":
                float(
                    tr[t]["area_ratio"]
                ),

            "registration_error":
                float(
                    tr[t][
                        "registration_error"
                    ]
                ),
        }

                                                          
                             
        for k, v in z.items():
            row[k] = v

        rows.append(
            row
        )


ext = pd.DataFrame(
    rows
)

EXT_OUT = (
    OUT /
    "fixed_cells_T0_T15_radius4.tsv.gz"
)

ext.to_csv(
    EXT_OUT,
    sep="\t",
    index=False,
    compression="gzip",
)


                                                              
              
                                                              

print()
print("=" * 80)
print("FIXED-CELL COMPLETENESS")
print("=" * 80)

comp = (
    ext.assign(
        valid=lambda q:
            q["ratio"].notna()
    )
    .groupby(
        [
            "genotype",
            "baseline_label",
        ]
    )["valid"]
    .sum()
    .reset_index(
        name="n_valid_T0_T15"
    )
)

print(
    comp.to_string(
        index=False
    )
)


                                                              
                                                             
 
                                                     
                                                
 
                                             
                                                              

summary_rows = []
trajectory_rows = []

for (
    genotype,
    scene,
    cell,
), g in ext.groupby(
    [
        "genotype",
        "scene",
        "baseline_label",
    ]
):

    g = g.sort_values(
        "time_index"
    )

                               
                 
                               

    early = (
        g[
            g["time_index"].isin(
                EARLY
            )
        ]
        .copy()
        .sort_values(
            "time_index"
        )
    )

    if (
        len(early) == 8
        and early["ratio"].notna().all()
    ):

        r1 = float(
            early.loc[
                early[
                    "time_index"
                ].eq(1),
                "ratio",
            ].iloc[0]
        )

        ey = np.log2(
            early[
                "ratio"
            ].to_numpy(
                dtype=float
            )
            /
            r1
        )

        et = early[
            "time_index"
        ].to_numpy(
            dtype=float
        )

        early_auc = float(
            np.trapezoid(
                ey,
                x=et,
            )
        )

        early_t8 = float(
            ey[-1]
        )

    else:
        early_auc = np.nan
        early_t8 = np.nan

                               
                           
     
                   
                                      
                               

    late = (
        g[
            g["time_index"].isin(
                LATE
            )
        ]
        .copy()
        .sort_values(
            "time_index"
        )
    )

    if (
        len(late) == 8
        and late["ratio"].notna().all()
    ):

        r8 = float(
            late.loc[
                late[
                    "time_index"
                ].eq(8),
                "ratio",
            ].iloc[0]
        )

        ly = np.log2(
            late[
                "ratio"
            ].to_numpy(
                dtype=float
            )
            /
            r8
        )

        lt = late[
            "time_index"
        ].to_numpy(
            dtype=float
        )

        late_auc = float(
            np.trapezoid(
                ly,
                x=lt,
            )
        )

        late_t15 = float(
            ly[-1]
        )

        ratio_t8 = float(
            late["ratio"].iloc[0]
        )

        ratio_t15 = float(
            late["ratio"].iloc[-1]
        )

    else:
        late_auc = np.nan
        late_t15 = np.nan
        ratio_t8 = np.nan
        ratio_t15 = np.nan

    summary_rows.append({
        "genotype":
            genotype,

        "scene":
            int(scene),

        "baseline_label":
            int(cell),

        "early_AUC_T1_T8":
            early_auc,

        "early_log2_T8_vs_T1":
            early_t8,

        "late_AUC_T8_T15":
            late_auc,

        "late_log2_T15_vs_T8":
            late_t15,

        "ratio_T8":
            ratio_t8,

        "ratio_T15":
            ratio_t15,
    })

                                                        
    gg = g.copy()

    r1_rows = gg[
        gg["time_index"].eq(1)
    ]

    r8_rows = gg[
        gg["time_index"].eq(8)
    ]

    r1 = (
        float(
            r1_rows["ratio"].iloc[0]
        )
        if (
            len(r1_rows)
            and pd.notna(
                r1_rows["ratio"].iloc[0]
            )
        )
        else np.nan
    )

    r8 = (
        float(
            r8_rows["ratio"].iloc[0]
        )
        if (
            len(r8_rows)
            and pd.notna(
                r8_rows["ratio"].iloc[0]
            )
        )
        else np.nan
    )

    for rr in gg.itertuples():

        ratio = float(
            rr.ratio
        ) if pd.notna(
            rr.ratio
        ) else np.nan

        trajectory_rows.append({
            "genotype":
                genotype,

            "scene":
                int(scene),

            "baseline_label":
                int(cell),

            "time_index":
                int(
                    rr.time_index
                ),

            "ratio":
                ratio,

            "log2_vs_T1":
                (
                    np.log2(
                        ratio / r1
                    )
                    if (
                        np.isfinite(ratio)
                        and np.isfinite(r1)
                        and ratio > 0
                        and r1 > 0
                    )
                    else np.nan
                ),

            "log2_vs_T8":
                (
                    np.log2(
                        ratio / r8
                    )
                    if (
                        np.isfinite(ratio)
                        and np.isfinite(r8)
                        and ratio > 0
                        and r8 > 0
                    )
                    else np.nan
                ),
        })


cells = pd.DataFrame(
    summary_rows
)

traj = pd.DataFrame(
    trajectory_rows
)

cells.to_csv(
    OUT /
    "fixed_cell_early_and_late_endpoints.tsv",
    sep="\t",
    index=False,
)

traj.to_csv(
    OUT /
    "fixed_cell_T0_T15_trajectory.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


                                                              
                                 
                                            
                                                              

group = (
    cells
    .groupby(
        "genotype"
    )
    .agg(
        n_cells=(
            "baseline_label",
            "size",
        ),

        mean_early_AUC=(
            "early_AUC_T1_T8",
            "mean",
        ),

        median_early_AUC=(
            "early_AUC_T1_T8",
            "median",
        ),

        mean_early_T8=(
            "early_log2_T8_vs_T1",
            "mean",
        ),

        mean_late_AUC=(
            "late_AUC_T8_T15",
            "mean",
        ),

        median_late_AUC=(
            "late_AUC_T8_T15",
            "median",
        ),

        mean_late_T15=(
            "late_log2_T15_vs_T8",
            "mean",
        ),

        median_late_T15=(
            "late_log2_T15_vs_T8",
            "median",
        ),

        mean_ratio_T8=(
            "ratio_T8",
            "mean",
        ),

        mean_ratio_T15=(
            "ratio_T15",
            "mean",
        ),
    )
    .reset_index()
)

group.to_csv(
    OUT /
    "resolved_genotype_early_late_summary.tsv",
    sep="\t",
    index=False,
)


mean_traj = (
    traj
    .groupby(
        [
            "genotype",
            "time_index",
        ]
    )
    .agg(
        n_cells=(
            "baseline_label",
            "nunique",
        ),

        mean_ratio=(
            "ratio",
            "mean",
        ),

        median_ratio=(
            "ratio",
            "median",
        ),

        mean_log2_vs_T1=(
            "log2_vs_T1",
            "mean",
        ),

        mean_log2_vs_T8=(
            "log2_vs_T8",
            "mean",
        ),
    )
    .reset_index()
)

mean_traj.to_csv(
    OUT /
    "resolved_genotype_T0_T15_trajectory.tsv",
    sep="\t",
    index=False,
)


                                                              
                                     
 
                   
                                                  
                                                              

try:

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    QC = OUT / "qc_selected_cells"

    QC.mkdir(
        exist_ok=True
    )

    qc_times = [
        0,
        8,
        12,
        15,
    ]

    for (
        genotype,
        scene,
        cell,
    ), g in ext.groupby(
        [
            "genotype",
            "scene",
            "baseline_label",
        ]
    ):

        fig, axes = plt.subplots(
            1,
            4,
            figsize=(
                14,
                4,
            ),
        )

        for ax, t in zip(
            axes,
            qc_times,
        ):

            rr = g[
                g["time_index"].eq(t)
            ]

            ax.set_title(
                f"T{t}"
            )

            ax.axis(
                "off"
            )

            if (
                len(rr) != 1
                or pd.isna(
                    rr.iloc[0]["current_label"]
                )
            ):
                ax.text(
                    0.5,
                    0.5,
                    "not tracked",
                    ha="center",
                    va="center",
                )
                continue

            row = rr.iloc[0]

            current_label = int(
                row["current_label"]
            )

            mask = np.load(
                mask_path(
                    "190418_WT_hsp104D",
                    int(scene),
                    t,
                )
            )

            cmask = (
                mask ==
                current_label
            )

            props = regionprops(
                cmask.astype(
                    np.uint8
                )
            )

            if not props:
                continue

            p = props[0]

            cy, cx = p.centroid

            half = 55

            y0 = max(
                0,
                int(cy) - half,
            )
            y1 = min(
                gfp.shape[2],
                int(cy) + half,
            )
            x0 = max(
                0,
                int(cx) - half,
            )
            x1 = min(
                gfp.shape[3],
                int(cx) + half,
            )

            crop = np.asarray(
                gfp[
                    int(scene),
                    t,
                    y0:y1,
                    x0:x1,
                ]
            )

            lo, hi = np.percentile(
                crop,
                [
                    1,
                    99.8,
                ],
            )

            ax.imshow(
                crop,
                cmap="gray",
                vmin=lo,
                vmax=hi,
            )

            cell_crop = cmask[
                y0:y1,
                x0:x1,
            ]

            ax.contour(
                cell_crop.astype(
                    float
                ),
                levels=[0.5],
                linewidths=0.8,
            )

            if (
                "focus_y" in row.index
                and
                pd.notna(
                    row["focus_y"]
                )
            ):

                fy = float(
                    row["focus_y"]
                ) - y0

                fx = float(
                    row["focus_x"]
                ) - x0

                ax.add_patch(
                    Circle(
                        (fx, fy),
                        PRIMARY_RADIUS,
                        fill=False,
                        linewidth=1.2,
                    )
                )

            if (
                "cytosol_y" in row.index
                and
                pd.notna(
                    row["cytosol_y"]
                )
            ):

                yy = float(
                    row["cytosol_y"]
                ) - y0

                xx = float(
                    row["cytosol_x"]
                ) - x0

                ax.add_patch(
                    Circle(
                        (xx, yy),
                        PRIMARY_RADIUS,
                        fill=False,
                        linewidth=1.0,
                        linestyle="--",
                    )
                )

            ratio = row.get(
                "ratio",
                np.nan,
            )

            ax.text(
                0.02,
                0.98,
                f"R={ratio:.2f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
            )

        fig.suptitle(
            f"{genotype}  S{scene}  "
            f"baseline label {cell}"
        )

        fig.tight_layout()

        fout = (
            QC /
            f"{genotype}_S{scene}_cell{cell}.png"
        )

        fig.savefig(
            fout,
            dpi=180,
            bbox_inches="tight",
        )

        plt.close(
            fig
        )

    qc_status = "QC panels written"

except Exception as e:

    qc_status = (
        "QC panel generation skipped: "
        f"{type(e).__name__}: {e}"
    )


                                                              
                   
                                                              

lookup = (
    group
    .set_index(
        "genotype"
    )
)

lines = []

lines.append(
    "HSP104 PHENOTYPE VALIDATION — FINAL"
)

lines.append(
    "=" * 72
)

lines.append(
    ""
)

lines.append(
    "Resolved acquisition mapping:"
)

lines.append(
    "  scene 0 = hsp104D"
)

lines.append(
    "  scene 3 = WT"
)

lines.append(
    "  scenes 1/2 excluded from genotype interpretation"
)

lines.append(
    ""
)

lines.append(
    "PRIMARY / PRE-EXISTING FROZEN PHENOTYPE:"
)

lines.append(
    "  T1:T8 AUC of log2(focus/cytosol ratio relative to T1)"
)

lines.append(
    "  Cell selection was performed at T0 only."
)

lines.append(
    "  Five cells per resolved scene."
)

lines.append(
    ""
)

if (
    "WT" in lookup.index
    and
    "hsp104D" in lookup.index
):

    wt_e = float(
        lookup.loc[
            "WT",
            "mean_early_AUC",
        ]
    )

    mut_e = float(
        lookup.loc[
            "hsp104D",
            "mean_early_AUC",
        ]
    )

    wt_l = float(
        lookup.loc[
            "WT",
            "mean_late_AUC",
        ]
    )

    mut_l = float(
        lookup.loc[
            "hsp104D",
            "mean_late_AUC",
        ]
    )

    wt_end = float(
        lookup.loc[
            "WT",
            "mean_late_T15",
        ]
    )

    mut_end = float(
        lookup.loc[
            "hsp104D",
            "mean_late_T15",
        ]
    )

    lines.append(
        f"  WT mean early AUC:       {wt_e:.6f}"
    )

    lines.append(
        f"  hsp104D mean early AUC:  {mut_e:.6f}"
    )

    lines.append(
        f"  WT - hsp104D:            "
        f"{wt_e-mut_e:.6f}"
    )

    lines.append(
        ""
    )

    lines.append(
        "SECONDARY / POST-PRIMARY LATE PHASE:"
    )

    lines.append(
        "  Same frozen cells; no reselection."
    )

    lines.append(
        "  Same tracker and radius-4 ROI."
    )

    lines.append(
        "  T8:T15 AUC normalized to each cell's T8."
    )

    lines.append(
        ""
    )

    lines.append(
        f"  WT mean late AUC:        {wt_l:.6f}"
    )

    lines.append(
        f"  hsp104D mean late AUC:   {mut_l:.6f}"
    )

    lines.append(
        f"  WT - hsp104D:            "
        f"{wt_l-mut_l:.6f}"
    )

    lines.append(
        ""
    )

    lines.append(
        f"  WT mean log2(T15/T8):       "
        f"{wt_end:.6f}"
    )

    lines.append(
        f"  hsp104D mean log2(T15/T8):  "
        f"{mut_end:.6f}"
    )

lines.append(
    ""
)

lines.append(
    "INFERENCE LIMIT:"
)

lines.append(
    "  One resolved scene per genotype."
)

lines.append(
    "  Cells are not independent biological replicates."
)

lines.append(
    "  No genotype P-value is claimed."
)

lines.append(
    ""
)

lines.append(
    qc_status
)

FINAL = (
    OUT /
    "FINAL_HSP104_VALIDATION.txt"
)

FINAL.write_text(
    "\n".join(lines)
    + "\n"
)


                                                              
       
                                                              

print()
print("=" * 80)
print("CELL ENDPOINTS")
print("=" * 80)

print(
    cells.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6f}",
    )
)

print()
print("=" * 80)
print("GENOTYPE SUMMARY")
print("=" * 80)

print(
    group.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6f}",
    )
)

print()
print("=" * 80)
print("FINAL FREEZE")
print("=" * 80)

print(
    FINAL.read_text()
)

print("OUTPUT:")
print(OUT)
