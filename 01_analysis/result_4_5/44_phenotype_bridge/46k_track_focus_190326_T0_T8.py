#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.ndimage import (
    binary_dilation,
    binary_erosion,
    gaussian_filter,
    convolve,
)
from scipy.optimize import linear_sum_assignment
from scipy.stats import linregress

from skimage.measure import regionprops
from skimage.morphology import disk
from skimage.registration import phase_cross_correlation


ROOT = Path("/")
OUT = ROOT / "44_phenotype_bridge"

SEG = (
    OUT /
    "46_phenotype_state" /
    "190326_dual_segmentation" /
    "yeast_BF_cp3"
)

DEST = OUT / "46_phenotype_state" / "190326_tracked_focus_T0_T8"
DEST.mkdir(parents=True, exist_ok=True)

CACHES = {
    "190326_ssa1-2DD":
        OUT / "projection_cache_190326_ssa1-2DD",
}

                                                   
TIMES = list(range(9))

                            
                                   
RADII = [3, 4, 5]

PRIMARY_RADIUS = 4

                                                                    
MAX_DISTANCE_PX = 30.0
MIN_AREA_RATIO = 0.40
MAX_AREA_RATIO = 2.50

                                              
                                       
                                                      
N_PER_SCENE = 5


def mask_path(genotype, scene, time_index):
    return (
        SEG /
        genotype /
        "masks" /
        f"{genotype}_S{scene}_T{time_index:02d}.npy"
    )


def get_background(img, mask):

    occupied = mask > 0

    blocked = binary_dilation(
        occupied,
        iterations=5,
    )

    pixels = img[~blocked]

    if pixels.size < 1000:
        pixels = img[~occupied]

    if pixels.size < 100:
        raise RuntimeError(
            "Too few extracellular background pixels"
        )

    return float(np.median(pixels))


def register_to_t0(dic0, dict_):

                                                                 
    a = gaussian_filter(
        dic0.astype(np.float32),
        sigma=1.0,
    )[::4, ::4]

    b = gaussian_filter(
        dict_.astype(np.float32),
        sigma=1.0,
    )[::4, ::4]

    shift, error, _ = phase_cross_correlation(
        a,
        b,
        upsample_factor=4,
    )

                                                        
                                                             
    shift = np.asarray(
        shift,
        dtype=float,
    ) * 4.0

    return shift, float(error)


def prop_table(mask):

    rows = []

    for p in regionprops(mask):

        rows.append({
            "label": int(p.label),
            "y": float(p.centroid[0]),
            "x": float(p.centroid[1]),
            "area": float(p.area),
        })

    return pd.DataFrame(rows)


def track_scene(dic, genotype, scene):
    """
    Sequential frame-to-frame tracking.

    Previous version matched every timepoint directly to T0.
    That is unnecessarily brittle across ~110 min of cell
    movement/budding. Here, identity propagates one frame at
    a time using:
      - DIC translational registration
      - centroid distance
      - area consistency
      - one-to-one Hungarian assignment

    A track is still considered complete only if it is
    successfully matched at every T0...T11 frame.
    """

    m0 = np.load(
        mask_path(
            genotype,
            scene,
            0,
        )
    )

    base = prop_table(m0)

    if len(base) == 0:
        raise RuntimeError(
            f"No baseline cells: {genotype} S{scene}"
        )

    tracks = {}

                                                            
                                         
    active = {}

    for r in base.itertuples():

        lab = int(r.label)

        tracks[lab] = {
            0: {
                "label": lab,
                "distance": 0.0,
                "area_ratio": 1.0,
                "registration_error": 0.0,
            }
        }

        active[lab] = {
            "label": lab,
            "y": float(r.y),
            "x": float(r.x),
            "area": float(r.area),
        }

    prev_dic = np.asarray(
        dic[scene, 0],
        dtype=np.float32,
    )

    for t in TIMES[1:]:

        mt = np.load(
            mask_path(
                genotype,
                scene,
                t,
            )
        )

        cur = prop_table(mt)

        if len(cur) == 0 or len(active) == 0:
            active = {}
            break

        curr_dic = np.asarray(
            dic[scene, t],
            dtype=np.float32,
        )

                                                          
                    
        shift, reg_error = register_to_t0(
            prev_dic,
            curr_dic,
        )

        active_ids = list(active.keys())

        prev_xy = np.array(
            [
                [
                    active[k]["y"],
                    active[k]["x"],
                ]
                for k in active_ids
            ],
            dtype=float,
        )

        prev_area = np.array(
            [
                active[k]["area"]
                for k in active_ids
            ],
            dtype=float,
        )

        cur_xy = cur[
            ["y", "x"]
        ].to_numpy(dtype=float)

        cur_area = cur[
            "area"
        ].to_numpy(dtype=float)

                                                         
                                                          
                                                           
                       
        cur_xy_aligned = (
            cur_xy +
            shift[None, :]
        )

        dy = (
            prev_xy[:, None, 0] -
            cur_xy_aligned[None, :, 0]
        )

        dx = (
            prev_xy[:, None, 1] -
            cur_xy_aligned[None, :, 1]
        )

        dist = np.sqrt(
            dy * dy +
            dx * dx
        )

        area_ratio = (
            cur_area[None, :] /
            prev_area[:, None]
        )

                                                      
                                           
        area_penalty = (
            8.0 *
            np.abs(
                np.log(
                    np.clip(
                        area_ratio,
                        1e-6,
                        None,
                    )
                )
            )
        )

        cost = (
            dist +
            area_penalty
        )

        ii, jj = linear_sum_assignment(
            cost
        )

        next_active = {}

        for i, j in zip(ii, jj):

            d = float(
                dist[i, j]
            )

            ar = float(
                area_ratio[i, j]
            )

            if d > MAX_DISTANCE_PX:
                continue

            if not (
                MIN_AREA_RATIO
                <= ar <=
                MAX_AREA_RATIO
            ):
                continue

            baseline_label = int(
                active_ids[i]
            )

            current_label = int(
                cur.iloc[j]["label"]
            )

            tracks[
                baseline_label
            ][t] = {
                "label":
                current_label,

                "distance":
                d,

                "area_ratio":
                ar,

                "registration_error":
                reg_error,
            }

                                                   
                                                          
                                        
            next_active[
                baseline_label
            ] = {
                "label":
                current_label,

                "y":
                float(
                    cur.iloc[j]["y"]
                ),

                "x":
                float(
                    cur.iloc[j]["x"]
                ),

                "area":
                float(
                    cur.iloc[j]["area"]
                ),
            }

        active = next_active
        prev_dic = curr_dic

    return tracks


def matched_roi_measure(
    img,
    mask,
    label,
    background,
    radius,
):

    cmask = (
        mask ==
        int(label)
    )

    if cmask.sum() < 100:
        return None

    kernel = disk(
        radius
    ).astype(bool)

                                                       
    centers = binary_erosion(
        cmask,
        structure=kernel,
        border_value=0,
    )

    if centers.sum() < 10:
        return None

    corr = (
        img.astype(np.float64) -
        background
    )

    smooth = gaussian_filter(
        corr,
        sigma=1.0,
    )

    coords = np.argwhere(
        centers
    )

    vals = smooth[
        centers
    ]

    k = int(
        np.argmax(vals)
    )

    fy, fx = map(
        int,
        coords[k],
    )

                          
    kfloat = kernel.astype(
        np.float64
    )

    local_sum = convolve(
        corr,
        kfloat,
        mode="constant",
        cval=0.0,
    )

    roi_area = int(
        kernel.sum()
    )

    focus_sum = float(
        local_sum[fy, fx]
    )

                                                     
    candidates = centers.copy()

    yy, xx = np.ogrid[
        :mask.shape[0],
        :mask.shape[1]
    ]

    exclude = (
        (yy - fy) ** 2 +
        (xx - fx) ** 2
        <=
        (3 * radius) ** 2
    )

    candidates[exclude] = False

    if candidates.sum() < 5:

                                              
        candidates = centers.copy()

        exclude = (
            (yy - fy) ** 2 +
            (xx - fx) ** 2
            <=
            (2 * radius) ** 2
        )

        candidates[exclude] = False

    if candidates.sum() < 5:
        return None

                             
                                                   
                                      
    cell_median = float(
        np.median(
            corr[cmask]
        )
    )

    ccoords = np.argwhere(
        candidates
    )

    candidate_means = (
        local_sum[candidates] /
        roi_area
    )

    j = int(
        np.argmin(
            np.abs(
                candidate_means -
                cell_median
            )
        )
    )

    cy, cx = map(
        int,
        ccoords[j],
    )

    cyt_sum = float(
        local_sum[cy, cx]
    )

    if not np.isfinite(
        focus_sum
    ):
        return None

    if not np.isfinite(
        cyt_sum
    ):
        return None

    if cyt_sum <= 0:
        return None

    ratio = (
        focus_sum /
        cyt_sum
    )

    return {
        "ratio": float(ratio),

        "focus_sum_ctcf":
        focus_sum,

        "cytosol_sum_ctcf":
        cyt_sum,

        "roi_area_px":
        roi_area,

        "focus_y":
        fy,

        "focus_x":
        fx,

        "cytosol_y":
        cy,

        "cytosol_x":
        cx,

        "cell_median_corr":
        cell_median,
    }


all_measurements = []
tracking_rows = []
selection_rows = []


for genotype, cache in CACHES.items():

    print()
    print("=" * 100)
    print("GENOTYPE:", genotype)
    print("=" * 100)

    gfp = np.load(
        cache /
        "egfp_max_projection.npy",
        mmap_mode="r",
    )

    dic = np.load(
        cache /
        "dic_midplane.npy",
        mmap_mode="r",
    )

    genotype_candidates = []

    scene_tracks = {}

                                                              
                          
                                                              

    for scene in range(4):

        tracks = track_scene(
            dic,
            genotype,
            scene,
        )

        scene_tracks[scene] = tracks

        complete = [
            lab
            for lab, tr in tracks.items()
            if all(
                t in tr
                for t in TIMES
            )
        ]

        print(
            f"S{scene}: "
            f"baseline={len(tracks)} "
            f"complete_tracks={len(complete)}"
        )

        for lab, tr in tracks.items():

            tracking_rows.append({
                "genotype":
                genotype,

                "scene":
                scene,

                "baseline_label":
                lab,

                "n_timepoints":
                len(tr),

                "complete":
                len(tr) == len(TIMES),

                "max_tracking_distance":
                max(
                    z["distance"]
                    for z in tr.values()
                ),

                "min_area_ratio":
                min(
                    z["area_ratio"]
                    for z in tr.values()
                ),

                "max_area_ratio":
                max(
                    z["area_ratio"]
                    for z in tr.values()
                ),
            })

                                                              
                                                           
                                                              
                                                  
                                                              

        mask0 = np.load(
            mask_path(
                genotype,
                scene,
                0,
            )
        )

        img0 = np.asarray(
            gfp[scene, 0],
            dtype=np.float32,
        )

        bg0 = get_background(
            img0,
            mask0,
        )

        candidates = []

        for lab in complete:

            m = matched_roi_measure(
                img0,
                mask0,
                lab,
                bg0,
                PRIMARY_RADIUS,
            )

            if m is None:
                continue

            if not np.isfinite(
                m["ratio"]
            ):
                continue

            candidates.append(
                (
                    float(m["ratio"]),
                    int(lab),
                )
            )

        candidates.sort(
            reverse=True
        )

        chosen = candidates[
            :N_PER_SCENE
        ]

        if len(chosen) < N_PER_SCENE:
            raise RuntimeError(
                f"{genotype} S{scene}: "
                f"only {len(chosen)} usable "
                f"complete baseline-focus tracks"
            )

        print(
            "  selected baseline ratios:",
            ", ".join(
                f"{r:.3f}"
                for r, _ in chosen
            ),
        )

        for rank, (
            baseline_ratio,
            lab
        ) in enumerate(
            chosen,
            start=1,
        ):

            genotype_candidates.append(
                (
                    scene,
                    lab,
                    baseline_ratio,
                )
            )

            selection_rows.append({
                "genotype":
                genotype,

                "scene":
                scene,

                "baseline_label":
                lab,

                "scene_rank":
                rank,

                "baseline_ratio_R4":
                baseline_ratio,
            })

    if len(genotype_candidates) != 20:
        raise RuntimeError(
            f"{genotype}: expected exactly "
            f"20 selected tracks, got "
            f"{len(genotype_candidates)}"
        )

                                                              
                                                        
                                                              

    for scene, base_lab, base_ratio in genotype_candidates:

        tr = scene_tracks[
            scene
        ][
            base_lab
        ]

        for t in TIMES:

            current_label = int(
                tr[t]["label"]
            )

            mask = np.load(
                mask_path(
                    genotype,
                    scene,
                    t,
                )
            )

            img = np.asarray(
                gfp[scene, t],
                dtype=np.float32,
            )

            bg = get_background(
                img,
                mask,
            )

            for radius in RADII:

                m = matched_roi_measure(
                    img,
                    mask,
                    current_label,
                    bg,
                    radius,
                )

                if m is None:
                    continue

                all_measurements.append({
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

                    "radius":
                    radius,

                    "frame_background":
                    bg,

                    "tracking_distance":
                    tr[t]["distance"],

                    "tracking_area_ratio":
                    tr[t]["area_ratio"],

                    **m,
                })


tracking = pd.DataFrame(
    tracking_rows
)

selection = pd.DataFrame(
    selection_rows
)

measurements = pd.DataFrame(
    all_measurements
)

tracking.to_csv(
    DEST /
    "44g_tracking_qc.tsv",
    sep="\t",
    index=False,
)

selection.to_csv(
    DEST /
    "44g_selected_tracks.tsv",
    sep="\t",
    index=False,
)

measurements.to_csv(
    DEST /
    "44g_focus_cytosol_measurements.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


                                                              
                    
                                                              

expected_rows = (
    len(CACHES) *
    4 *
    N_PER_SCENE *
    len(TIMES) *
    len(RADII)
)

print()
print("=" * 100)
print("MEASUREMENT COMPLETENESS")
print("=" * 100)

print(
    "expected:",
    expected_rows
)

print(
    "observed:",
    len(measurements)
)

if len(measurements) != expected_rows:

    print(
        "WARNING: some ROI measurements "
        "were not obtainable."
    )


                                                              
                         
                                                              

trajectory = (
    measurements
    .groupby(
        [
            "genotype",
            "radius",
            "time_index",
        ],
        as_index=False,
    )
    .agg(
        n=("ratio", "size"),
        mean_ratio=("ratio", "mean"),
        median_ratio=("ratio", "median"),
        sd_ratio=("ratio", "std"),
    )
)

trajectory.to_csv(
    DEST /
    "44g_trajectory.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                         
                                 
                                                              

reg_rows = []

for (
    genotype,
    radius
), g in trajectory.groupby(
    ["genotype", "radius"]
):

    g = g.sort_values(
        "time_index"
    )

    fit = linregress(
        g["time_index"],
        g["mean_ratio"],
    )

    reg_rows.append({
        "genotype":
        genotype,

        "radius":
        radius,

        "slope_per_frame":
        float(fit.slope),

        "intercept":
        float(fit.intercept),

        "r_squared":
        float(
            fit.rvalue ** 2
        ),

        "p_linear_slope":
        float(fit.pvalue),

        "ratio_T0":
        float(
            g.loc[
                g["time_index"] == 0,
                "mean_ratio",
            ].iloc[0]
        ),

        "ratio_T8":
        float(
            g.loc[
                g["time_index"] == TIMES[-1],
                "mean_ratio",
            ].iloc[0]
        ),

        "delta_T8_T0":
        float(
            g.loc[
                g["time_index"] == TIMES[-1],
                "mean_ratio",
            ].iloc[0]
            -
            g.loc[
                g["time_index"] == 0,
                "mean_ratio",
            ].iloc[0]
        ),
    })

reg = pd.DataFrame(
    reg_rows
)

reg.to_csv(
    DEST /
    "44g_regression_summary.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                 
                                                              

scene_rows = []

scene_traj = (
    measurements
    .groupby(
        [
            "genotype",
            "scene",
            "radius",
            "time_index",
        ],
        as_index=False,
    )
    ["ratio"]
    .mean()
)

for (
    genotype,
    scene,
    radius
), g in scene_traj.groupby(
    [
        "genotype",
        "scene",
        "radius",
    ]
):

    if g["time_index"].nunique() < len(TIMES):
        continue

    fit = linregress(
        g["time_index"],
        g["ratio"],
    )

    scene_rows.append({
        "genotype":
        genotype,

        "scene":
        scene,

        "radius":
        radius,

        "slope_per_frame":
        float(fit.slope),

        "r_squared":
        float(
            fit.rvalue ** 2
        ),
    })

scene_slopes = pd.DataFrame(
    scene_rows
)

scene_slopes.to_csv(
    DEST /
    "44g_scene_slopes.tsv",
    sep="\t",
    index=False,
)


                                                              
                       
                                                              

print()
print("=" * 100)
print("PRIMARY TRACKED FOCUS/CYTOSOL RESULT")
print("=" * 100)

print(
    reg.to_string(
        index=False,
        float_format=lambda z:
        f"{z:.6g}",
    )
)

print()
print("=" * 100)
print("PRIMARY RADIUS = 4 TRAJECTORY")
print("=" * 100)

r4 = trajectory[
    trajectory["radius"] ==
    PRIMARY_RADIUS
]

wide = r4.pivot(
    index="time_index",
    columns="genotype",
    values="mean_ratio",
)

print(
    wide.to_string(
        float_format=lambda z:
        f"{z:.4f}"
    )
)

print()
print("=" * 100)
print("RADIUS=4 SCENE SLOPES")
print("=" * 100)

print(
    scene_slopes[
        scene_slopes["radius"] ==
        PRIMARY_RADIUS
    ].to_string(
        index=False,
        float_format=lambda z:
        f"{z:.6g}",
    )
)

(DEST / "44g_COMPLETE.txt").write_text(
    "Tracked matched-ROI focus/cytosol analysis complete.\n"
)
