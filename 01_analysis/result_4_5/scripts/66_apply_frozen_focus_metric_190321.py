#!/usr/bin/env python3

from pathlib import Path
import hashlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import (
    binary_dilation,
    convolve,
    distance_transform_edt,
    gaussian_filter,
)

from skimage.morphology import disk
from aicspylibczi import CziFile


ROOT = Path(
    "/"
)

BASE = (
    ROOT
    / "44_phenotype_bridge"
)

DATA = (
    BASE
    / "downloads"
    / "190321_WT"
)

MASK_ROOT = (
    BASE
    / "segmentation"
    / "190321_WT"
    / "all_scenes_halfscale"
)

TRACK_ROOT = (
    BASE
    / "tracking"
    / "190321_WT"
)

FREEZE_FILE = (
    BASE
    / "focus_metric_freeze"
    / "FOCUS_METRIC_V3_PREVALIDATION_FREEZE.txt"
)

V3_SOURCE = (
    ROOT
    / "scripts"
    / "64_transport_focus_metric_190320_v3.py"
)

OUT = (
    BASE
    / "focus_metric_validation"
    / "190321_WT"
    / "frozen_v3"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                                                          
                                          
                                                              

EXPECTED_FREEZE_SHA256 = (
    "2a12ddb34b6521f35acf147c89276ca0dc2815a6f3eeb8ed640c4101da4a8415"
)

EXPECTED_V3_SOURCE_SHA256 = (
    "d15675c82c93076a42e14b763a253e724fa2d54e68b0d7ef54fe14897854e6d7"
)


def sha256(path):

    return hashlib.sha256(
        path.read_bytes()
    ).hexdigest()


observed_freeze_sha = sha256(
    FREEZE_FILE
)

observed_v3_sha = sha256(
    V3_SOURCE
)


if observed_freeze_sha != EXPECTED_FREEZE_SHA256:
    raise RuntimeError(
        "Frozen phenotype specification has changed.\n"
        f"Expected: {EXPECTED_FREEZE_SHA256}\n"
        f"Observed: {observed_freeze_sha}"
    )


if observed_v3_sha != EXPECTED_V3_SOURCE_SHA256:
    raise RuntimeError(
        "Frozen V3 development source has changed.\n"
        f"Expected: {EXPECTED_V3_SOURCE_SHA256}\n"
        f"Observed: {observed_v3_sha}"
    )


print("Frozen specification SHA256: VERIFIED")
print("Frozen V3 source SHA256:    VERIFIED")


                                                              
                   
 
                                                  
                                
                                                              

SCENES = [
    0,
    1,
    2,
    3,
]

                           
                          
TIMES = list(
    range(
        1,
        13,
    )
)

ZS = list(
    range(
        10
    )
)

GFP_CHANNEL = 1

BACKGROUND_DILATION_PX = 5

ROI_RADIUS_PX = 4

SMALL_SIGMA_PX = 1.0
BROAD_SIGMA_PX = 6.0

MIN_FOCUS_CYT_DISTANCE_PX = (
    2 * ROI_RADIUS_PX
    + 2
)

FOCUS_POSITIVE_Z_THRESHOLD = 5.0

EPS = 1e-9


                                                              
     
                                                              

files = sorted(
    DATA.glob(
        "*.czi"
    )
)

if len(files) != 1:
    raise RuntimeError(
        f"Expected exactly one held-out WT CZI; found {len(files)}"
    )

czi_path = files[0]

print()
print("Held-out CZI:")
print(czi_path)

czi = CziFile(
    czi_path
)


                                                              
                
                                                              

kernel = disk(
    ROI_RADIUS_PX
).astype(
    np.float32
)

kernel_area = float(
    kernel.sum()
)

if int(
    kernel_area
) != 49:
    raise RuntimeError(
        f"Unexpected ROI area: {kernel_area}"
    )


                                                              
         
                                                              

def read_gfp(
    scene,
    t,
    z,
):

    img, _ = czi.read_image(
        S=scene,
        T=t,
        C=GFP_CHANNEL,
        Z=z,
    )

    x = np.squeeze(
        np.asarray(
            img
        )
    )

    if x.ndim != 2:
        raise RuntimeError(
            f"Scene {scene} T={t} Z={z}: "
            f"expected YX; got {x.shape}"
        )

    return x.astype(
        np.float32,
        copy=False,
    )


                                                              
                            
 
                                     
                                                              
 
              
                                    
                                                              

def make_projection(
    scene,
    mask,
    t,
):

    stack = np.stack(
        [
            read_gfp(
                scene=scene,
                t=t,
                z=z,
            )
            for z in ZS
        ],
        axis=0,
    )

    raw_proj = np.max(
        stack,
        axis=0,
    )


    occupied = (
        mask > 0
    )

    protected = binary_dilation(
        occupied,
        iterations=BACKGROUND_DILATION_PX,
    )

    bg_mask = (
        ~protected
    )

    bg_values = raw_proj[
        bg_mask
    ]

    if bg_values.size == 0:
        raise RuntimeError(
            f"Scene {scene} T={t}: no extracellular background"
        )


    background = float(
        np.median(
            bg_values
        )
    )


    corrected = (
        raw_proj
        - background
    ).astype(
        np.float32,
        copy=False,
    )


    return (
        corrected,
        background,
    )


                                                              
             
                                                              

def make_maps(
    projection,
):

    raw_disk_mean = (
        convolve(
            projection,
            kernel,
            mode="constant",
            cval=0.0,
        )
        / kernel_area
    )


    fine = gaussian_filter(
        projection,
        sigma=SMALL_SIGMA_PX,
    )

    broad = gaussian_filter(
        projection,
        sigma=BROAD_SIGMA_PX,
    )

    localized = (
        fine
        - broad
    )


    local_disk_mean = (
        convolve(
            localized,
            kernel,
            mode="constant",
            cval=0.0,
        )
        / kernel_area
    )


    return (
        raw_disk_mean,
        local_disk_mean,
    )


                                                              
                                            
                                                              

def candidate_centres(
    mask,
    label,
):

    cell = (
        mask == label
    )

    if not np.any(
        cell
    ):
        return None


    dist = distance_transform_edt(
        cell
    )

    valid = (
        dist
        >= (
            ROI_RADIUS_PX
            + 0.5
        )
    )

    yy, xx = np.nonzero(
        valid
    )

    if len(
        yy
    ) < 2:
        return None


    return (
        yy.astype(
            int
        ),
        xx.astype(
            int
        ),
    )


                                                              
                                
                                                              

def robust_prominence(
    local_values,
    selected_value,
):

    med = float(
        np.median(
            local_values
        )
    )

    mad = float(
        np.median(
            np.abs(
                local_values
                - med
            )
        )
    )

    sigma = (
        1.4826
        * mad
    )

    if sigma <= EPS:
        return np.nan


    return (
        float(
            selected_value
        )
        - med
    ) / sigma


                                                              
                                   
                                                              

eligible_rows = []
init_rows = []
long_rows = []
failure_rows = []
cohort_rows = []


                                                              
                                   
                                                              

for scene in SCENES:

    print()
    print("=" * 110)
    print(
        f"SCENE {scene}"
    )
    print("=" * 110)


    mask_dir = (
        MASK_ROOT
        / f"scene{scene}"
    )

    track_file = (
        TRACK_ROOT
        / f"scene{scene}"
        / "CELL_TRACKS_LONG.tsv"
    )


    tracks = pd.read_csv(
        track_file,
        sep="\t",
    )


                                                              
                                    
     
                                                              
                                                             
     
                                      
                                     
                                                              

    track_summary = (
        tracks.groupby(
            "track_id"
        )
        .agg(
            first_time=(
                "time_index",
                "min",
            ),
            last_time=(
                "time_index",
                "max",
            ),
            n_frames_total=(
                "time_index",
                "nunique",
            ),
        )
        .reset_index()
    )


    tsets = (
        tracks.groupby(
            "track_id"
        )[
            "time_index"
        ]
        .apply(
            lambda x:
                set(
                    int(v)
                    for v in x
                )
        )
    )


    required_times = set(
        TIMES
    )


    eligible_ids = []

    for _, r in track_summary.iterrows():

        track_id = int(
            r[
                "track_id"
            ]
        )

        if int(
            r[
                "first_time"
            ]
        ) != 0:
            continue

        if not required_times.issubset(
            tsets.loc[
                track_id
            ]
        ):
            continue

        eligible_ids.append(
            track_id
        )


    eligible_ids = sorted(
        eligible_ids
    )


    print(
        "T0-origin tracks complete through T1-T12:",
        len(
            eligible_ids
        ),
    )


    track_lookup = (
        tracks.set_index(
            [
                "track_id",
                "time_index",
            ]
        )
    )


                                                              
                       
                                                              

    t = 1

    mask = np.load(
        mask_dir
        / f"T{t:02d}_masks.npy"
    )

    projection, background = make_projection(
        scene=scene,
        mask=mask,
        t=t,
    )

    raw_map, local_map = make_maps(
        projection
    )


                                                              
                               
                                                              

    states = {}


    for track_id in eligible_ids:

        tr = track_lookup.loc[
            (
                track_id,
                1,
            )
        ]

        label = int(
            tr[
                "label"
            ]
        )


        cand = candidate_centres(
            mask,
            label,
        )

        if cand is None:

            failure_rows.append({
                "scene":
                    scene,

                "track_id":
                    track_id,

                "time_index":
                    1,

                "failure":
                    "no_valid_T1_ROI_centres",
            })

            continue


        yy, xx = cand

        local_values = local_map[
            yy,
            xx,
        ]

        raw_values = raw_map[
            yy,
            xx,
        ]


                                           
        fi = int(
            np.argmax(
                local_values
            )
        )

        fy = int(
            yy[
                fi
            ]
        )

        fx = int(
            xx[
                fi
            ]
        )


        initial_z = robust_prominence(
            local_values,
            local_map[
                fy,
                fx,
            ],
        )


                                          
        sep = np.hypot(
            yy - fy,
            xx - fx,
        )

        cyt_ok = (
            sep
            >= MIN_FOCUS_CYT_DISTANCE_PX
        )


        if not np.any(
            cyt_ok
        ):

            failure_rows.append({
                "scene":
                    scene,

                "track_id":
                    track_id,

                "time_index":
                    1,

                "failure":
                    "no_valid_T1_cytosol_ROI",
            })

            continue


        eligible_cyt = np.where(
            cyt_ok
        )[0]


        target_cyt_mean = float(
            np.median(
                raw_values[
                    cyt_ok
                ]
            )
        )


        ci = int(
            eligible_cyt[
                np.argmin(
                    np.abs(
                        raw_values[
                            cyt_ok
                        ]
                        - target_cyt_mean
                    )
                )
            ]
        )


        cy = int(
            yy[
                ci
            ]
        )

        cx = int(
            xx[
                ci
            ]
        )


        focus_mean = float(
            raw_map[
                fy,
                fx,
            ]
        )

        cyt_mean = float(
            raw_map[
                cy,
                cx,
            ]
        )


        ratio = (
            focus_mean
            / cyt_mean
            if cyt_mean > EPS
            else np.nan
        )


        focus_positive = bool(
            np.isfinite(
                initial_z
            )
            and
            (
                initial_z
                >= FOCUS_POSITIVE_Z_THRESHOLD
            )
        )


        eligible_rows.append({
            "scene":
                scene,

            "track_id":
                track_id,

            "label_T1":
                label,

            "initial_focus_prominence_z":
                initial_z,

            "initial_focus_to_cytosol_ratio":
                ratio,

            "focus_positive_z_ge_5":
                focus_positive,
        })


        if not focus_positive:
            continue


        init_rows.append({
            "scene":
                scene,

            "track_id":
                track_id,

            "label_T1":
                label,

            "focus_y_T1":
                fy,

            "focus_x_T1":
                fx,

            "cytosol_y_T1":
                cy,

            "cytosol_x_T1":
                cx,

            "initial_focus_prominence_z":
                initial_z,

            "initial_focus_to_cytosol_ratio":
                ratio,

            "focus_mean_T1":
                focus_mean,

            "cytosol_mean_T1":
                cyt_mean,
        })


        states[
            track_id
        ] = {
            "focus_y":
                float(
                    fy
                ),

            "focus_x":
                float(
                    fx
                ),

            "cytosol_y":
                float(
                    cy
                ),

            "cytosol_x":
                float(
                    cx
                ),

            "prev_centroid_y":
                float(
                    tr[
                        "centroid_y"
                    ]
                ),

            "prev_centroid_x":
                float(
                    tr[
                        "centroid_x"
                    ]
                ),

            "initial_z":
                float(
                    initial_z
                ),

            "initial_ratio":
                float(
                    ratio
                ),
        }


    print(
        "Frozen z>=5 focus-positive tracks:",
        len(
            states
        ),
    )


    cohort_rows.append({
        "scene":
            scene,

        "n_T0_origin_complete_T1_T12":
            len(
                eligible_ids
            ),

        "n_T1_geometry_valid":
            sum(
                1
                for r in eligible_rows
                if r[
                    "scene"
                ] == scene
            ),

        "n_focus_positive_z_ge_5":
            len(
                states
            ),
    })


                                                              
                                               
     
                    
                                                         
                                                              

    active_states = dict(
        states
    )


    for t in TIMES:

        print(
            f"Scene {scene} T={t:02d} "
            f"active={len(active_states)}"
        )


        mask = np.load(
            mask_dir
            / f"T{t:02d}_masks.npy"
        )

        projection, background = make_projection(
            scene=scene,
            mask=mask,
            t=t,
        )

        raw_map, local_map = make_maps(
            projection
        )


        failed_now = []


        for track_id, state in active_states.items():

            tr = track_lookup.loc[
                (
                    track_id,
                    t,
                )
            ]

            label = int(
                tr[
                    "label"
                ]
            )

            centroid_y = float(
                tr[
                    "centroid_y"
                ]
            )

            centroid_x = float(
                tr[
                    "centroid_x"
                ]
            )


            cand = candidate_centres(
                mask,
                label,
            )

            if cand is None:

                failure_rows.append({
                    "scene":
                        scene,

                    "track_id":
                        track_id,

                    "time_index":
                        t,

                    "failure":
                        "no_valid_longitudinal_ROI_centres",
                })

                failed_now.append(
                    track_id
                )

                continue


            yy, xx = cand


                                                              
                                      
                                                              

            if t == 1:

                pred_fy = state[
                    "focus_y"
                ]

                pred_fx = state[
                    "focus_x"
                ]

                pred_cy = state[
                    "cytosol_y"
                ]

                pred_cx = state[
                    "cytosol_x"
                ]


            else:

                dy = (
                    centroid_y
                    - state[
                        "prev_centroid_y"
                    ]
                )

                dx = (
                    centroid_x
                    - state[
                        "prev_centroid_x"
                    ]
                )


                pred_fy = (
                    state[
                        "focus_y"
                    ]
                    + dy
                )

                pred_fx = (
                    state[
                        "focus_x"
                    ]
                    + dx
                )

                pred_cy = (
                    state[
                        "cytosol_y"
                    ]
                    + dy
                )

                pred_cx = (
                    state[
                        "cytosol_x"
                    ]
                    + dx
                )


                                                              
                             
                                                      
                                                              

            df = np.hypot(
                yy - pred_fy,
                xx - pred_fx,
            )

            fi = int(
                np.argmin(
                    df
                )
            )

            fy = int(
                yy[
                    fi
                ]
            )

            fx = int(
                xx[
                    fi
                ]
            )

            focus_snap_error = float(
                df[
                    fi
                ]
            )


                                                              
                      
                                                   
                                         
                                                              

            sep = np.hypot(
                yy - fy,
                xx - fx,
            )

            cyt_ok = (
                sep
                >= MIN_FOCUS_CYT_DISTANCE_PX
            )


            if np.any(
                cyt_ok
            ):

                eligible_cyt = np.where(
                    cyt_ok
                )[0]

                dc = np.hypot(
                    yy[
                        eligible_cyt
                    ]
                    - pred_cy,
                    xx[
                        eligible_cyt
                    ]
                    - pred_cx,
                )

                ci = int(
                    eligible_cyt[
                        np.argmin(
                            dc
                        )
                    ]
                )

            else:

                                                   
                ci = int(
                    np.argmax(
                        sep
                    )
                )


            cy = int(
                yy[
                    ci
                ]
            )

            cx = int(
                xx[
                    ci
                ]
            )

            cyt_snap_error = float(
                np.hypot(
                    cy - pred_cy,
                    cx - pred_cx,
                )
            )


                                                              
                            
                                                              

            focus_mean = float(
                raw_map[
                    fy,
                    fx,
                ]
            )

            cyt_mean = float(
                raw_map[
                    cy,
                    cx,
                ]
            )


            ratio = (
                focus_mean
                / cyt_mean
                if cyt_mean > EPS
                else np.nan
            )


            difference = (
                focus_mean
                - cyt_mean
            )


            denom = (
                focus_mean
                + cyt_mean
            )

            bounded_contrast = (
                difference
                / denom
                if abs(
                    denom
                ) > EPS
                else np.nan
            )


            current_z = robust_prominence(
                local_map[
                    yy,
                    xx,
                ],
                local_map[
                    fy,
                    fx,
                ],
            )


            long_rows.append({
                "scene":
                    scene,

                "track_id":
                    track_id,

                "time_index":
                    t,

                "biological_min":
                    int(
                        t * 10
                    ),

                "label":
                    label,

                "centroid_y":
                    centroid_y,

                "centroid_x":
                    centroid_x,

                "focus_y":
                    fy,

                "focus_x":
                    fx,

                "cytosol_y":
                    cy,

                "cytosol_x":
                    cx,

                "focus_snap_error_px":
                    focus_snap_error,

                "cytosol_snap_error_px":
                    cyt_snap_error,

                "focus_mean_bg_corrected":
                    focus_mean,

                "cytosol_mean_bg_corrected":
                    cyt_mean,

                "focus_to_cytosol_ratio":
                    ratio,

                "focus_minus_cytosol":
                    difference,

                "bounded_focus_contrast":
                    bounded_contrast,

                "current_focus_prominence_z":
                    current_z,

                "initial_focus_prominence_z":
                    state[
                        "initial_z"
                    ],

                "initial_focus_to_cytosol_ratio":
                    state[
                        "initial_ratio"
                    ],

                "projection_background":
                    background,
            })


            state[
                "focus_y"
            ] = float(
                fy
            )

            state[
                "focus_x"
            ] = float(
                fx
            )

            state[
                "cytosol_y"
            ] = float(
                cy
            )

            state[
                "cytosol_x"
            ] = float(
                cx
            )

            state[
                "prev_centroid_y"
            ] = centroid_y

            state[
                "prev_centroid_x"
            ] = centroid_x


                                                          
                                                               
        for track_id in failed_now:
            active_states.pop(
                track_id,
                None,
            )


                                                              
                                                                
                                                              

eligible = pd.DataFrame(
    eligible_rows
)

initial = pd.DataFrame(
    init_rows
)

cohort = pd.DataFrame(
    cohort_rows
)

failures = pd.DataFrame(
    failure_rows
)


eligible.to_csv(
    OUT
    / "HELDOUT_T1_ELIGIBLE_INITIALIZATION.tsv",
    sep="\t",
    index=False,
)

initial.to_csv(
    OUT
    / "HELDOUT_T1_FOCUS_POSITIVE.tsv",
    sep="\t",
    index=False,
)

cohort.to_csv(
    OUT
    / "HELDOUT_COHORT_COUNTS.tsv",
    sep="\t",
    index=False,
)

failures.to_csv(
    OUT
    / "HELDOUT_GEOMETRY_FAILURES.tsv",
    sep="\t",
    index=False,
)


                                                              
            
                                                              

long = pd.DataFrame(
    long_rows
)


                                                           
long[
    "cell_id"
] = (
    "S"
    + long[
        "scene"
    ].astype(
        str
    )
    + "_T"
    + long[
        "track_id"
    ].astype(
        str
    )
)


                   
t1 = (
    long[
        long[
            "time_index"
        ] == 1
    ][
        [
            "cell_id",
            "focus_to_cytosol_ratio",
        ]
    ]
    .rename(
        columns={
            "focus_to_cytosol_ratio":
                "T1_ratio"
        }
    )
)


long = long.merge(
    t1,
    on="cell_id",
    how="left",
    validate="m:1",
)


long[
    "ratio_relative_to_T1"
] = (
    long[
        "focus_to_cytosol_ratio"
    ]
    / long[
        "T1_ratio"
    ]
)


long.to_csv(
    OUT
    / "HELDOUT_FOCUS_LONG.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                              
 
                              
                                 
                                                              

complete_metric = (
    long.groupby(
        "cell_id"
    )[
        "time_index"
    ]
    .nunique()
)

complete_ids = set(
    complete_metric[
        complete_metric
        == len(
            TIMES
        )
    ].index
)


long_complete = long[
    long[
        "cell_id"
    ].isin(
        complete_ids
    )
].copy()


long_complete.to_csv(
    OUT
    / "HELDOUT_FOCUS_LONG_COMPLETE_T1_T12.tsv",
    sep="\t",
    index=False,
)


                                                              
                    
 
                                                      
                               
                                                              

pooled = (
    long_complete.groupby(
        [
            "time_index",
            "biological_min",
        ]
    )
    .agg(
        n_cells=(
            "cell_id",
            "nunique",
        ),

        n_scenes=(
            "scene",
            "nunique",
        ),

        median_ratio=(
            "focus_to_cytosol_ratio",
            "median",
        ),

        mean_ratio=(
            "focus_to_cytosol_ratio",
            "mean",
        ),

        q25_ratio=(
            "focus_to_cytosol_ratio",
            lambda x:
                x.quantile(
                    0.25
                ),
        ),

        q75_ratio=(
            "focus_to_cytosol_ratio",
            lambda x:
                x.quantile(
                    0.75
                ),
        ),

        median_ratio_relative_T1=(
            "ratio_relative_to_T1",
            "median",
        ),

        median_focus_minus_cytosol=(
            "focus_minus_cytosol",
            "median",
        ),

        median_current_z=(
            "current_focus_prominence_z",
            "median",
        ),

        median_focus_snap_error_px=(
            "focus_snap_error_px",
            "median",
        ),

        q95_focus_snap_error_px=(
            "focus_snap_error_px",
            lambda x:
                x.quantile(
                    0.95
                ),
        ),
    )
    .reset_index()
)


pooled.to_csv(
    OUT
    / "HELDOUT_POOLED_TRAJECTORY.tsv",
    sep="\t",
    index=False,
)


                                                              
                       
                                                              

per_scene = (
    long_complete.groupby(
        [
            "scene",
            "time_index",
            "biological_min",
        ]
    )
    .agg(
        n_cells=(
            "cell_id",
            "nunique",
        ),

        median_ratio=(
            "focus_to_cytosol_ratio",
            "median",
        ),

        median_ratio_relative_T1=(
            "ratio_relative_to_T1",
            "median",
        ),

        median_focus_minus_cytosol=(
            "focus_minus_cytosol",
            "median",
        ),

        median_current_z=(
            "current_focus_prominence_z",
            "median",
        ),
    )
    .reset_index()
)


per_scene.to_csv(
    OUT
    / "HELDOUT_PER_SCENE_TRAJECTORY.tsv",
    sep="\t",
    index=False,
)


                                                              
                    
                                                              

endpoint = (
    long_complete[
        long_complete[
            "time_index"
        ].isin(
            [
                1,
                12,
            ]
        )
    ][
        [
            "cell_id",
            "scene",
            "track_id",
            "time_index",
            "focus_to_cytosol_ratio",
            "ratio_relative_to_T1",
            "current_focus_prominence_z",
        ]
    ]
    .pivot(
        index=[
            "cell_id",
            "scene",
            "track_id",
        ],
        columns="time_index",
    )
)

endpoint.columns = [
    f"{a}_T{b}"
    for a, b in endpoint.columns
]

endpoint = (
    endpoint.reset_index()
)


endpoint[
    "T12_lower_than_T1"
] = (
    endpoint[
        "focus_to_cytosol_ratio_T12"
    ]
    < endpoint[
        "focus_to_cytosol_ratio_T1"
    ]
)


endpoint.to_csv(
    OUT
    / "HELDOUT_ENDPOINT_BY_CELL.tsv",
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
    pooled[
        "biological_min"
    ],
    pooled[
        "median_ratio"
    ],
    marker="o",
)

ax.fill_between(
    pooled[
        "biological_min"
    ],
    pooled[
        "q25_ratio"
    ],
    pooled[
        "q75_ratio"
    ],
    alpha=0.2,
)

ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xlabel(
    "Recovery time (min)"
)

ax.set_ylabel(
    "Original-focus / transported cytosol ratio"
)

ax.set_title(
    "190321 WT — frozen automated held-out phenotype"
)

fig.tight_layout()

fig.savefig(
    OUT
    / "HELDOUT_AUTOMATED_TRAJECTORY.png",
    dpi=180,
)

plt.close(
    fig
)


fig, ax = plt.subplots(
    figsize=(
        8,
        5,
    )
)

for scene, g in per_scene.groupby(
    "scene"
):

    ax.plot(
        g[
            "biological_min"
        ],
        g[
            "median_ratio"
        ],
        marker="o",
        label=f"scene {scene}",
    )

ax.axhline(
    1.0,
    linestyle="--",
)

ax.set_xlabel(
    "Recovery time (min)"
)

ax.set_ylabel(
    "Median original-focus / cytosol ratio"
)

ax.set_title(
    "190321 WT — held-out scene consistency"
)

ax.legend()

fig.tight_layout()

fig.savefig(
    OUT
    / "HELDOUT_PER_SCENE_TRAJECTORIES.png",
    dpi=180,
)

plt.close(
    fig
)


                                                              
                      
                                                              

manifest = f"""FROZEN V3 HELD-OUT APPLICATION

Held-out dataset:
190321 WT

Formal validation window:
T1-T12 = 10-120 min

Frozen focus-positive threshold:
initial T1 prominence z >= {FOCUS_POSITIVE_Z_THRESHOLD:.1f}

ROI radius:
{ROI_RADIUS_PX} px

ROI area:
{int(kernel_area)} px

Background dilation:
{BACKGROUND_DILATION_PX} px

Localized signal:
Gaussian sigma {SMALL_SIGMA_PX} minus sigma {BROAD_SIGMA_PX}

Post-T1 position rule:
geometry-only transport using tracked cell-centroid displacement.
No fluorescence-guided reacquisition.

Frozen specification SHA256:
{observed_freeze_sha}

Frozen V3 development source SHA256:
{observed_v3_sha}

IMPORTANT:
This script does not read the Andersson S10 manual reference.
All files in this directory were generated before automated/manual
trajectory comparison.
"""

(
    OUT
    / "HELDOUT_APPLICATION_MANIFEST.txt"
).write_text(
    manifest
)


                                                              
                 
                                                              

print()
print("=" * 120)
print("HELD-OUT COHORT COUNTS")
print("=" * 120)

print(
    cohort.to_string(
        index=False
    )
)


print()
print(
    "Total focus-positive initialized cells:",
    initial.shape[
        0
    ],
)

print(
    "Complete automated T1-T12 cells:",
    len(
        complete_ids
    ),
)


print()
print("=" * 120)
print("POOLED HELD-OUT AUTOMATED TRAJECTORY")
print("=" * 120)

print(
    pooled.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("PER-SCENE AUTOMATED TRAJECTORY")
print("=" * 120)

print(
    per_scene.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("ENDPOINT")
print("=" * 120)

print(
    "Cells with T12 < T1:",
    int(
        endpoint[
            "T12_lower_than_T1"
        ].sum()
    ),
    "/",
    len(
        endpoint
    ),
)

print(
    "Fraction:",
    float(
        endpoint[
            "T12_lower_than_T1"
        ].mean()
    ),
)

print(
    "Median individual T12/T1:",
    float(
        endpoint[
            "ratio_relative_to_T1_T12"
        ].median()
    ),
)


                                                              
                                                           
                                                              

hash_targets = [
    OUT
    / "HELDOUT_T1_ELIGIBLE_INITIALIZATION.tsv",

    OUT
    / "HELDOUT_T1_FOCUS_POSITIVE.tsv",

    OUT
    / "HELDOUT_COHORT_COUNTS.tsv",

    OUT
    / "HELDOUT_GEOMETRY_FAILURES.tsv",

    OUT
    / "HELDOUT_FOCUS_LONG.tsv",

    OUT
    / "HELDOUT_FOCUS_LONG_COMPLETE_T1_T12.tsv",

    OUT
    / "HELDOUT_POOLED_TRAJECTORY.tsv",

    OUT
    / "HELDOUT_PER_SCENE_TRAJECTORY.tsv",

    OUT
    / "HELDOUT_ENDPOINT_BY_CELL.tsv",

    OUT
    / "HELDOUT_APPLICATION_MANIFEST.txt",

    ROOT
    / "scripts"
    / "66_apply_frozen_focus_metric_190321.py",
]


with (
    OUT
    / "HELDOUT_AUTOMATED_PRECOMPARISON_SHA256SUMS.txt"
).open(
    "w"
) as fh:

    for p in hash_targets:

        h = hashlib.sha256(
            p.read_bytes()
        ).hexdigest()

        fh.write(
            f"{h}  {p}\n"
        )


print()
print(
    "Output:",
    OUT
)

print(
    "Automated held-out outputs hashed BEFORE S10 comparison."
)
