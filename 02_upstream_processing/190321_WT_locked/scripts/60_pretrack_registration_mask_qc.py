#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.ndimage import shift as ndi_shift


ROOT = Path(
    "/"
)

MASK_ROOT = (
    ROOT
    / "44_phenotype_bridge"
    / "segmentation"
    / "190321_WT"
    / "all_scenes_halfscale"
)

REG_FILE = (
    ROOT
    / "44_phenotype_bridge"
    / "qc"
    / "190321_WT"
    / "registration_focus"
    / "REGISTRATION_FOCUS_QC.tsv"
)

OUT = (
    ROOT
    / "44_phenotype_bridge"
    / "qc"
    / "190321_WT"
    / "pretrack_registration"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


reg = pd.read_csv(
    REG_FILE,
    sep="\t",
)

SCENES = sorted(
    reg["scene"].unique()
)


                                                              
                                
                                                              

chain_rows = []

for scene in SCENES:

    g = (
        reg[
            reg["scene"] == scene
        ]
        .sort_values(
            "time_index"
        )
        .copy()
    )

    prev_direct_dy = None
    prev_direct_dx = None

    for _, r in g.iterrows():

        t = int(
            r["time_index"]
        )

        cum_dy = float(
            r["cumulative_shift_dy_px"]
        )

        cum_dx = float(
            r["cumulative_shift_dx_px"]
        )

        direct_dy = float(
            r["direct_T0_shift_dy_px"]
        )

        direct_dx = float(
            r["direct_T0_shift_dx_px"]
        )


        chain_resid = np.hypot(
            cum_dy - direct_dy,
            cum_dx - direct_dx,
        )


        if prev_direct_dy is None:

            direct_delta_dy = 0.0
            direct_delta_dx = 0.0
            pair_delta_resid = 0.0

        else:

            direct_delta_dy = (
                direct_dy
                - prev_direct_dy
            )

            direct_delta_dx = (
                direct_dx
                - prev_direct_dx
            )

            pair_delta_resid = np.hypot(
                float(
                    r["pair_shift_dy_px"]
                )
                - direct_delta_dy,

                float(
                    r["pair_shift_dx_px"]
                )
                - direct_delta_dx,
            )


        chain_rows.append({
            "scene":
                scene,

            "time_index":
                t,

            "pair_dy":
                float(
                    r["pair_shift_dy_px"]
                ),

            "pair_dx":
                float(
                    r["pair_shift_dx_px"]
                ),

            "direct_delta_dy":
                direct_delta_dy,

            "direct_delta_dx":
                direct_delta_dx,

            "pair_vs_direct_delta_residual_px":
                pair_delta_resid,

            "cumulative_dy":
                cum_dy,

            "cumulative_dx":
                cum_dx,

            "direct_T0_dy":
                direct_dy,

            "direct_T0_dx":
                direct_dx,

            "cumulative_vs_direct_residual_px":
                chain_resid,
        })


        prev_direct_dy = direct_dy
        prev_direct_dx = direct_dx


chain = pd.DataFrame(
    chain_rows
)

chain.to_csv(
    OUT
    / "REGISTRATION_CHAIN_CONSISTENCY.tsv",
    sep="\t",
    index=False,
)


                                                              
                     
                                                              

def best_iou_stats(
    prev,
    curr,
):

    prev = np.asarray(
        prev,
        dtype=np.int32,
    )

    curr = np.asarray(
        curr,
        dtype=np.int32,
    )

    max_prev = int(
        prev.max()
    )

    max_curr = int(
        curr.max()
    )


    prev_area = np.bincount(
        prev.ravel(),
        minlength=max_prev + 1,
    )

    curr_area = np.bincount(
        curr.ravel(),
        minlength=max_curr + 1,
    )


    best_prev = np.zeros(
        max_prev + 1,
        dtype=np.float64,
    )

    best_curr = np.zeros(
        max_curr + 1,
        dtype=np.float64,
    )


    valid = (
        (prev > 0)
        &
        (curr > 0)
    )


    if np.any(
        valid
    ):

        base = (
            max_curr
            + 1
        )

        codes = (
            prev[
                valid
            ].astype(
                np.int64
            )
            * base
            + curr[
                valid
            ].astype(
                np.int64
            )
        )


        unique_codes, intersections = np.unique(
            codes,
            return_counts=True,
        )


        pids = (
            unique_codes
            // base
        ).astype(
            int
        )

        cids = (
            unique_codes
            % base
        ).astype(
            int
        )


        unions = (
            prev_area[
                pids
            ]
            + curr_area[
                cids
            ]
            - intersections
        )


        ious = (
            intersections
            / unions
        )


        np.maximum.at(
            best_prev,
            pids,
            ious,
        )

        np.maximum.at(
            best_curr,
            cids,
            ious,
        )


    prev_values = best_prev[
        1:
    ]

    curr_values = best_curr[
        1:
    ]


    return {
        "n_prev":
            max_prev,

        "n_curr":
            max_curr,

        "prev_median_best_iou":
            float(
                np.median(
                    prev_values
                )
            ),

        "prev_q05_best_iou":
            float(
                np.quantile(
                    prev_values,
                    0.05,
                )
            ),

        "prev_frac_iou_ge_020":
            float(
                np.mean(
                    prev_values
                    >= 0.20
                )
            ),

        "prev_frac_iou_ge_050":
            float(
                np.mean(
                    prev_values
                    >= 0.50
                )
            ),

        "curr_median_best_iou":
            float(
                np.median(
                    curr_values
                )
            ),

        "curr_frac_iou_ge_020":
            float(
                np.mean(
                    curr_values
                    >= 0.20
                )
            ),
    }


                                                              
                                      
                                         
                                                              

alignment_rows = []


for scene in SCENES:

    scene_reg = (
        reg[
            reg["scene"] == scene
        ]
        .set_index(
            "time_index"
        )
    )


    times = sorted(
        int(x)
        for x in scene_reg.index
    )


    for t in times[1:]:

        prev_t = (
            t - 1
        )


        prev_path = (
            MASK_ROOT
            / f"scene{scene}"
            / f"T{prev_t:02d}_masks.npy"
        )

        curr_path = (
            MASK_ROOT
            / f"scene{scene}"
            / f"T{t:02d}_masks.npy"
        )


        prev = np.load(
            prev_path
        )

        curr = np.load(
            curr_path
        )


        dy = float(
            scene_reg.loc[
                t,
                "pair_shift_dy_px",
            ]
        )

        dx = float(
            scene_reg.loc[
                t,
                "pair_shift_dx_px",
            ]
        )


        curr_aligned = ndi_shift(
            curr,
            shift=(
                dy,
                dx,
            ),
            order=0,
            mode="constant",
            cval=0,
            prefilter=False,
        ).astype(
            curr.dtype
        )


        raw_stats = best_iou_stats(
            prev,
            curr,
        )

        aligned_stats = best_iou_stats(
            prev,
            curr_aligned,
        )


        alignment_rows.append({
            "scene":
                scene,

            "time_index":
                t,

            "pair_shift_dy_px":
                dy,

            "pair_shift_dx_px":
                dx,

            "pair_shift_magnitude_px":
                float(
                    np.hypot(
                        dy,
                        dx,
                    )
                ),

            "raw_prev_median_best_iou":
                raw_stats[
                    "prev_median_best_iou"
                ],

            "aligned_prev_median_best_iou":
                aligned_stats[
                    "prev_median_best_iou"
                ],

            "raw_prev_frac_iou_ge_020":
                raw_stats[
                    "prev_frac_iou_ge_020"
                ],

            "aligned_prev_frac_iou_ge_020":
                aligned_stats[
                    "prev_frac_iou_ge_020"
                ],

            "aligned_prev_frac_iou_ge_050":
                aligned_stats[
                    "prev_frac_iou_ge_050"
                ],

            "aligned_prev_q05_best_iou":
                aligned_stats[
                    "prev_q05_best_iou"
                ],

            "aligned_curr_median_best_iou":
                aligned_stats[
                    "curr_median_best_iou"
                ],

            "aligned_curr_frac_iou_ge_020":
                aligned_stats[
                    "curr_frac_iou_ge_020"
                ],
        })


alignment = pd.DataFrame(
    alignment_rows
)

alignment.to_csv(
    OUT
    / "MASK_ALIGNMENT_QC.tsv",
    sep="\t",
    index=False,
)


                                                              
               
                                                              

summary_rows = []


for scene in SCENES:

    a = alignment[
        alignment[
            "scene"
        ] == scene
    ]

    c = chain[
        (
            chain[
                "scene"
            ] == scene
        )
        &
        (
            chain[
                "time_index"
            ] > 0
        )
    ]


    worst_idx = (
        a[
            "aligned_prev_frac_iou_ge_020"
        ].idxmin()
    )

    worst = alignment.loc[
        worst_idx
    ]


    summary_rows.append({
        "scene":
            scene,

        "max_pair_shift_px":
            float(
                a[
                    "pair_shift_magnitude_px"
                ].max()
            ),

        "median_aligned_best_iou":
            float(
                a[
                    "aligned_prev_median_best_iou"
                ].median()
            ),

        "worst_transition_frac_iou_ge_020":
            float(
                worst[
                    "aligned_prev_frac_iou_ge_020"
                ]
            ),

        "worst_transition_time_index":
            int(
                worst[
                    "time_index"
                ]
            ),

        "worst_transition_median_iou":
            float(
                worst[
                    "aligned_prev_median_best_iou"
                ]
            ),

        "max_pair_vs_direct_delta_residual_px":
            float(
                c[
                    "pair_vs_direct_delta_residual_px"
                ].max()
            ),

        "max_cumulative_vs_direct_residual_px":
            float(
                c[
                    "cumulative_vs_direct_residual_px"
                ].max()
            ),
    })


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUT
    / "PRETRACK_REGISTRATION_SUMMARY.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 120)
print("PRE-TRACK REGISTRATION SUMMARY")
print("=" * 120)

print(
    summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("SCENE 0 TRANSITIONS")
print("=" * 120)

cols = [
    "time_index",
    "pair_shift_dy_px",
    "pair_shift_dx_px",
    "pair_shift_magnitude_px",
    "raw_prev_median_best_iou",
    "aligned_prev_median_best_iou",
    "aligned_prev_frac_iou_ge_020",
    "aligned_prev_frac_iou_ge_050",
]

print(
    alignment[
        alignment[
            "scene"
        ] == 0
    ][
        cols
    ].to_string(
        index=False
    )
)


print()
print("=" * 120)
print("LARGEST CHAIN/DIRECT DISAGREEMENTS")
print("=" * 120)

print(
    chain[
        chain[
            "time_index"
        ] > 0
    ]
    .sort_values(
        "pair_vs_direct_delta_residual_px",
        ascending=False,
    )[
        [
            "scene",
            "time_index",
            "pair_dy",
            "pair_dx",
            "direct_delta_dy",
            "direct_delta_dx",
            "pair_vs_direct_delta_residual_px",
            "cumulative_vs_direct_residual_px",
        ]
    ]
    .head(
        20
    )
    .to_string(
        index=False
    )
)


print()
print(
    "Output:",
    OUT
)
