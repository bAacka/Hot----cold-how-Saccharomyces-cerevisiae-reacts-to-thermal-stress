#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import shift as ndi_shift
from scipy.optimize import linear_sum_assignment
from skimage.measure import regionprops_table


ROOT = Path(
    "/"
)

MASK_DIR = (
    ROOT
    / "44_phenotype_bridge"
    / "segmentation"
    / "190321_WT"
    / "all_scenes_halfscale"
    / "scene2"
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
    / "tracking"
    / "190321_WT"
    / "scene2"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


SCENE = 2

MIN_IOU = 0.20
MAX_CENTROID_DISTANCE = 20.0

TIMES = list(range(16))


                                                              
              
 
                                            
                                                                 
                                                              

reg = pd.read_csv(
    REG_FILE,
    sep="\t",
)

reg = reg[
    reg["scene"] == SCENE
].copy()

reg = reg.set_index(
    "time_index"
)


                                                              
            
                                                              

masks = {}

for t in TIMES:

    p = (
        MASK_DIR
        / f"T{t:02d}_masks.npy"
    )

    if not p.exists():
        raise FileNotFoundError(p)

    masks[t] = np.load(p)


                                                              
                
                                                              

def props(mask):

    x = pd.DataFrame(
        regionprops_table(
            mask,
            properties=(
                "label",
                "area",
                "centroid",
            ),
        )
    )

    return x.set_index(
        "label"
    )


props_by_t = {
    t: props(masks[t])
    for t in TIMES
}


                                                              
                       
                                                              

next_track_id = 1

                                    
label_to_track = {}

tracks = []

events = []


                                                              
               
                                                              

for label, row in props_by_t[0].iterrows():

    track_id = next_track_id
    next_track_id += 1

    label_to_track[
        (0, int(label))
    ] = track_id

    tracks.append({
        "track_id":
            track_id,

        "time_index":
            0,

        "label":
            int(label),

        "centroid_y":
            float(row["centroid-0"]),

        "centroid_x":
            float(row["centroid-1"]),

        "area_px":
            float(row["area"]),

        "origin":
            "T0",
    })


print(
    "T0 tracks:",
    len(props_by_t[0])
)


                                                              
                         
                                                              

for t in TIMES[1:]:

    prev_t = t - 1

    prev_mask = masks[
        prev_t
    ]

    curr_mask = masks[
        t
    ]

    prev_props = props_by_t[
        prev_t
    ]

    curr_props = props_by_t[
        t
    ]


                                                              
                                                          
                                   
                                                              

    dy = float(
        reg.loc[
            t,
            "pair_shift_dy_px"
        ]
    )

    dx = float(
        reg.loc[
            t,
            "pair_shift_dx_px"
        ]
    )

    curr_aligned = ndi_shift(
        curr_mask,
        shift=(dy, dx),
        order=0,
        mode="constant",
        cval=0,
        prefilter=False,
    ).astype(
        curr_mask.dtype
    )


    prev_ids = np.array(
        sorted(
            int(x)
            for x in np.unique(prev_mask)
            if x != 0
        ),
        dtype=int,
    )

    curr_ids = np.array(
        sorted(
            int(x)
            for x in np.unique(curr_mask)
            if x != 0
        ),
        dtype=int,
    )


    prev_index = {
        label: i
        for i, label in enumerate(prev_ids)
    }

    curr_index = {
        label: j
        for j, label in enumerate(curr_ids)
    }


                                                              
           
                                                              

    prev_area = {
        int(label):
            int(
                np.sum(
                    prev_mask == label
                )
            )
        for label in prev_ids
    }

    curr_aligned_area = {
        int(label):
            int(
                np.sum(
                    curr_aligned == label
                )
            )
        for label in curr_ids
    }


                                                              
                
                                                              

    iou = np.zeros(
        (
            len(prev_ids),
            len(curr_ids),
        ),
        dtype=np.float32,
    )


    for pid in prev_ids:

        p_region = (
            prev_mask == pid
        )

        candidate_ids, intersections = np.unique(
            curr_aligned[
                p_region
            ],
            return_counts=True,
        )

        for cid, intersection in zip(
            candidate_ids,
            intersections,
        ):

            if cid == 0:
                continue

            cid = int(cid)

            if cid not in curr_index:
                continue

            union = (
                prev_area[int(pid)]
                + curr_aligned_area[cid]
                - int(intersection)
            )

            if union <= 0:
                continue

            iou[
                prev_index[int(pid)],
                curr_index[cid],
            ] = (
                float(intersection)
                / float(union)
            )


                                                              
                                                     
                                                              

    distance = np.full(
        iou.shape,
        np.inf,
        dtype=np.float32,
    )


    for pid in prev_ids:

        py = float(
            prev_props.loc[
                pid,
                "centroid-0"
            ]
        )

        px = float(
            prev_props.loc[
                pid,
                "centroid-1"
            ]
        )

        for cid in curr_ids:

            cy = float(
                curr_props.loc[
                    cid,
                    "centroid-0"
                ]
                + dy
            )

            cx = float(
                curr_props.loc[
                    cid,
                    "centroid-1"
                ]
                + dx
            )

            distance[
                prev_index[int(pid)],
                curr_index[int(cid)],
            ] = np.hypot(
                py - cy,
                px - cx,
            )


                                                              
                     
     
                    
                                     
                                                              

    cost = (
        1.0 - iou
    )

    cost += (
        np.minimum(
            distance,
            MAX_CENTROID_DISTANCE,
        )
        / MAX_CENTROID_DISTANCE
        * 0.10
    )


    rows, cols = linear_sum_assignment(
        cost
    )


    matched_prev = set()
    matched_curr = set()


    for r, c in zip(
        rows,
        cols,
    ):

        pid = int(
            prev_ids[r]
        )

        cid = int(
            curr_ids[c]
        )

        this_iou = float(
            iou[r, c]
        )

        this_distance = float(
            distance[r, c]
        )


                                              
                                               
        accept = (
            this_iou >= MIN_IOU
            or
            (
                this_iou > 0
                and
                this_distance
                <= MAX_CENTROID_DISTANCE
            )
        )

        if not accept:
            continue


        track_id = label_to_track[
            (
                prev_t,
                pid,
            )
        ]

        label_to_track[
            (
                t,
                cid,
            )
        ] = track_id


        row = curr_props.loc[
            cid
        ]

        tracks.append({
            "track_id":
                track_id,

            "time_index":
                t,

            "label":
                cid,

            "centroid_y":
                float(
                    row[
                        "centroid-0"
                    ]
                ),

            "centroid_x":
                float(
                    row[
                        "centroid-1"
                    ]
                ),

            "area_px":
                float(
                    row[
                        "area"
                    ]
                ),

            "origin":
                "continued",
        })


        events.append({
            "time_index":
                t,

            "event":
                "match",

            "track_id":
                track_id,

            "previous_label":
                pid,

            "current_label":
                cid,

            "iou":
                this_iou,

            "centroid_distance_px":
                this_distance,
        })


        matched_prev.add(
            pid
        )

        matched_curr.add(
            cid
        )


                                                              
                                            
                                                              

    births = 0

    for cid in curr_ids:

        cid = int(cid)

        if cid in matched_curr:
            continue


        track_id = next_track_id
        next_track_id += 1

        label_to_track[
            (
                t,
                cid,
            )
        ] = track_id


        row = curr_props.loc[
            cid
        ]

        tracks.append({
            "track_id":
                track_id,

            "time_index":
                t,

            "label":
                cid,

            "centroid_y":
                float(
                    row[
                        "centroid-0"
                    ]
                ),

            "centroid_x":
                float(
                    row[
                        "centroid-1"
                    ]
                ),

            "area_px":
                float(
                    row[
                        "area"
                    ]
                ),

            "origin":
                "birth",
        })


        events.append({
            "time_index":
                t,

            "event":
                "birth",

            "track_id":
                track_id,

            "previous_label":
                np.nan,

            "current_label":
                cid,

            "iou":
                np.nan,

            "centroid_distance_px":
                np.nan,
        })


        births += 1


    deaths = (
        len(prev_ids)
        - len(matched_prev)
    )


    print(
        f"T={prev_t:02d}->{t:02d}",
        f"prev={len(prev_ids)}",
        f"curr={len(curr_ids)}",
        f"matched={len(matched_curr)}",
        f"births={births}",
        f"unmatched_prev={deaths}",
        f"shift=({dy:+.1f},{dx:+.1f})",
    )


                                                              
        
                                                              

tracks = pd.DataFrame(
    tracks
)

events = pd.DataFrame(
    events
)


tracks.to_csv(
    OUT
    / "CELL_TRACKS_LONG.tsv",
    sep="\t",
    index=False,
)

events.to_csv(
    OUT
    / "TRACKING_EVENTS.tsv",
    sep="\t",
    index=False,
)


                                                              
               
                                                              

summary = (
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

        n_frames=(
            "time_index",
            "size",
        ),

        median_area_px=(
            "area_px",
            "median",
        ),
    )
    .reset_index()
)


summary[
    "is_T0_cohort"
] = (
    summary[
        "first_time"
    ] == 0
)


summary[
    "survives_to_T15"
] = (
    summary[
        "last_time"
    ] == 15
)


summary.to_csv(
    OUT
    / "TRACK_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
                     
                                                              

t0 = summary[
    summary[
        "is_T0_cohort"
    ]
].copy()


retention_rows = []

for t in TIMES:

    present_ids = set(
        tracks.loc[
            tracks[
                "time_index"
            ] == t,
            "track_id"
        ]
    )

    n_present = sum(
        int(track_id in present_ids)
        for track_id in t0[
            "track_id"
        ]
    )

    retention_rows.append({
        "time_index":
            t,

        "T0_tracks_present":
            n_present,

        "T0_tracks_total":
            len(t0),

        "fraction_retained":
            n_present
            / len(t0),
    })


retention = pd.DataFrame(
    retention_rows
)


retention.to_csv(
    OUT
    / "T0_COHORT_RETENTION.tsv",
    sep="\t",
    index=False,
)


                                                              
                             
                                                              

match_quality = (
    events[
        events[
            "event"
        ] == "match"
    ]
    .groupby(
        "time_index"
    )
    .agg(
        n_matches=(
            "track_id",
            "size",
        ),

        median_IoU=(
            "iou",
            "median",
        ),

        q05_IoU=(
            "iou",
            lambda x:
                x.quantile(
                    0.05
                ),
        ),

        median_distance_px=(
            "centroid_distance_px",
            "median",
        ),

        q95_distance_px=(
            "centroid_distance_px",
            lambda x:
                x.quantile(
                    0.95
                ),
        ),
    )
    .reset_index()
)


match_quality.to_csv(
    OUT
    / "MATCH_QUALITY_BY_TRANSITION.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.plot(
    retention[
        "time_index"
    ],
    retention[
        "fraction_retained"
    ],
    marker="o",
)

ax.set_ylim(
    0,
    1.05,
)

ax.set_xlabel(
    "CZI time index"
)

ax.set_ylabel(
    "Fraction of T0 cohort retained"
)

fig.tight_layout()

fig.savefig(
    OUT
    / "T0_COHORT_RETENTION.png",
    dpi=160,
)

plt.close(
    fig
)


fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.plot(
    match_quality[
        "time_index"
    ],
    match_quality[
        "median_IoU"
    ],
    marker="o",
)

ax.set_xlabel(
    "Current frame"
)

ax.set_ylabel(
    "Median adjacent-frame IoU"
)

fig.tight_layout()

fig.savefig(
    OUT
    / "MATCH_IOU_BY_TRANSITION.png",
    dpi=160,
)

plt.close(
    fig
)


                                                              
                 
                                                              

print()
print("=" * 120)
print("T0 COHORT RETENTION")
print("=" * 120)

print(
    retention.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("MATCH QUALITY")
print("=" * 120)

print(
    match_quality.to_string(
        index=False
    )
)


print()
print("=" * 120)

print(
    "Total persistent tracks:",
    len(summary)
)

print(
    "T0 tracks:",
    len(t0)
)

print(
    "T0 surviving through T15:",
    int(
        t0[
            "survives_to_T15"
        ].sum()
    )
)

print(
    "Fraction T0 -> T15:",
    float(
        t0[
            "survives_to_T15"
        ].mean()
    )
)

print()
print(
    "Output:",
    OUT
)
