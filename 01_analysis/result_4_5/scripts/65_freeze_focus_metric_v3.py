#!/usr/bin/env python3

from pathlib import Path
import hashlib

import numpy as np
import pandas as pd


ROOT = Path(
    "/"
)

DEV = (
    ROOT
    / "44_phenotype_bridge"
    / "focus_metric_development"
    / "190320_WT"
    / "scene1_v3_transport"
)

INIT_FILE = (
    ROOT
    / "44_phenotype_bridge"
    / "focus_metric_development"
    / "190320_WT"
    / "scene1_v2_tracked"
    / "T1_INITIALIZATION.tsv"
)

LONG_FILE = (
    DEV
    / "TRANSPORTED_FOCUS_LONG_WITH_T1_NORMALIZATION.tsv"
)

OUT = (
    ROOT
    / "44_phenotype_bridge"
    / "focus_metric_freeze"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                                    
 
                                                  
                                                     
                                                              

PRIMARY_Z_THRESHOLD = 5.0

SENSITIVITY_THRESHOLDS = [
    3.0,
    4.0,
    5.0,
    6.0,
]


init = pd.read_csv(
    INIT_FILE,
    sep="\t",
)

long = pd.read_csv(
    LONG_FILE,
    sep="\t",
)


                                                              
                                                
                                                              

rows = []


for threshold in SENSITIVITY_THRESHOLDS:

    ids = set(
        init.loc[
            init[
                "initial_focus_prominence_z"
            ] >= threshold,
            "track_id",
        ].astype(int)
    )

    x = long[
        long[
            "track_id"
        ].isin(ids)
    ].copy()


    def vals(t, col):

        return x.loc[
            x[
                "time_index"
            ] == t,
            col,
        ].dropna()


    t1 = vals(
        1,
        "focus_to_cytosol_ratio",
    )

    t12 = vals(
        12,
        "focus_to_cytosol_ratio",
    )

    t14 = vals(
        14,
        "focus_to_cytosol_ratio",
    )

    rel12 = vals(
        12,
        "ratio_relative_to_T1",
    )

    rel14 = vals(
        14,
        "ratio_relative_to_T1",
    )


    wide = (
        x[
            x[
                "time_index"
            ].isin(
                [
                    1,
                    12,
                    14,
                ]
            )
        ][
            [
                "track_id",
                "time_index",
                "focus_to_cytosol_ratio",
            ]
        ]
        .pivot(
            index="track_id",
            columns="time_index",
            values="focus_to_cytosol_ratio",
        )
    )


    frac_down_12 = float(
        (
            wide[
                12
            ]
            < wide[
                1
            ]
        ).mean()
    )

    frac_down_14 = float(
        (
            wide[
                14
            ]
            < wide[
                1
            ]
        ).mean()
    )


    rows.append({
        "z_threshold":
            threshold,

        "n_focus_positive":
            len(ids),

        "fraction_of_116":
            len(ids)
            / 116.0,

        "median_ratio_T1":
            float(
                t1.median()
            ),

        "median_ratio_T12":
            float(
                t12.median()
            ),

        "median_ratio_T14":
            float(
                t14.median()
            ),

        "median_T12_over_T1":
            float(
                rel12.median()
            ),

        "median_T14_over_T1":
            float(
                rel14.median()
            ),

        "fraction_tracks_lower_T12_than_T1":
            frac_down_12,

        "fraction_tracks_lower_T14_than_T1":
            frac_down_14,
    })


audit = pd.DataFrame(
    rows
)

audit.to_csv(
    OUT
    / "DEVELOPMENT_THRESHOLD_SENSITIVITY.tsv",
    sep="\t",
    index=False,
)


                                                              
                       
                                                              

primary = init[
    init[
        "initial_focus_prominence_z"
    ]
    >= PRIMARY_Z_THRESHOLD
].copy()

primary = primary.sort_values(
    [
        "initial_focus_prominence_z",
        "track_id",
    ],
    ascending=[
        False,
        True,
    ],
)

primary.to_csv(
    OUT
    / "190320_DEVELOPMENT_FOCUS_POSITIVE_T1.tsv",
    sep="\t",
    index=False,
)


                                                              
                      
                                                              

freeze_text = f"""AUTOMATED GUK1-7-GFP RECOVERY PHENOTYPE -- PRE-VALIDATION FREEZE

DEVELOPMENT DATA
190320 WT
Scene 1 only.

HELD-OUT VALIDATION DATA
190321 WT.
The automated V3 phenotype has NOT yet been applied to 190321
at the time of this freeze.

BIOLOGICAL QUANTITY
Loss of GFP enrichment at the original T1 focus location during
post-heat-shock recovery.

CELL SEGMENTATION
Cellpose cpsam_v2.
DIC Z=5.
0.5x inference with categorical nearest-neighbour restoration.

CELL TRACKING
Registration-aware adjacent-frame Hungarian assignment.
MIN_IOU = 0.20.
MAX_CENTROID_DISTANCE = 20 px.

FLUORESCENCE IMAGE
Raw EGFP 10-plane Z maximum projection.
No spatial registration/resampling of GFP.

BACKGROUND CORRECTION
One per-frame extracellular median background subtraction from
the raw GFP maximum projection.
Background pixels are outside a 5-px dilation of all cell masks.
No clipping of negative background-corrected pixel deviations.

T1 FOCUS INITIALIZATION
Localized signal:
Gaussian sigma 1 px minus Gaussian sigma 6 px.

Focus ROI:
radius 4 px.
area 49 px.
full disk required inside cell.

Cytosol ROI:
same area as focus ROI.
At T1 selected from eligible intracellular locations with centre
at least 10 px from the focus, using the candidate whose raw
background-corrected local mean is nearest the median eligible
cytosolic local mean.

LONGITUDINAL RULE
After T1, fluorescence intensity has ZERO influence on ROI
position.

Focus and cytosol ROIs are transported using tracked cell-centroid
displacement and snapped to the geometrically nearest valid
full-disk intracellular position.
The cytosol ROI remains constrained not to overlap the focus ROI.

PRIMARY PER-CELL PHENOTYPE
focus_mean_bg_corrected / cytosol_mean_bg_corrected

SECONDARY PHENOTYPES
focus_minus_cytosol
bounded_focus_contrast
current_focus_prominence_z
within-cell focus/cytosol ratio relative to T1

FOCUS-POSITIVE INCLUSION
initial T1 robust localized-focus prominence z >= {PRIMARY_Z_THRESHOLD:.1f}

RATIONALE FOR THRESHOLD
Selected conservatively from visual QC of the 190320 development
movie before applying the automated method to held-out 190321.
The visible transition from convincing compact puncta to ambiguous
or diffuse signal occurred around z=4-5; z>=5 was chosen as the
high-confidence boundary.

VALIDATION WINDOW
For 190321:
T1 through T12 correspond to 10 through 120 min.
T0 is not part of the published S10 comparison.

VALIDATION LEVEL
Population/distribution and trajectory validation.
Do NOT claim one-to-one correspondence between automated cells
and the 20 manually scored S10 cells because manual ROI coordinates
are unavailable.

STATUS
FROZEN BEFORE AUTOMATED HELD-OUT VALIDATION.
"""


(
    OUT
    / "FOCUS_METRIC_V3_PREVALIDATION_FREEZE.txt"
).write_text(
    freeze_text
)


print()
print("=" * 120)
print("DEVELOPMENT THRESHOLD SENSITIVITY")
print("=" * 120)

print(
    audit.to_string(
        index=False
    )
)

print()
print(
    "PRIMARY z threshold:",
    PRIMARY_Z_THRESHOLD,
)

print(
    "Primary development focus-positive n:",
    len(
        primary
    ),
)

print()
print(freeze_text)


                                                              
                                  
                                                              

files_to_hash = [
    OUT
    / "DEVELOPMENT_THRESHOLD_SENSITIVITY.tsv",

    OUT
    / "190320_DEVELOPMENT_FOCUS_POSITIVE_T1.tsv",

    OUT
    / "FOCUS_METRIC_V3_PREVALIDATION_FREEZE.txt",

    ROOT
    / "scripts"
    / "64_transport_focus_metric_190320_v3.py",
]


with (
    OUT
    / "FOCUS_METRIC_V3_PREVALIDATION_SHA256SUMS.txt"
).open(
    "w"
) as fh:

    for p in files_to_hash:

        h = hashlib.sha256(
            p.read_bytes()
        ).hexdigest()

        fh.write(
            f"{h}  {p}\n"
        )


print(
    "Freeze:",
    OUT
)
