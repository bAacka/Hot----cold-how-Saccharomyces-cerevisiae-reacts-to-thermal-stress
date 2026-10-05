#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(
    "/"
)

INFILE = (
    ROOT
    / "44_phenotype_bridge"
    / "dataset_audit"
    / "190321_WT"
    / "timestamps"
    / "FRAME_TIMESTAMPS.tsv"
)

OUT = INFILE.parent


d = pd.read_csv(
    INFILE,
    sep="\t",
)


                                                              
                                                       
                                                              

summary = (
    d.groupby(
        "time_index"
    )[
        "elapsed_min_from_scene_first"
    ]
    .agg(
        actual_min_mean="mean",
        actual_min_median="median",
        actual_min_min="min",
        actual_min_max="max",
    )
    .reset_index()
)


summary[
    "biological_min"
] = (
    summary[
        "time_index"
    ]
    * 10
)


summary[
    "manual_S10_timepoint"
] = (
    summary[
        "biological_min"
    ].between(
        10,
        120,
    )
)


summary[
    "manual_S10_min"
] = np.where(
    summary[
        "manual_S10_timepoint"
    ],
    summary[
        "biological_min"
    ],
    np.nan,
)


summary[
    "abs_error_from_nominal_min"
] = np.abs(
    summary[
        "actual_min_median"
    ]
    - summary[
        "biological_min"
    ]
)


summary.to_csv(
    OUT
    / "TIME_MAPPING_LOCKED.tsv",
    sep="\t",
    index=False,
)


manual = summary[
    summary[
        "manual_S10_timepoint"
    ]
].copy()


if manual[
    "time_index"
].tolist() != list(
    range(
        1,
        13,
    )
):
    raise RuntimeError(
        "Unexpected manual validation frame mapping"
    )


if manual[
    "biological_min"
].tolist() != list(
    range(
        10,
        121,
        10,
    )
):
    raise RuntimeError(
        "Unexpected biological minute mapping"
    )


max_error = float(
    summary[
        "abs_error_from_nominal_min"
    ].max()
)


freeze = f"""190321 WT CZI time-axis freeze

CZI:
16 timepoints, T0 through T15.

Per-frame acquisition timestamps were recovered independently
from CZI subblock metadata for all four scenes.

Observed frame cadence:
approximately 9.99 minutes.

Locked biological mapping:
T0  =   0 min
T1  =  10 min
T2  =  20 min
...
T12 = 120 min
T13 = 130 min
T14 = 140 min
T15 = 150 min

Published Andersson S10 WT manual measurements:
10 through 120 min.

Therefore held-out manual-validation comparison uses:
CZI T1 through T12
versus
S10 10 through 120 min.

CZI T0 is an additional 0-min recovery frame and is NOT part
of the published S10 manual reference.

CZI T13-T15 are additional 130-150 min frames and are NOT part
of the published S10 manual reference.

Maximum absolute deviation of median actual scene elapsed time
from nominal 10-minute biological time:
{max_error:.6f} min.

Do not shift the S10 trajectory by one frame during validation.
"""


(
    OUT
    / "TIME_MAPPING_FREEZE.txt"
).write_text(
    freeze
)


print("=" * 110)
print("LOCKED TIME MAPPING")
print("=" * 110)

print(
    summary.to_string(
        index=False
    )
)

print()
print(freeze)
