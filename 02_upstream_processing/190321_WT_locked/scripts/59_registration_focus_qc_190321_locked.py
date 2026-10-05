#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import laplace
from skimage.filters import sobel
from skimage.registration import phase_cross_correlation

from aicspylibczi import CziFile


ROOT = Path(
    "/"
)

DATA = (
    ROOT
    / "44_phenotype_bridge"
    / "downloads"
    / "190321_WT"
)

OUT = (
    ROOT
    / "44_phenotype_bridge"
    / "qc"
    / "190321_WT"
    / "registration_focus"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
            
                                                              

files = sorted(
    DATA.rglob("*.czi")
)

if len(files) != 1:
    raise RuntimeError(
        f"Expected exactly one CZI; found {len(files)}"
    )

path = files[0]

czi = CziFile(
    path
)

shape = czi.get_dims_shape()[0]


def indices(dim):
    start, n = shape[dim]

    return list(
        range(
            int(start),
            int(start + n),
        )
    )


SCENES = indices("S")
TIMES = indices("T")
ZS = indices("Z")

FIXED_DIC_Z = ZS[
    len(ZS) // 2
]


print("CZI:", path)
print("Scenes:", SCENES)
print("Times:", TIMES)
print("Z:", ZS)
print("Fixed registration DIC Z:", FIXED_DIC_Z)


                                                              
                 
                                                              

def read_plane(
    scene,
    time,
    channel,
    z,
):

    img, dims = czi.read_image(
        S=scene,
        T=time,
        C=channel,
        Z=z,
    )

    x = np.squeeze(
        np.asarray(img)
    )

    if x.ndim != 2:
        raise RuntimeError(
            f"Expected YX; got {x.shape}"
        )

    return x.astype(
        np.float32,
        copy=False,
    )


                                                              
                                                            
 
                                                             
                                                 
                                                              

def robust_normalize(x):

    lo, hi = np.percentile(
        x,
        [2, 98],
    )

    if hi <= lo:
        hi = lo + 1

    y = (
        x - lo
    ) / (
        hi - lo
    )

    return np.clip(
        y,
        0,
        1,
    )


def registration_image(x):

    y = robust_normalize(
        x
    )

                                                  
                                  
    y = sobel(
        y
    )

                                                         
    margin = 50

    if (
        y.shape[0] > 2 * margin
        and
        y.shape[1] > 2 * margin
    ):
        y = y[
            margin:-margin,
            margin:-margin,
        ]

    return y


def focus_score(x):

    y = robust_normalize(
        x
    )

                            
                                              
    l = laplace(
        y
    )

    return float(
        np.var(
            l
        )
    )


                                                              
                  
                                                              

rows = []


for scene in SCENES:

    print()
    print("=" * 100)
    print("SCENE", scene)
    print("=" * 100)


                                                              
                                                    
                                                              

    reg_images = {}

    for t in TIMES:

        dic = read_plane(
            scene=scene,
            time=t,
            channel=0,
            z=FIXED_DIC_Z,
        )

        reg_images[t] = (
            registration_image(
                dic
            )
        )


                                                              
                                                
                                                              

    cumulative_y = 0.0
    cumulative_x = 0.0


    for ti, t in enumerate(
        TIMES
    ):

                                                  
                                                  
                                                  

        focus_by_z = []

        for z in ZS:

            dic_z = read_plane(
                scene=scene,
                time=t,
                channel=0,
                z=z,
            )

            fs = focus_score(
                dic_z
            )

            focus_by_z.append(
                (
                    z,
                    fs,
                )
            )


        best_z, best_score = max(
            focus_by_z,
            key=lambda q:
                q[1],
        )


        fixed_score = dict(
            focus_by_z
        )[
            FIXED_DIC_Z
        ]


                                                  
                      
                                                  

        if ti == 0:

            pair_dy = 0.0
            pair_dx = 0.0
            pair_error = 0.0

            direct_dy = 0.0
            direct_dx = 0.0
            direct_error = 0.0

        else:

            prev_t = TIMES[
                ti - 1
            ]


            pair_shift, pair_error, _ = (
                phase_cross_correlation(
                    reg_images[
                        prev_t
                    ],
                    reg_images[
                        t
                    ],
                    upsample_factor=10,
                    normalization="phase",
                )
            )


            pair_dy = float(
                pair_shift[0]
            )

            pair_dx = float(
                pair_shift[1]
            )


            cumulative_y += pair_dy
            cumulative_x += pair_dx


            direct_shift, direct_error, _ = (
                phase_cross_correlation(
                    reg_images[
                        TIMES[0]
                    ],
                    reg_images[
                        t
                    ],
                    upsample_factor=10,
                    normalization="phase",
                )
            )


            direct_dy = float(
                direct_shift[0]
            )

            direct_dx = float(
                direct_shift[1]
            )


        rec = {
            "scene":
                scene,

            "time_index":
                t,

            "best_focus_z":
                int(
                    best_z
                ),

            "best_focus_score":
                best_score,

            "fixed_z":
                FIXED_DIC_Z,

            "fixed_z_focus_score":
                fixed_score,

            "pair_shift_dy_px":
                pair_dy,

            "pair_shift_dx_px":
                pair_dx,

            "pair_registration_error":
                float(
                    pair_error
                ),

            "cumulative_shift_dy_px":
                cumulative_y,

            "cumulative_shift_dx_px":
                cumulative_x,

            "cumulative_shift_magnitude_px":
                float(
                    np.hypot(
                        cumulative_y,
                        cumulative_x,
                    )
                ),

            "direct_T0_shift_dy_px":
                direct_dy,

            "direct_T0_shift_dx_px":
                direct_dx,

            "direct_T0_shift_magnitude_px":
                float(
                    np.hypot(
                        direct_dy,
                        direct_dx,
                    )
                ),

            "direct_T0_registration_error":
                float(
                    direct_error
                ),
        }


        rows.append(
            rec
        )


        print(
            f"T={t:02d}",
            f"bestZ={best_z}",
            f"focus={best_score:.6g}",
            f"pair=({pair_dy:+.2f},{pair_dx:+.2f})",
            f"cum=({cumulative_y:+.2f},{cumulative_x:+.2f})",
            f"|cum|={rec['cumulative_shift_magnitude_px']:.2f}",
        )


qc = pd.DataFrame(
    rows
)


qc.to_csv(
    OUT
    / "REGISTRATION_FOCUS_QC.tsv",
    sep="\t",
    index=False,
)


                                                              
                   
                                                              

summary = (
    qc.groupby(
        "scene",
        as_index=False,
    )
    .agg(
        n_timepoints=(
            "time_index",
            "size",
        ),

        best_z_min=(
            "best_focus_z",
            "min",
        ),

        best_z_max=(
            "best_focus_z",
            "max",
        ),

        median_best_z=(
            "best_focus_z",
            "median",
        ),

        max_pair_shift_px=(
            "pair_shift_dy_px",
            lambda x:
                np.nan,
        ),

        max_cumulative_shift_px=(
            "cumulative_shift_magnitude_px",
            "max",
        ),

        max_direct_T0_shift_px=(
            "direct_T0_shift_magnitude_px",
            "max",
        ),
    )
)


                                                  
pair_mag = (
    np.hypot(
        qc[
            "pair_shift_dy_px"
        ],
        qc[
            "pair_shift_dx_px"
        ],
    )
)

qc[
    "pair_shift_magnitude_px"
] = pair_mag


pair_summary = (
    qc.groupby(
        "scene"
    )[
        "pair_shift_magnitude_px"
    ]
    .max()
)


summary[
    "max_pair_shift_px"
] = (
    summary[
        "scene"
    ].map(
        pair_summary
    )
)


summary.to_csv(
    OUT
    / "REGISTRATION_FOCUS_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                              
qc.to_csv(
    OUT
    / "REGISTRATION_FOCUS_QC.tsv",
    sep="\t",
    index=False,
)


                                                              
                  
                                                              

for metric, ylabel, fn in [
    (
        "cumulative_shift_magnitude_px",
        "Cumulative XY shift (pixels)",
        "cumulative_drift.png",
    ),
    (
        "best_focus_z",
        "Best-focus DIC Z",
        "best_focus_z.png",
    ),
    (
        "best_focus_score",
        "DIC focus score",
        "focus_score.png",
    ),
]:

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    for scene in SCENES:

        d = qc[
            qc["scene"]
            == scene
        ]

        ax.plot(
            d[
                "time_index"
            ],
            d[
                metric
            ],
            marker="o",
            label=f"Scene {scene}",
        )


    ax.set_xlabel(
        "CZI time index"
    )

    ax.set_ylabel(
        ylabel
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        OUT
        / fn,
        dpi=160,
    )

    plt.close(
        fig
    )


print()
print("=" * 120)
print("SCENE SUMMARY")
print("=" * 120)

print(
    summary.to_string(
        index=False
    )
)

print()
print(
    "Output:",
    OUT,
)
