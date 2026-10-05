#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import math
import sys

import av
import cv2
import numpy as np
import pandas as pd
from aicspylibczi import CziFile


                                                                       
       
                                                                       

ROOT = Path("/")
OUT = ROOT / "44_phenotype_bridge"
RES = OUT / "46_phenotype_state/46r_hsp104_genotype_resolution"
MOV = RES / "published_movies"

CACHE = RES / "46r2_raw_zmax_cache"
CACHE.mkdir(parents=True, exist_ok=True)

            
                                                                 
 
                             
                                                   
                                                                       

CANDIDATES = {
    "A": (
        OUT
        / "downloads/190417_WT_and_hsp104D/extracted/"
          "190417 wt hsp104D guk1-7-GFP timelaps/"
          "190417 WT hsp104D pguk1-7-GFP recovery timelapse -05.czi"
    ),
    "B": (
        OUT
        / "downloads/190418_WT_and_hsp104D/extracted/"
          "190418 WT hsp104D guk1-7-GFP timelaps/"
          "190418 WT hsp104D pguk1-7-GFP recovery timelapse -08.czi"
    ),
}

MOVIES = {
    "M1": MOV / "S1_WT_guk1-7-GFP.avi",
    "M2": MOV / "S4_hsp104D_guk1-7-GFP.avi",
}


                                                                       
                                 
 
                                                           
                                                                       

IMPL_FREEZE = RES / "46r2_IMPLEMENTATION_FREEZE_PRE_MATCH.txt"

IMPL_FREEZE.write_text(
"""HSP104 RAW-CZI ↔ MOVIE MATCH — IMPLEMENTATION FREEZE
=======================================================

This implementation is frozen before calculating any CZI-to-movie
similarity score.

BLINDING
--------
Raw acquisitions:
    A = 190417 long acquisition
    B = 190418 long acquisition

Published movies:
    M1 = published movie 1
    M2 = published movie 2

Genotype labels are not used during score computation.

RAW REPRESENTATION
------------------
Exactly as previously validated:

    C = 1
    all Z=0..9
    pixelwise Z maximum
    S=0..3

MOVIE REPRESENTATION
--------------------
12 grayscale AVI frames.

BIOLOGICAL QUANTITIES NOT USED
------------------------------
No segmentation.
No focus count.
No focus persistence.
No fluorescence-summary comparison.
No recovery-direction comparison.

SPATIAL SEARCH
--------------
For every movie:
    every raw candidate
    every scene
    every contiguous 12-frame raw time window

Publication scale is searched by an isotropic downsampling factor.

Coarse scale grid:
    1.0 through 5.7 inclusive
    step 0.2

The first movie frame determines the global crop location at a
particular scale.

Frames 6 and 12 (zero-based indices 5 and 11) must match near the
same crop location.

Allowed local movement:
    +/- 12 pixels in published-image coordinates.

Coarse hypothesis score:
    median normalized spatial correlation across movie frames
    0, 5, 11.

REFINEMENT
----------
Top coarse hypotheses for each movie are refined over:

    coarse_scale +/- 0.20
    step 0.02

The complete 12-frame series is then scored using:

    same raw acquisition
    same raw scene
    same contiguous time window
    same isotropic scale
    crop constrained to +/-12 pixels around the anchor location.

PRIMARY FINAL SCORE
-------------------
Median normalized spatial correlation across all 12 frames.

Supporting diagnostics:
    minimum correlation
    25th percentile correlation
    maximum correlation
    crop-coordinate variability

GENOTYPE
--------
No genotype assignment is performed by this script.
Only blinded M1/M2 versus A/B identity scores are written.

""",
    encoding="utf-8",
)

freeze_sha = hashlib.sha256(
    IMPL_FREEZE.read_bytes()
).hexdigest()

print("implementation freeze sha256:", freeze_sha)


                                                                       
                                                              
                                  
                                                                       

cv2.setNumThreads(1)

                                   
_test_big = np.arange(
    100 * 120,
    dtype=np.float32
).reshape(100, 120)

_test_small = cv2.resize(
    _test_big,
    (60, 50),
    interpolation=cv2.INTER_AREA,
)

_test_tpl = _test_small[10:30, 15:40].copy()

_test_res = cv2.matchTemplate(
    _test_small,
    _test_tpl,
    cv2.TM_CCOEFF_NORMED,
)

if not np.isfinite(_test_res).all():
    raise RuntimeError(
        "OpenCV matchTemplate smoke test failed."
    )

print("OpenCV image-operation smoke test: PASS")


                                                                       
              
                                                                       

def dim_sizes(czi):
    dims = czi.get_dims_shape()

    if len(dims) != 1:
        raise RuntimeError(
            f"Expected one CZI dimension block; got {dims}"
        )

    d = dims[0]

    out = {}

    for name, (a, b) in d.items():
        out[name] = b - a

    return out


def read_selected(czi, **kwargs):
    arr, shape_desc = czi.read_image(
        **kwargs
    )

    names = [
        x[0]
        for x in shape_desc
    ]

    sizes = [
        x[1]
        for x in shape_desc
    ]

    arr = np.asarray(arr)

    if tuple(arr.shape) != tuple(sizes):
        raise RuntimeError(
            f"Shape mismatch: "
            f"ndarray={arr.shape}; "
            f"desc={shape_desc}"
        )

    index = []
    kept = []

    for name, size in zip(
        names,
        sizes,
    ):

        if name in {"Z", "Y", "X"}:
            index.append(
                slice(None)
            )
            kept.append(
                name
            )

        else:
            if size != 1:
                raise RuntimeError(
                    f"Unexpected non-singleton "
                    f"selected axis {name}={size}; "
                    f"kwargs={kwargs}"
                )

            index.append(
                0
            )

    arr = arr[
        tuple(index)
    ]

    return arr, kept


def as_zyx(czi, **kwargs):
    arr, names = read_selected(
        czi,
        **kwargs,
    )

    if set(names) != {
        "Z",
        "Y",
        "X",
    }:
        raise RuntimeError(
            f"Expected Z/Y/X; "
            f"got names={names}, "
            f"shape={arr.shape}"
        )

    order = [
        names.index(k)
        for k in (
            "Z",
            "Y",
            "X",
        )
    ]

    return np.transpose(
        arr,
        order,
    )


                                                                       
                      
 
                                                           
                                                                       

def cache_candidate(alias, path):

    if not path.exists():
        raise FileNotFoundError(
            path
        )

    cache_path = (
        CACHE
        / f"{alias}_EGFP_C1_ZMAX.npy"
    )

    meta_path = (
        CACHE
        / f"{alias}_EGFP_C1_ZMAX.json"
    )

    czi = CziFile(
        path
    )

    dims = dim_sizes(
        czi
    )

    required = {
        "S": 4,
        "C": 2,
        "Z": 10,
        "Y": 1104,
        "X": 1376,
    }

    for k, expected in required.items():
        got = dims.get(k)

        if got != expected:
            raise RuntimeError(
                f"{alias}: unexpected {k}: "
                f"{got} != {expected}"
            )

    nt = dims["T"]

    expected_nt = {
        "A": 19,
        "B": 16,
    }[alias]

    if nt != expected_nt:
        raise RuntimeError(
            f"{alias}: expected T={expected_nt}; "
            f"got T={nt}"
        )

    expected_shape = (
        4,
        nt,
        1104,
        1376,
    )

    if (
        cache_path.exists()
        and meta_path.exists()
    ):

        z = np.load(
            cache_path,
            mmap_mode="r",
        )

        if z.shape != expected_shape:
            raise RuntimeError(
                f"Existing cache has bad shape: "
                f"{z.shape}"
            )

        print(
            f"{alias}: using existing cache "
            f"{cache_path}"
        )

        return cache_path

    print()
    print("=" * 80)
    print(
        f"CACHING {alias}"
    )
    print("=" * 80)

    mm = np.lib.format.open_memmap(
        cache_path,
        mode="w+",
        dtype=np.uint16,
        shape=expected_shape,
    )

    for s in range(4):

        for t in range(nt):

            stack = as_zyx(
                czi,
                S=s,
                T=t,
                C=1,
            )

            if stack.shape != (
                10,
                1104,
                1376,
            ):
                raise RuntimeError(
                    f"{alias} S={s} T={t}: "
                    f"unexpected stack "
                    f"{stack.shape}"
                )

            mx = np.max(
                stack,
                axis=0,
            )

            mm[
                s,
                t,
            ] = mx

            print(
                f"{alias}: "
                f"S={s} "
                f"T={t:02d}/{nt-1:02d}",
                flush=True,
            )

    mm.flush()

    metadata = {
        "alias": alias,
        "source_czi": str(path),
        "shape": list(expected_shape),
        "dtype": "uint16",
        "channel": 1,
        "z_operation": "pixelwise_max",
        "z_planes": list(
            range(10)
        ),
    }

    meta_path.write_text(
        json.dumps(
            metadata,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return cache_path


                                                                       
               
                                                                       

def load_movie(path):

    container = av.open(
        str(path)
    )

    stream = (
        container
        .streams
        .video[0]
    )

    frames = []

    for frame in container.decode(
        stream
    ):

        x = frame.to_ndarray(
            format="gray"
        )

        x = np.asarray(
            x,
            dtype=np.float32,
        )

        frames.append(
            x
        )

    container.close()

    if len(frames) != 12:
        raise RuntimeError(
            f"{path}: expected 12 frames; "
            f"got {len(frames)}"
        )

    if any(
        x.shape != (186, 240)
        for x in frames
    ):
        raise RuntimeError(
            f"{path}: unexpected AVI frame shapes"
        )

    return frames


                                                                       
                
                                                                       

MOV_H = 186
MOV_W = 240

LOCAL_RADIUS = 12


def resize_raw(raw, scale):

    h, w = raw.shape

    new_w = int(
        round(
            w / scale
        )
    )

    new_h = int(
        round(
            h / scale
        )
    )

    if (
        new_w < MOV_W
        or new_h < MOV_H
    ):
        return None

    out = cv2.resize(
        raw.astype(
            np.float32
        ),
        (
            new_w,
            new_h,
        ),
        interpolation=cv2.INTER_AREA,
    )

    return out


def global_match(raw_scaled, movie_frame):

    result = cv2.matchTemplate(
        raw_scaled,
        movie_frame,
        cv2.TM_CCOEFF_NORMED,
    )

    _, maxval, _, maxloc = (
        cv2.minMaxLoc(
            result
        )
    )

    return (
        float(maxval),
        int(maxloc[0]),
        int(maxloc[1]),
    )


def local_match(
    raw_scaled,
    movie_frame,
    anchor_x,
    anchor_y,
    radius=LOCAL_RADIUS,
):

    H, W = raw_scaled.shape

    x0 = max(
        0,
        anchor_x - radius,
    )

    y0 = max(
        0,
        anchor_y - radius,
    )

    x1 = min(
        W - MOV_W,
        anchor_x + radius,
    )

    y1 = min(
        H - MOV_H,
        anchor_y + radius,
    )

    if (
        x1 < x0
        or y1 < y0
    ):
        return (
            np.nan,
            np.nan,
            np.nan,
        )

    region = raw_scaled[
        y0:y1 + MOV_H + 1,
        x0:x1 + MOV_W + 1,
    ]

    if (
        region.shape[0] < MOV_H
        or region.shape[1] < MOV_W
    ):
        return (
            np.nan,
            np.nan,
            np.nan,
        )

    result = cv2.matchTemplate(
        region,
        movie_frame,
        cv2.TM_CCOEFF_NORMED,
    )

    _, maxval, _, maxloc = (
        cv2.minMaxLoc(
            result
        )
    )

    return (
        float(maxval),
        int(x0 + maxloc[0]),
        int(y0 + maxloc[1]),
    )


                                                                       
              
                                                                       

cache_paths = {}

for alias, path in CANDIDATES.items():
    cache_paths[alias] = (
        cache_candidate(
            alias,
            path,
        )
    )


raw_cache = {
    alias:
        np.load(
            path,
            mmap_mode="r",
        )

    for alias, path
    in cache_paths.items()
}


movies = {
    alias:
        load_movie(
            path
        )

    for alias, path
    in MOVIES.items()
}


                                                                       
               
                                                                       

COARSE_SCALES = np.round(
    np.arange(
        1.0,
        5.7001,
        0.2,
    ),
    2,
)

ANCHORS = [
    0,
    5,
    11,
]

coarse_rows = []


for movie_alias, movie in movies.items():

    print()
    print("#" * 100)
    print(
        f"COARSE SEARCH {movie_alias}"
    )
    print("#" * 100)

    for candidate_alias, cache in raw_cache.items():

        ns, nt, ny, nx = (
            cache.shape
        )

        nstarts = (
            nt - 12 + 1
        )

        for scene in range(ns):

            for start in range(
                nstarts
            ):

                for scale in COARSE_SCALES:

                                   
                                         
                    raw0 = resize_raw(
                        cache[
                            scene,
                            start + ANCHORS[0],
                        ],
                        float(scale),
                    )

                    if raw0 is None:
                        continue

                    score0, x0, y0 = (
                        global_match(
                            raw0,
                            movie[
                                ANCHORS[0]
                            ],
                        )
                    )

                    scores = [
                        score0
                    ]

                    xs = [
                        x0
                    ]

                    ys = [
                        y0
                    ]

                    valid = True

                                        
                                             
                    for k in ANCHORS[1:]:

                        rawk = resize_raw(
                            cache[
                                scene,
                                start + k,
                            ],
                            float(scale),
                        )

                        if rawk is None:
                            valid = False
                            break

                        sc, xx, yy = (
                            local_match(
                                rawk,
                                movie[k],
                                x0,
                                y0,
                            )
                        )

                        if not np.isfinite(
                            sc
                        ):
                            valid = False
                            break

                        scores.append(
                            sc
                        )

                        xs.append(
                            xx
                        )

                        ys.append(
                            yy
                        )

                    if not valid:
                        continue

                    coarse_rows.append({
                        "movie":
                            movie_alias,

                        "candidate":
                            candidate_alias,

                        "scene":
                            scene,

                        "start":
                            start,

                        "scale":
                            float(scale),

                        "anchor_score_median":
                            float(
                                np.median(
                                    scores
                                )
                            ),

                        "anchor_score_min":
                            float(
                                np.min(
                                    scores
                                )
                            ),

                        "anchor_score_0":
                            float(
                                scores[0]
                            ),

                        "anchor_score_5":
                            float(
                                scores[1]
                            ),

                        "anchor_score_11":
                            float(
                                scores[2]
                            ),

                        "anchor_x":
                            x0,

                        "anchor_y":
                            y0,

                        "x_range":
                            int(
                                max(xs)
                                - min(xs)
                            ),

                        "y_range":
                            int(
                                max(ys)
                                - min(ys)
                            ),
                    })


coarse = pd.DataFrame(
    coarse_rows
)

if coarse.empty:
    raise RuntimeError(
        "No valid coarse-match hypotheses."
    )


coarse = coarse.sort_values(
    [
        "movie",
        "anchor_score_median",
        "anchor_score_min",
    ],
    ascending=[
        True,
        False,
        False,
    ],
)


coarse.to_csv(
    RES
    / "46r2_COARSE_BLINDED_MATCHES.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 100)
print("TOP COARSE HYPOTHESES")
print("=" * 100)

for movie_alias in sorted(
    movies
):

    q = (
        coarse[
            coarse[
                "movie"
            ].eq(
                movie_alias
            )
        ]
        .head(12)
    )

    print()
    print(movie_alias)

    print(
        q[
            [
                "candidate",
                "scene",
                "start",
                "scale",
                "anchor_score_median",
                "anchor_score_min",
                "anchor_x",
                "anchor_y",
            ]
        ].to_string(
            index=False
        )
    )


                                                                       
            
 
                                                     
                                                       
                                                                       

TOP_K = 20

refined_rows = []


for movie_alias, movie in movies.items():

    seed = (
        coarse[
            coarse[
                "movie"
            ].eq(
                movie_alias
            )
        ]
        .head(
            TOP_K
        )
    )

    for _, r in seed.iterrows():

        candidate_alias = (
            r[
                "candidate"
            ]
        )

        scene = int(
            r[
                "scene"
            ]
        )

        start = int(
            r[
                "start"
            ]
        )

        coarse_scale = float(
            r[
                "scale"
            ]
        )

        cache = raw_cache[
            candidate_alias
        ]

        scales = np.round(
            np.arange(
                max(
                    1.0,
                    coarse_scale - 0.20,
                ),
                min(
                    5.7,
                    coarse_scale + 0.20,
                )
                + 0.0001,
                0.02,
            ),
            3,
        )

        for scale in scales:

            raw0 = resize_raw(
                cache[
                    scene,
                    start,
                ],
                float(scale),
            )

            if raw0 is None:
                continue

            score0, x0, y0 = (
                global_match(
                    raw0,
                    movie[0],
                )
            )

            frame_scores = []
            frame_x = []
            frame_y = []

            ok = True

            for k in range(12):

                rawk = resize_raw(
                    cache[
                        scene,
                        start + k,
                    ],
                    float(scale),
                )

                if rawk is None:
                    ok = False
                    break

                if k == 0:
                    sc = score0
                    xx = x0
                    yy = y0

                else:
                    sc, xx, yy = (
                        local_match(
                            rawk,
                            movie[k],
                            x0,
                            y0,
                        )
                    )

                if not np.isfinite(
                    sc
                ):
                    ok = False
                    break

                frame_scores.append(
                    float(sc)
                )

                frame_x.append(
                    int(xx)
                )

                frame_y.append(
                    int(yy)
                )

            if not ok:
                continue

            refined_rows.append({
                "movie":
                    movie_alias,

                "candidate":
                    candidate_alias,

                "scene":
                    scene,

                "start":
                    start,

                "scale":
                    float(scale),

                "score_median":
                    float(
                        np.median(
                            frame_scores
                        )
                    ),

                "score_q25":
                    float(
                        np.quantile(
                            frame_scores,
                            0.25,
                        )
                    ),

                "score_min":
                    float(
                        np.min(
                            frame_scores
                        )
                    ),

                "score_max":
                    float(
                        np.max(
                            frame_scores
                        )
                    ),

                "crop_x_median":
                    float(
                        np.median(
                            frame_x
                        )
                    ),

                "crop_y_median":
                    float(
                        np.median(
                            frame_y
                        )
                    ),

                "crop_x_range":
                    int(
                        max(frame_x)
                        - min(frame_x)
                    ),

                "crop_y_range":
                    int(
                        max(frame_y)
                        - min(frame_y)
                    ),

                **{
                    f"frame_{i:02d}_score":
                        frame_scores[i]

                    for i in range(12)
                },
            })


refined = pd.DataFrame(
    refined_rows
)

if refined.empty:
    raise RuntimeError(
        "No valid refined hypotheses."
    )


refined = refined.sort_values(
    [
        "movie",
        "score_median",
        "score_q25",
        "score_min",
    ],
    ascending=[
        True,
        False,
        False,
        False,
    ],
)


refined_path = (
    RES
    / "46r2_REFINED_BLINDED_MATCHES.tsv"
)

refined.to_csv(
    refined_path,
    sep="\t",
    index=False,
)


                                                                       
                 
 
                                      
                                    
                                                                       

summary_rows = []

for movie_alias in sorted(
    movies
):

    q = refined[
        refined[
            "movie"
        ].eq(
            movie_alias
        )
    ]

    for candidate_alias in sorted(
        CANDIDATES
    ):

        z = q[
            q[
                "candidate"
            ].eq(
                candidate_alias
            )
        ]

        if z.empty:
            continue

        best = z.iloc[0]

        summary_rows.append({
            "movie":
                movie_alias,

            "candidate":
                candidate_alias,

            "scene":
                int(
                    best[
                        "scene"
                    ]
                ),

            "start":
                int(
                    best[
                        "start"
                    ]
                ),

            "scale":
                float(
                    best[
                        "scale"
                    ]
                ),

            "score_median":
                float(
                    best[
                        "score_median"
                    ]
                ),

            "score_q25":
                float(
                    best[
                        "score_q25"
                    ]
                ),

            "score_min":
                float(
                    best[
                        "score_min"
                    ]
                ),

            "score_max":
                float(
                    best[
                        "score_max"
                    ]
                ),

            "crop_x_range":
                int(
                    best[
                        "crop_x_range"
                    ]
                ),

            "crop_y_range":
                int(
                    best[
                        "crop_y_range"
                    ]
                ),
        })


summary = pd.DataFrame(
    summary_rows
)

summary_path = (
    RES
    / "46r2_BLINDED_IDENTITY_SUMMARY.tsv"
)

summary.to_csv(
    summary_path,
    sep="\t",
    index=False,
)


print()
print("=" * 100)
print("BLINDED IDENTITY SUMMARY")
print("=" * 100)

print(
    summary.to_string(
        index=False
    )
)


print()
print("=" * 100)
print("TOP 10 REFINED HYPOTHESES PER MOVIE")
print("=" * 100)

for movie_alias in sorted(
    movies
):

    q = (
        refined[
            refined[
                "movie"
            ].eq(
                movie_alias
            )
        ]
        .head(10)
    )

    print()
    print(movie_alias)

    print(
        q[
            [
                "candidate",
                "scene",
                "start",
                "scale",
                "score_median",
                "score_q25",
                "score_min",
                "score_max",
                "crop_x_range",
                "crop_y_range",
            ]
        ].to_string(
            index=False
        )
    )


                                                                       
            
                                                                       

prov = {
    "implementation_freeze_sha256":
        freeze_sha,

    "candidate_aliases": {
        k: str(v)
        for k, v in CANDIDATES.items()
    },

    "movie_aliases": {
        "M1":
            "published_movie_1",

        "M2":
            "published_movie_2",
    },

    "genotype_used_in_scoring":
        False,

    "phenotype_metric_used":
        False,

    "focus_detection_used":
        False,

    "recovery_direction_used":
        False,
}


(
    RES
    / "46r2_PROVENANCE.json"
).write_text(
    json.dumps(
        prov,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)


print()
print("WRITTEN:")
print(
    RES
    / "46r2_COARSE_BLINDED_MATCHES.tsv"
)
print(
    refined_path
)
print(
    summary_path
)
print(
    RES
    / "46r2_PROVENANCE.json"
)
