#!/usr/bin/env python3

from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from skimage.measure import regionprops_table
from skimage.segmentation import find_boundaries
from skimage.transform import resize

from aicspylibczi import CziFile
from cellpose import models


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
    / "segmentation"
    / "190321_WT"
    / "all_scenes_halfscale"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                   
 
                                                             
                                             
                                                              

DIC_CHANNEL = 0
Z_PLANE = 5
SCALE = 0.5

MIN_SIZE_NATIVE = 200
FLOW_THRESHOLD = 0.4
CELLPROB_THRESHOLD = 0.0

SCENES = range(4)
TIMES = range(16)


                                                              
       
                                                              

files = sorted(
    DATA.glob("*.czi")
)

if len(files) != 1:
    raise RuntimeError(
        f"Expected one WT CZI; found {len(files)}"
    )

path = files[0]

czi = CziFile(
    path
)


print("CZI:", path)
print()
print("LOCKED SETTINGS")
print("DIC channel:", DIC_CHANNEL)
print("Z plane:", Z_PLANE)
print("scale:", SCALE)
print("native min_size:", MIN_SIZE_NATIVE)
print("flow_threshold:", FLOW_THRESHOLD)
print("cellprob_threshold:", CELLPROB_THRESHOLD)


                                                              
                       
                                                              

print()
print("Loading cpsam_v2...")

model = models.CellposeModel(
    gpu=False,
    pretrained_model="cpsam_v2",
)


summary_rows = []
object_rows = []


                                                              
         
                                                              

for scene in SCENES:

    scene_dir = (
        OUT
        / f"scene{scene}"
    )

    scene_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    for t in TIMES:

        print()
        print("=" * 100)
        print(
            f"Scene {scene} | T={t} "
            f"| frame {scene * 16 + t + 1}/64"
        )
        print("=" * 100)


        img, _ = czi.read_image(
            S=scene,
            T=t,
            C=DIC_CHANNEL,
            Z=Z_PLANE,
        )

        raw = np.squeeze(
            np.asarray(img)
        ).astype(
            np.float32
        )

        if raw.ndim != 2:
            raise RuntimeError(
                f"Scene {scene} T={t}: "
                f"expected YX; got {raw.shape}"
            )


                                                              
                                                            
                                                              

        lo, hi = np.percentile(
            raw,
            [1, 99],
        )

        norm = np.clip(
            (
                raw
                - lo
            )
            / max(
                float(
                    hi - lo
                ),
                1.0,
            ),
            0,
            1,
        )


        small_shape = (
            int(
                round(
                    raw.shape[0]
                    * SCALE
                )
            ),
            int(
                round(
                    raw.shape[1]
                    * SCALE
                )
            ),
        )


        small = resize(
            norm,
            small_shape,
            order=1,
            preserve_range=True,
            anti_aliasing=True,
        ).astype(
            np.float32
        )


        start = perf_counter()

        masks_small, flows, styles = model.eval(
            small,
            normalize=True,
            diameter=None,
            min_size=int(
                round(
                    MIN_SIZE_NATIVE
                    * SCALE
                    * SCALE
                )
            ),
            flow_threshold=FLOW_THRESHOLD,
            cellprob_threshold=CELLPROB_THRESHOLD,
        )

        elapsed = (
            perf_counter()
            - start
        )


        masks = resize(
            np.asarray(
                masks_small
            ),
            raw.shape,
            order=0,
            preserve_range=True,
            anti_aliasing=False,
        ).astype(
            np.int32
        )


        n_masks = int(
            masks.max()
        )


        np.save(
            scene_dir
            / f"T{t:02d}_masks.npy",
            masks,
        )


        if n_masks > 0:

            props = pd.DataFrame(
                regionprops_table(
                    masks,
                    properties=(
                        "label",
                        "area",
                        "centroid",
                        "eccentricity",
                        "solidity",
                        "bbox",
                    ),
                )
            )

            props.insert(
                0,
                "scene",
                scene,
            )

            props.insert(
                1,
                "time_index",
                t,
            )

            object_rows.append(
                props
            )

            median_area = float(
                props[
                    "area"
                ].median()
            )

            q05_area = float(
                props[
                    "area"
                ].quantile(
                    0.05
                )
            )

            q95_area = float(
                props[
                    "area"
                ].quantile(
                    0.95
                )
            )

        else:

            median_area = np.nan
            q05_area = np.nan
            q95_area = np.nan


        summary_rows.append({
            "scene":
                scene,

            "time_index":
                t,

            "n_masks":
                n_masks,

            "median_area_px":
                median_area,

            "q05_area_px":
                q05_area,

            "q95_area_px":
                q95_area,

            "inference_seconds":
                elapsed,
        })


        print(
            "masks:",
            n_masks,
        )

        print(
            "seconds:",
            f"{elapsed:.1f}",
        )


                                                              
             
                                                              

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUT
    / "SEGMENTATION_SUMMARY.tsv",
    sep="\t",
    index=False,
)


objects = pd.concat(
    object_rows,
    ignore_index=True,
)

objects.to_csv(
    OUT
    / "SEGMENTATION_OBJECTS.tsv",
    sep="\t",
    index=False,
)


                                                              
                              
 
                    
                                                              

for scene in SCENES:

    fig, axes = plt.subplots(
        4,
        4,
        figsize=(16, 15),
        constrained_layout=True,
    )

    axes = axes.ravel()


    for ax, t in zip(
        axes,
        TIMES,
    ):

        img, _ = czi.read_image(
            S=scene,
            T=t,
            C=DIC_CHANNEL,
            Z=Z_PLANE,
        )

        raw = np.squeeze(
            np.asarray(img)
        ).astype(
            np.float32
        )

        lo, hi = np.percentile(
            raw,
            [1, 99],
        )

        norm = np.clip(
            (
                raw
                - lo
            )
            / max(
                float(
                    hi - lo
                ),
                1.0,
            ),
            0,
            1,
        )


        mask = np.load(
            OUT
            / f"scene{scene}"
            / f"T{t:02d}_masks.npy"
        )

        boundaries = find_boundaries(
            mask,
            mode="outer",
        )


        rgb = np.stack(
            [
                norm,
                norm,
                norm,
            ],
            axis=-1,
        )


                           
        rgb[
            boundaries
        ] = 1.0


        n = int(
            mask.max()
        )

        ax.imshow(
            rgb
        )

        ax.set_title(
            f"T={t} | n={n}"
        )

        ax.axis(
            "off"
        )


    fig.suptitle(
        (
            f"190321 WT | Scene {scene} | "
            "LOCKED Cellpose cpsam_v2 | "
            "DIC Z=5 | 0.5x"
        ),
        fontsize=15,
    )


    fig.savefig(
        OUT
        / f"SCENE{scene}_SEGMENTATION_QC.png",
        dpi=150,
    )

    plt.close(
        fig
    )


                                                              
         
                                                              

print()
print("=" * 120)
print("SEGMENTATION SUMMARY")
print("=" * 120)

print(
    summary.to_string(
        index=False
    )
)

print()
print(
    "Total inference hours:",
    f"{summary['inference_seconds'].sum() / 3600:.3f}"
)

print()
print(
    "Output:",
    OUT
)
