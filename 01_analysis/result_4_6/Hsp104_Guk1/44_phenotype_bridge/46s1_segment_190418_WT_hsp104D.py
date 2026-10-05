#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import time

import numpy as np

from cellpose import models
from skimage.measure import regionprops


def summarize(mask):

    ny, nx = mask.shape
    props = regionprops(mask)

    keep = []

    for p in props:

        minr, minc, maxr, maxc = p.bbox

        border = (
            minr <= 0 or
            minc <= 0 or
            maxr >= ny or
            maxc >= nx
        )

        plausible = 300 <= p.area <= 10000

        if (not border) and plausible:
            keep.append(p)

    areas = np.array(
        [p.area for p in keep],
        dtype=float,
    )

    if len(areas):
        med = float(np.median(areas))
        q10 = float(np.percentile(areas, 10))
        q90 = float(np.percentile(areas, 90))
    else:
        med = q10 = q90 = np.nan

    return props, keep, areas, med, q10, q90


def filter_and_relabel(mask, keep):

    out = np.zeros(
        mask.shape,
        dtype=np.uint16,
    )

    for new_id, p in enumerate(keep, 1):
        out[mask == p.label] = new_id

    return out


def run_dataset(
    model,
    model_name,
    genotype,
    cache,
    dest,
):

    dic = np.load(
        cache / "dic_midplane.npy",
        mmap_mode="r",
    )

    ns, nt, ny, nx = dic.shape

    masks_dir = (
        dest /
        model_name /
        genotype /
        "masks"
    )
    masks_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    table_dir = (
        dest /
        model_name /
        genotype
    )
    table_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    table = table_dir / "segmentation.tsv"

    fields = [
        "model",
        "genotype",
        "scene",
        "time_index",
        "n_raw",
        "n_kept",
        "median_area_px",
        "q10_area_px",
        "q90_area_px",
        "foreground_fraction",
        "seconds",
    ]

                                                   
    rows = []

    for s in range(ns):
        for t in range(nt):

            maskfile = (
                masks_dir /
                f"{genotype}_S{s}_T{t:02d}.npy"
            )

            start = time.time()

            if maskfile.exists():

                filtered = np.load(maskfile)

                props = regionprops(filtered)

                areas = np.array(
                    [p.area for p in props],
                    dtype=float,
                )

                n_raw = len(props)
                n_kept = len(props)

                if len(areas):
                    med = float(np.median(areas))
                    q10 = float(np.percentile(areas, 10))
                    q90 = float(np.percentile(areas, 90))
                else:
                    med = q10 = q90 = np.nan

                foreground = float(
                    np.mean(filtered > 0)
                )

                elapsed = 0.0

                print(
                    f"SKIP {model_name} {genotype} "
                    f"S={s} T={t:02d} "
                    f"cells={n_kept}",
                    flush=True,
                )

            else:

                img = np.asarray(dic[s, t])

                mask, flows, styles = model.eval(
                    img,
                    channels=[0, 0],
                    diameter=60,
                    flow_threshold=0.4,
                    cellprob_threshold=0.0,
                    min_size=100,
                )

                mask = np.asarray(mask)

                if mask.shape != (ny, nx):
                    raise RuntimeError(
                        f"Unexpected mask shape "
                        f"{mask.shape}"
                    )

                props, keep, areas, med, q10, q90 = \
                    summarize(mask)

                n_raw = len(props)
                n_kept = len(keep)

                filtered = filter_and_relabel(
                    mask,
                    keep,
                )

                np.save(
                    maskfile,
                    filtered,
                )

                foreground = float(
                    np.mean(filtered > 0)
                )

                elapsed = time.time() - start

                print(
                    f"DONE {model_name:14s} "
                    f"{genotype:9s} "
                    f"S={s} T={t:02d} "
                    f"raw={n_raw:4d} "
                    f"kept={n_kept:4d} "
                    f"area={med:8.1f} "
                    f"fg={foreground:.4f} "
                    f"sec={elapsed:6.2f}",
                    flush=True,
                )

            rows.append({
                "model": model_name,
                "genotype": genotype,
                "scene": s,
                "time_index": t,
                "n_raw": n_raw,
                "n_kept": n_kept,
                "median_area_px": med,
                "q10_area_px": q10,
                "q90_area_px": q90,
                "foreground_fraction": foreground,
                "seconds": elapsed,
            })

    with table.open("w", newline="") as f:

        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
        )

        w.writeheader()
        w.writerows(rows)

    complete = table_dir / "COMPLETE.txt"

    complete.write_text(
        f"model={model_name}\n"
        f"genotype={genotype}\n"
        f"scenes={ns}\n"
        f"timepoints={nt}\n"
        f"frames={ns*nt}\n"
    )


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--out",
        required=True,
    )

    args = ap.parse_args()

    root = Path(
        "/"
    )

    source = root / "44_phenotype_bridge"

    dest = Path(args.out)
    dest.mkdir(
        parents=True,
        exist_ok=True,
    )

    datasets = [
        (
            "190418_WT_hsp104D",
            source / "projection_cache_190418_WT_hsp104D",
        ),
    ]

    model_names = [
        "yeast_BF_cp3",
        "yeast_PhC_cp3",
    ]

    global_start = time.time()

    for model_name in model_names:

        print()
        print("=" * 80)
        print("LOAD MODEL:", model_name)
        print("=" * 80)

        start = time.time()

        model = models.CellposeModel(
            gpu=False,
            model_type=model_name,
        )

        print(
            f"model_load_seconds="
            f"{time.time()-start:.2f}",
            flush=True,
        )

        for genotype, cache in datasets:

            print()
            print("=" * 80)
            print(model_name, genotype)
            print("=" * 80)

            run_dataset(
                model=model,
                model_name=model_name,
                genotype=genotype,
                cache=cache,
                dest=dest,
            )

    (dest / "44d1_COMPLETE.txt").write_text(
        "dual yeast segmentation complete\n"
    )

    print()
    print("=" * 80)
    print("ALL COMPLETE")
    print("=" * 80)
    print(
        f"total_seconds="
        f"{time.time()-global_start:.2f}"
    )


if __name__ == "__main__":
    main()
