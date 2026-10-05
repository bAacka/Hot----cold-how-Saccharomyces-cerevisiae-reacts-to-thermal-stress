#!/usr/bin/env python3

from pathlib import Path
import hashlib
import math

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle


                                                              
       
                                                              

ROOT = Path("/")

TRAJ = ROOT / "14_response_programs/response_class_trajectories_long.tsv.gz"
CLASS = ROOT / "14_response_programs/response_classification.tsv"
SUMMARY = ROOT / "14_response_programs/response_class_summary.tsv"

OUTROOT = ROOT / "54_figures"
CURATED = OUTROOT / "02_curated_inputs" / "Fig4_1"
FINAL = OUTROOT / "05_final"
LEGENDS = OUTROOT / "06_legends"

for d in [CURATED, FINAL, LEGENDS]:
    d.mkdir(parents=True, exist_ok=True)


                                                              
                             
                                                              

N_EXPECTED = 347
TIMES = list(range(0, 29, 2))

CLASS_ORDER = [
    "Heat",
    "Bidirectional",
    "Temp_Independent",
    "Cold",
]

EXPECTED_COUNTS = {
    "Heat": 163,
    "Bidirectional": 124,
    "Temp_Independent": 56,
    "Cold": 4,
}

CLASS_COLORS = {
    "Heat": "#D95F4B",
    "Bidirectional": "#8C6BB1",
    "Temp_Independent": "#4DAF8A",
    "Cold": "#4C78A8",
}

HEAT_COLOR = "#C84B31"
COLD_COLOR = "#3478A8"
ZERO_COLOR = "#666666"
MISSING_COLOR = "#E8E8E8"

REPRESENTATIVES = [
    ("RGA2_S733", "Rga2-S733"),
    ("HSP42_S223", "Hsp42-S223"),
    ("LEU1_T345", "Leu1-T345"),
    ("NUP60_S10", "Nup60-S10"),
    ("MSN2_S633", "Msn2-S633"),
]


                                                              
         
                                                              

def fail(msg):
    raise RuntimeError(msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def round_up_half(x):
    return math.ceil(float(x) * 2.0) / 2.0


                                                              
      
                                                              

print("===== READ INPUTS =====")

traj = pd.read_csv(TRAJ, sep="\t", low_memory=False)
cls = pd.read_csv(CLASS, sep="\t", low_memory=False)
summary = pd.read_csv(SUMMARY, sep="\t", low_memory=False)

print("trajectory rows :", len(traj))
print("classification  :", len(cls))
print("summary rows    :", len(summary))


                                                              
                   
                                                              

required_traj = {
    "site_id",
    "gene",
    "psite_group",
    "response_class",
    "response_semantics",
    "condition",
    "time_min",
    "log2FC",
}

required_cls = {
    "site_id",
    "gene_final",
    "psite_only",
    "response_class_final",
    "response_semantics",
    "heat_signed_auc",
    "cold_signed_auc",
    "heat_cold_paired_correlation",
}

for name, df, required in [
    ("trajectory", traj, required_traj),
    ("classification", cls, required_cls),
]:
    missing = required - set(df.columns)
    if missing:
        fail(f"{name}: missing columns: {sorted(missing)}")


if len(cls) != N_EXPECTED:
    fail(f"Expected 347 classification rows; found {len(cls)}")

if cls["site_id"].nunique() != N_EXPECTED:
    fail("Classification site_id values are not unique.")

expected_long = N_EXPECTED * 2 * len(TIMES)

if len(traj) != expected_long:
    fail(
        f"Expected {expected_long} trajectory rows; "
        f"found {len(traj)}"
    )

if traj["site_id"].nunique() != N_EXPECTED:
    fail("Trajectory table does not contain exactly 347 site_ids.")

if set(traj["site_id"]) != set(cls["site_id"]):
    fail("Trajectory and classification site_id sets differ.")

conditions = sorted(
    traj["condition"]
    .astype(str)
    .str.lower()
    .unique()
)

if conditions != ["cold", "heat"]:
    fail(f"Unexpected conditions: {conditions}")

times_found = sorted(
    pd.to_numeric(
        traj["time_min"],
        errors="raise"
    )
    .astype(int)
    .unique()
)

if times_found != TIMES:
    fail(f"Unexpected time points: {times_found}")

duplicates = traj.duplicated(
    ["site_id", "condition", "time_min"]
).sum()

if duplicates:
    fail(f"Duplicate long-table observations: {duplicates}")


                     
check = (
    traj[
        ["site_id", "response_class"]
    ]
    .drop_duplicates()
    .merge(
        cls[
            ["site_id", "response_class_final"]
        ],
        on="site_id",
        validate="one_to_one",
    )
)

bad = check[
    check["response_class"]
    != check["response_class_final"]
]

if len(bad):
    print(bad.to_string(index=False))
    fail("response_class differs from response_class_final.")


observed_counts = (
    cls["response_class_final"]
    .value_counts()
    .to_dict()
)

if observed_counts != EXPECTED_COUNTS:
    fail(
        f"Class counts differ.\n"
        f"Expected: {EXPECTED_COUNTS}\n"
        f"Observed: {observed_counts}"
    )

print("Validation: PASS")


                                                              
                               
                                                              

meta = cls[
    [
        "site_id",
        "gene_final",
        "psite_only",
        "response_class_final",
        "response_semantics",
        "heat_signed_auc",
        "cold_signed_auc",
        "heat_cold_paired_correlation",
    ]
].copy()

meta["delta_auc"] = (
    meta["heat_signed_auc"]
    - meta["cold_signed_auc"]
)

class_rank = {
    c: i for i, c in enumerate(CLASS_ORDER)
}

meta["_class_rank"] = (
    meta["response_class_final"]
    .map(class_rank)
)

meta = (
    meta
    .sort_values(
        ["_class_rank", "delta_auc", "site_id"],
        ascending=[True, True, True],
        na_position="last",
    )
    .reset_index(drop=True)
)

row_order = meta["site_id"].tolist()


                                                              
                             
                                                              

long = traj.copy()

long["condition"] = (
    long["condition"]
    .astype(str)
    .str.lower()
)

long["time_min"] = pd.to_numeric(
    long["time_min"],
    errors="raise",
).astype(int)

long["log2FC"] = pd.to_numeric(
    long["log2FC"],
    errors="coerce",
)


def make_matrix(condition):
    x = (
        long[
            long["condition"] == condition
        ]
        .pivot(
            index="site_id",
            columns="time_min",
            values="log2FC",
        )
        .reindex(
            index=row_order,
            columns=TIMES,
        )
    )

    if x.shape != (347, 15):
        fail(
            f"{condition}: unexpected matrix shape {x.shape}"
        )

    return x


heat_df = make_matrix("heat")
cold_df = make_matrix("cold")

heat_matrix = heat_df.to_numpy(float)
cold_matrix = cold_df.to_numpy(float)

                                                
gap = np.full((N_EXPECTED, 1), np.nan)

heatmap_matrix = np.concatenate(
    [
        heat_matrix,
        gap,
        cold_matrix,
    ],
    axis=1,
)


                                                              
               
                                                              

finite_heatmap = np.abs(
    heatmap_matrix[
        np.isfinite(heatmap_matrix)
    ]
)

vlim = round_up_half(
    np.quantile(
        finite_heatmap,
        0.99,
    )
)

if vlim <= 0:
    vlim = 1.0

print(f"Heatmap colour range: ±{vlim:.1f}")


                                                              
                  
                                                              

boundaries = {}

start = 0

for c in CLASS_ORDER:
    end = start + EXPECTED_COUNTS[c]
    boundaries[c] = (start, end)
    start = end

if start != N_EXPECTED:
    fail("Class boundaries do not sum to 347.")


                                                              
                              
                                                              

missing_reps = [
    x
    for x, _ in REPRESENTATIVES
    if x not in set(meta["site_id"])
]

if missing_reps:
    fail(
        "Missing representative sites: "
        + ", ".join(missing_reps)
    )

rep_abs = []

for site_id, _ in REPRESENTATIVES:
    x = np.concatenate(
        [
            heat_df.loc[site_id].to_numpy(float),
            cold_df.loc[site_id].to_numpy(float),
        ]
    )
    rep_abs.extend(
        np.abs(
            x[np.isfinite(x)]
        )
    )

rep_ylim = round_up_half(
    max(rep_abs) * 1.08
)

if rep_ylim <= 0:
    rep_ylim = 1.0


                                                              
                       
                                                              

auc_values = np.concatenate(
    [
        meta["heat_signed_auc"].to_numpy(float),
        meta["cold_signed_auc"].to_numpy(float),
    ]
)

auc_values = auc_values[
    np.isfinite(auc_values)
]

auc_lim = np.max(
    np.abs(auc_values)
) * 1.08

auc_lim = math.ceil(
    auc_lim / 10.0
) * 10.0

if auc_lim <= 0:
    auc_lim = 10


                                                              
            
                                                              

manifest = pd.DataFrame(
    [
        {
            "role": "dynamic_trajectories",
            "path": str(TRAJ),
            "sha256": sha256(TRAJ),
            "size_bytes": TRAJ.stat().st_size,
        },
        {
            "role": "final_response_classification",
            "path": str(CLASS),
            "sha256": sha256(CLASS),
            "size_bytes": CLASS.stat().st_size,
        },
        {
            "role": "response_class_summary",
            "path": str(SUMMARY),
            "sha256": sha256(SUMMARY),
            "size_bytes": SUMMARY.stat().st_size,
        },
    ]
)

manifest.to_csv(
    CURATED / "Fig4_1_source_manifest.tsv",
    sep="\t",
    index=False,
)

meta.drop(
    columns=["_class_rank"]
).to_csv(
    CURATED / "Fig4_1_row_order.tsv",
    sep="\t",
    index=False,
)


                                                              
               
                                                              

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.2,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "legend.fontsize": 8.0,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


                                                              
                 
                                                              

fig = plt.figure(
    figsize=(16.0, 11.5),
    layout="constrained",
)

outer = fig.add_gridspec(
    nrows=3,
    ncols=1,
    height_ratios=[
        2.25,
        5.2,
        2.55,
    ],
)


                                                              
                
                                                              

top = outer[0].subgridspec(
    1,
    2,
    width_ratios=[
        1.0,
        2.45,
    ],
    wspace=0.18,
)


                                                              
                  
                                                              

axA = fig.add_subplot(
    top[0, 0]
)

summary_ix = (
    summary
    .set_index("response_class_final")
)

n_groups = [
    int(
        summary_ix.loc[c, "n_groups"]
    )
    for c in CLASS_ORDER
]

n_genes = [
    int(
        summary_ix.loc[c, "n_genes"]
    )
    for c in CLASS_ORDER
]

y = np.arange(
    len(CLASS_ORDER)
)

bars = axA.barh(
    y,
    n_groups,
    height=0.62,
    color=[
        CLASS_COLORS[c]
        for c in CLASS_ORDER
    ],
    edgecolor="none",
)

axA.set_yticks(y)
axA.set_yticklabels(CLASS_ORDER)
axA.invert_yaxis()

axA.set_xlabel(
    "Брой фосфопротеомни групи"
)

axA.set_xlim(
    0,
    max(n_groups) * 1.38,
)

for bar, n, ng in zip(
    bars,
    n_groups,
    n_genes,
):
    axA.text(
        bar.get_width()
        + max(n_groups) * 0.025,
        bar.get_y()
        + bar.get_height() / 2,
        f"{n} ({ng} гена)",
        ha="left",
        va="center",
        fontsize=8.5,
    )

axA.text(
    0.98,
    0.04,
    "Общо: 347 групи",
    transform=axA.transAxes,
    ha="right",
    va="bottom",
    fontsize=8.7,
)

axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

axA.text(
    -0.16,
    1.10,
    "A",
    transform=axA.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                          
                                                              

axC = fig.add_subplot(
    top[0, 1]
)

axC.axhline(
    0,
    color="#888888",
    lw=0.8,
    zorder=0,
)

axC.axvline(
    0,
    color="#888888",
    lw=0.8,
    zorder=0,
)

axC.plot(
    [-auc_lim, auc_lim],
    [-auc_lim, auc_lim],
    ls="--",
    lw=0.9,
    color="#777777",
    alpha=0.7,
    zorder=0,
)

for c in CLASS_ORDER:
    sub = meta[
        meta["response_class_final"] == c
    ]

    axC.scatter(
        sub["heat_signed_auc"],
        sub["cold_signed_auc"],
        s=24 if c != "Cold" else 34,
        alpha=0.78,
        color=CLASS_COLORS[c],
        edgecolors="white",
        linewidths=0.35,
        label=f"{c} (n={len(sub)})",
        zorder=2,
    )

axC.set_xlim(
    -auc_lim,
    auc_lim,
)

axC.set_ylim(
    -auc_lim,
    auc_lim,
)

axC.set_aspect(
    "equal",
    adjustable="box",
)

axC.set_xlabel(
    "Знакова AUC при топлинен стрес"
)

axC.set_ylabel(
    "Знакова AUC при студов стрес"
)

axC.legend(
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(1.01, 1.0),
    borderaxespad=0,
    handletextpad=0.4,
)

axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

axC.text(
    -0.09,
    1.10,
    "C",
    transform=axC.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                    
                                                              

middle = outer[1].subgridspec(
    1,
    2,
    width_ratios=[
        1.75,
        12.0,
    ],
    wspace=0.02,
)

axClass = fig.add_subplot(
    middle[0, 0]
)

axB = fig.add_subplot(
    middle[0, 1]
)


                                                              
                            
                                                              

axClass.set_xlim(
    0,
    1,
)

axClass.set_ylim(
    N_EXPECTED,
    0,
)

axClass.set_xticks([])
axClass.set_yticks([])

for spine in axClass.spines.values():
    spine.set_visible(False)

for c in CLASS_ORDER:

    start, end = boundaries[c]
    height = end - start
    mid = (start + end) / 2

    axClass.add_patch(
        Rectangle(
            (0.82, start),
            0.16,
            height,
            facecolor=CLASS_COLORS[c],
            edgecolor="none",
        )
    )

    axClass.text(
        0.76,
        mid,
        c,
        ha="right",
        va="center",
        fontsize=8.3,
    )


                                                              
         
                                                              

cmap = plt.get_cmap(
    "RdBu_r"
).copy()

cmap.set_bad(
    MISSING_COLOR
)

im = axB.imshow(
    heatmap_matrix,
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=-vlim,
    vmax=vlim,
)

                       
separator_x = len(TIMES) - 0.5

axB.axvline(
    separator_x + 0.5,
    color="white",
    lw=7,
    zorder=3,
)


                  
for c in CLASS_ORDER[:-1]:

    end = boundaries[c][1]

    axB.axhline(
        end - 0.5,
        color="black",
        lw=0.65,
        alpha=0.65,
    )


            
tick_times = [
    0, 4, 8, 12, 16, 20, 24, 28
]

heat_ticks = [
    TIMES.index(t)
    for t in tick_times
]

cold_offset = len(TIMES) + 1

cold_ticks = [
    cold_offset + TIMES.index(t)
    for t in tick_times
]

axB.set_xticks(
    heat_ticks + cold_ticks
)

axB.set_xticklabels(
    [str(t) for t in tick_times]
    + [str(t) for t in tick_times]
)

axB.set_yticks([])

axB.set_xlabel(
    "Време след температурното въздействие (min)"
)

heat_center = (
    heat_ticks[0]
    + heat_ticks[-1]
) / 2

cold_center = (
    cold_ticks[0]
    + cold_ticks[-1]
) / 2

axB.text(
    heat_center,
    -13,
    "Топлинен стрес",
    ha="center",
    va="bottom",
    color=HEAT_COLOR,
    fontsize=10.5,
    fontweight="bold",
)

axB.text(
    cold_center,
    -13,
    "Студов стрес",
    ha="center",
    va="bottom",
    color=COLD_COLOR,
    fontsize=10.5,
    fontweight="bold",
)

cbar = fig.colorbar(
    im,
    ax=axB,
    fraction=0.018,
    pad=0.012,
)

cbar.set_label(
    "log$_2$FC"
)

axClass.text(
    -0.03,
    1.055,
    "B",
    transform=axClass.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                 
                                                              

bottom = outer[2].subgridspec(
    1,
    5,
    wspace=0.22,
)

axesD = []

for i, (site_id, label) in enumerate(
    REPRESENTATIVES
):

    ax = fig.add_subplot(
        bottom[0, i]
    )

    axesD.append(ax)

    heat_y = (
        heat_df
        .loc[site_id]
        .to_numpy(float)
    )

    cold_y = (
        cold_df
        .loc[site_id]
        .to_numpy(float)
    )

    site_class = (
        meta
        .set_index("site_id")
        .loc[
            site_id,
            "response_class_final"
        ]
    )

    ax.axhline(
        0,
        color=ZERO_COLOR,
        lw=0.7,
        alpha=0.65,
    )

    ax.plot(
        TIMES,
        heat_y,
        color=HEAT_COLOR,
        marker="o",
        markersize=3.0,
        lw=1.5,
        label="Топлинен стрес",
    )

    ax.plot(
        TIMES,
        cold_y,
        color=COLD_COLOR,
        marker="o",
        markersize=3.0,
        lw=1.5,
        label="Студов стрес",
    )

    ax.set_xlim(
        0,
        28,
    )

    ax.set_ylim(
        -rep_ylim,
        rep_ylim,
    )

    ax.set_xticks(
        [0, 8, 16, 24, 28]
    )

    ax.set_xlabel(
        "Време (min)"
    )

    ax.set_title(
        f"{label}\n{site_class}",
        pad=7,
    )

    if i == 0:
        ax.set_ylabel(
            "log$_2$FC"
        )

        ax.text(
            -0.18,
            1.14,
            "D",
            transform=ax.transAxes,
            fontsize=15,
            fontweight="bold",
            va="top",
        )
    else:
        ax.tick_params(
            axis="y",
            labelleft=False,
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


                                             
axesD[-1].legend(
    frameon=False,
    loc="upper right",
    bbox_to_anchor=(1.0, -0.23),
    ncol=1,
)


                                                              
        
                                                              

OUTBASE = (
    FINAL
    / "Fig4_1_dynamic_phosphoproteomic_response_v2"
)

fig.savefig(
    str(OUTBASE) + ".png",
    dpi=600,
)

fig.savefig(
    str(OUTBASE) + ".pdf",
)

fig.savefig(
    str(OUTBASE) + ".svg",
)

plt.close(fig)


                                                              
        
                                                              

legend = f"""Фигура 4.1. Динамичен фосфопротеомен отговор при топлинен и студов стрес.

(A) Разпределение на 347-те динамични фосфопротеомни групи в четирите окончателни класа на температурния отговор: Heat (n=163), Bidirectional (n=124), Temp_Independent (n=56) и Cold (n=4). В скоби е показан броят на представените гени.

(B) Сдвоена топлинна карта на фосфорилационните траектории през първите 28 min след топлинен и студов стрес. Всеки ред представя една и съща експериментална фосфопротеомна група при двете условия. Групите са подредени първо според окончателния биологичен клас и след това според разликата между знаковата AUC при топлинен и студов стрес. Цветовата скала показва log₂FC и е ограничена за визуализация до 99-ия персентил на абсолютните измерени стойности (±{vlim:.1f} log₂FC). Сивите клетки обозначават липсващи измервания.

(C) Съпоставка между знаковата AUC при топлинен и студов стрес за всички 347 динамични фосфопротеомни групи. Цветът показва окончателния клас на температурния отговор. Прекъснатата диагонална линия обозначава равен интегрален отговор при двете температурни условия.

(D) Представителни високоамплитудни времеви профили за Rga2-S733, Hsp42-S223, Leu1-T345, Nup60-S10 и Msn2-S633.

Класовете на температурния отговор са базирани на окончателната биологична класификация и не са еквивалентни на първоначалните кинетични ClusterID означения.
"""

legend_path = (
    LEGENDS
    / "Fig4_1_legend_bg_v2.txt"
)

legend_path.write_text(
    legend,
    encoding="utf-8",
)


print()
print("===== FIGURE 4.1 V2 COMPLETE =====")
print("Class counts:", EXPECTED_COUNTS)
print(f"Heatmap scale: ±{vlim:.1f} log2FC")
print(f"AUC range: ±{auc_lim:.0f}")
print()
print(str(OUTBASE) + ".png")
print(str(OUTBASE) + ".pdf")
print(str(OUTBASE) + ".svg")
print(legend_path)
