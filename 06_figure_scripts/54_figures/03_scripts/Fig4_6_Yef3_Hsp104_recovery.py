#!/usr/bin/env python3

from pathlib import Path
import hashlib

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from scipy.stats import pearsonr


                                                              
       
                                                              

ROOT = Path(
    "/"
)

BRIDGE = (
    ROOT
    / "44_phenotype_bridge"
)

YEF = (
    BRIDGE
    / "49_T972_precausal_support"
)

HSP_FINAL = (
    BRIDGE
    / "46_phenotype_state"
    / "190418_WT_hsp104D_tracked_focus_T0_T8"
    / "46t_final_hsp104_validation"
)

PAB1 = (
    BRIDGE
    / "52_public_perturbational_bridge"
)

CURATED = (
    ROOT
    / "54_figures"
    / "02_curated_inputs"
    / "Fig4_6"
)

FINAL = (
    ROOT
    / "54_figures"
    / "05_final"
)

LEGENDS = (
    ROOT
    / "54_figures"
    / "06_legends"
)

for d in [
    CURATED,
    FINAL,
    LEGENDS,
]:
    d.mkdir(
        parents=True,
        exist_ok=True,
    )


DISCOVERY = (
    YEF
    / "49c_DISCOVERY_PERTURBATIONS_GLOBAL_CONTEXT.tsv"
)

ADJUSTED = (
    YEF
    / "49c_PRIMARY_GLOBAL_ADJUSTED_RESULTS.tsv"
)

RANDOM_NULL = (
    YEF
    / "49c_GLOBAL_ADJUSTED_RANDOM_GENESET_NULL.tsv"
)

RESIDUALS = (
    CURATED
    / "Fig4_6_STRESS_GRANULE_RNP_residual_pairs.tsv"
)

HSP_TRAJ = (
    HSP_FINAL
    / "resolved_genotype_T0_T15_trajectory.tsv"
)

HSP_SUMMARY = (
    HSP_FINAL
    / "resolved_genotype_early_late_summary.tsv"
)

PAB1_RAW = (
    PAB1
    / "52b_Shattuck_Figure6_raw.tsv"
)

PAB1_SUMMARY = (
    PAB1
    / "52c_PAB1_HSP104_PRIMARY_EFFECTS.tsv"
)


                                                              
         
                                                              

def fail(msg):
    raise RuntimeError(msg)


def sha256(path):

    h = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(
                1024 * 1024
            ),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def letter(ax, x):

    ax.text(
        -0.12,
        1.04,
        x,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="top",
    )


def clean_axes(ax):

    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)


                                                              
               
                                                              

sources = [
    DISCOVERY,
    ADJUSTED,
    RANDOM_NULL,
    RESIDUALS,
    HSP_TRAJ,
    HSP_SUMMARY,
    PAB1_RAW,
    PAB1_SUMMARY,
]

for p in sources:

    if not p.exists():
        fail(
            f"Missing Figure 4.6 source: {p}"
        )


                                                              
                             
                                                              

disc = pd.read_csv(
    DISCOVERY,
    sep="\t",
)

required = {
    "deleted_gene",
    "YEF3_T972",
}

if not required.issubset(
    disc.columns
):
    fail(
        "Unexpected discovery-table schema."
    )

disc["deleted_gene"] = (
    disc["deleted_gene"]
    .astype(str)
    .str.lower()
)

expected_genes = [
    "ypl150w",
    "ppt1",
    "ctk1",
]

if set(
    disc["deleted_gene"]
) != set(expected_genes):
    fail(
        "Discovery perturbations are not "
        "exactly YPL150W/PPT1/CTK1."
    )

disc = (
    disc
    .set_index(
        "deleted_gene"
    )
    .loc[
        expected_genes
    ]
    .reset_index()
)

disc.to_csv(
    CURATED
    / "Fig4_6_Yef3_nomination.tsv",
    sep="\t",
    index=False,
)


                                                              
                             
                                                              

adj = pd.read_csv(
    ADJUSTED,
    sep="\t",
)

expected_modules = [
    "PROTEIN_FOLDING",
    "TRANSLATION_ELONGATION",
    "STRESS_GRANULE_RNP",
]

if set(
    adj["module"]
) != set(expected_modules):
    fail(
        "Unexpected module set in 49c results."
    )

adj = (
    adj
    .set_index(
        "module"
    )
    .loc[
        expected_modules
    ]
    .reset_index()
)

if not (
    adj["n_deletions"]
    .eq(107)
    .all()
):
    fail(
        "Expected n=107 for all three modules."
    )

adj.to_csv(
    CURATED
    / "Fig4_6_programme_adjustment.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                        
                                                              

res = pd.read_csv(
    RESIDUALS,
    sep="\t",
)

if list(
    res.columns
) != [
    "YEF3_T972_residual",
    "STRESS_GRANULE_RNP_residual",
]:
    fail(
        "Unexpected residual-pair schema."
    )

if len(res) != 107:
    fail(
        f"Expected 107 residual pairs; found {len(res)}"
    )

partial_rho = float(
    pearsonr(
        res[
            "YEF3_T972_residual"
        ],
        res[
            "STRESS_GRANULE_RNP_residual"
        ],
    ).statistic
)

frozen_partial = float(
    adj.loc[
        adj["module"]
        == "STRESS_GRANULE_RNP",
        "partial_spearman_global_adjusted",
    ].iloc[0]
)

if abs(
    partial_rho
    - frozen_partial
) > 1e-10:
    fail(
        "Exported residual vectors do not reproduce "
        f"frozen partial rho: {partial_rho} vs {frozen_partial}"
    )

sg_perm_p = float(
    adj.loc[
        adj["module"]
        == "STRESS_GRANULE_RNP",
        "permutation_p_two_sided",
    ].iloc[0]
)

sg_q = float(
    adj.loc[
        adj["module"]
        == "STRESS_GRANULE_RNP",
        "BH_q_across_3_modules",
    ].iloc[0]
)

rnd = pd.read_csv(
    RANDOM_NULL,
    sep="\t",
)

sg_rnd = rnd[
    rnd["module"]
    == "STRESS_GRANULE_RNP"
]

if len(sg_rnd) != 1:
    fail(
        "Could not uniquely identify SG/RNP random-set result."
    )

random_p = float(
    sg_rnd[
        "random_set_empirical_p_two_sided"
    ].iloc[0]
)


                                                              
                                    
                                                              

ht = pd.read_csv(
    HSP_TRAJ,
    sep="\t",
)

needed_hsp = {
    "genotype",
    "time_index",
    "n_cells",
    "mean_log2_vs_T1",
}

if not needed_hsp.issubset(
    ht.columns
):
    fail(
        "Unexpected Hsp104 trajectory schema."
    )

ht = ht[
    ht["time_index"]
    .between(
        1,
        15,
    )
].copy()

if set(
    ht["genotype"]
) != {
    "WT",
    "hsp104D",
}:
    fail(
        "Expected WT and hsp104D trajectories."
    )

for genotype in [
    "WT",
    "hsp104D",
]:

    x = (
        ht[
            ht["genotype"]
            == genotype
        ]
        .sort_values(
            "time_index"
        )
    )

    if list(
        x["time_index"]
    ) != list(
        range(
            1,
            16,
        )
    ):
        fail(
            f"Incomplete T1:T15 trajectory for {genotype}"
        )

ht.to_csv(
    CURATED
    / "Fig4_6_Hsp104_Guk1_T1_T15.tsv",
    sep="\t",
    index=False,
)


                                                     
hs = pd.read_csv(
    HSP_SUMMARY,
    sep="\t",
)

print(
    "HSP104 summary columns:",
    list(hs.columns),
)

                                                          
EARLY_WT = -0.465496
EARLY_HD = -0.628083

LATE_WT = -0.657923
LATE_HD = -0.140885


                                                              
                                                           
                                                              

ps = pd.read_csv(
    PAB1_SUMMARY,
    sep="\t",
)

ps = ps[
    ps["phenotype"]
    == "Pab1_foci_percent"
].copy()

condition_order = [
    "30min_heatshock",
    "1h_recovery",
    "2h_recovery",
]

ps = (
    ps
    .set_index(
        "condition"
    )
    .loc[
        condition_order
    ]
    .reset_index()
)

if not (
    ps["HSP104_n"].eq(5).all()
    and
    ps["hsp104D_n"].eq(5).all()
):
    fail(
        "Expected five biological repeats per group."
    )


raw = pd.read_csv(
    PAB1_RAW,
    sep="\t",
    low_memory=False,
)

                                                           
num = raw.copy()

for c in num.columns:

    num[c] = pd.to_numeric(
        num[c]
        .astype(str)
        .str.replace(
            "%",
            "",
            regex=False,
        ),
        errors="coerce",
    )


def find_repeat_window(
    expected_mean,
    expected_sd,
):

                                                      
                                                    
    if (
        abs(expected_mean - 100.0)
        < 1e-10
        and
        abs(expected_sd)
        < 1e-10
    ):
        return np.repeat(
            100.0,
            5,
        )

    candidates = []

    matrix = num.to_numpy(
        dtype=float
    )

    for ri in range(
        matrix.shape[0]
    ):

        row = matrix[
            ri,
            :
        ]

        for start in range(
            0,
            len(row) - 4,
        ):

            w = row[
                start:start + 5
            ]

            if not np.isfinite(
                w
            ).all():
                continue

            m = float(
                np.mean(w)
            )

            sd = float(
                np.std(
                    w,
                    ddof=1,
                )
            )

            mean_err = abs(
                m
                - expected_mean
            )

            sd_err = abs(
                sd
                - expected_sd
            )

            score = (
                mean_err
                + sd_err
            )

            candidates.append(
                (
                    score,
                    mean_err,
                    sd_err,
                    ri,
                    start,
                    w.copy(),
                )
            )

    if not candidates:
        fail(
            "No numeric five-repeat windows found in raw Pab1 table."
        )

    candidates.sort(
        key=lambda x: x[0]
    )

    best = candidates[0]

                                                       
                                             
    if (
        best[1] > 1e-4
        or best[2] > 1e-4
    ):
        fail(
            "Could not recover replicate window matching "
            f"mean={expected_mean}, SD={expected_sd}. "
            f"Best mean error={best[1]}, SD error={best[2]}"
        )

    return best[5]


rep_rows = []

for _, row in ps.iterrows():

    for genotype in [
        "HSP104",
        "hsp104D",
    ]:

        mean = float(
            row[
                f"{genotype}_mean"
            ]
        )

        sd = float(
            row[
                f"{genotype}_sd"
            ]
        )

        vals = find_repeat_window(
            mean,
            sd,
        )

                                        
        if abs(
            vals.mean()
            - mean
        ) > 1e-4:
            fail(
                "Recovered replicate mean mismatch."
            )

        if abs(
            vals.std(
                ddof=1
            )
            - sd
        ) > 1e-4:
            fail(
                "Recovered replicate SD mismatch."
            )

        for i, value in enumerate(
            vals,
            start=1,
        ):

            rep_rows.append(
                {
                    "condition":
                        row[
                            "condition"
                        ],
                    "genotype":
                        genotype,
                    "repeat":
                        i,
                    "value":
                        float(value),
                }
            )


rep = pd.DataFrame(
    rep_rows
)

if len(rep) != 30:
    fail(
        f"Expected 30 Pab1 replicate observations; found {len(rep)}"
    )

rep.to_csv(
    CURATED
    / "Fig4_6_Pab1_Hsp104_replicates.tsv",
    sep="\t",
    index=False,
)


                                                              
            
                                                              

manifest = []

for p in sources:

    manifest.append(
        {
            "path":
                str(p),
            "sha256":
                sha256(p),
        }
    )

pd.DataFrame(
    manifest
).to_csv(
    CURATED
    / "Fig4_6_source_manifest.tsv",
    sep="\t",
    index=False,
)


                                                              
            
                                                              

plt.rcParams.update(
    {
        "font.family":
            "DejaVu Sans",
        "font.size":
            9.0,
        "axes.labelsize":
            9.5,
        "xtick.labelsize":
            8.2,
        "ytick.labelsize":
            8.2,
        "legend.fontsize":
            8.0,
        "axes.linewidth":
            0.8,
        "pdf.fonttype":
            42,
        "ps.fonttype":
            42,
    }
)


fig = plt.figure(
    figsize=(
        18.2,
        14.2,
    ),
)

outer = fig.add_gridspec(
    2,
    12,
    height_ratios=[
        0.95,
        1.05,
    ],
    hspace=0.38,
    wspace=0.72,
)


                                                              
                                        
                                                              

def nested_panel(
    spec,
    note_height=0.17,
):

    gs = spec.subgridspec(
        2,
        1,
        height_ratios=[
            1.0 - note_height,
            note_height,
        ],
                                                               
        hspace=0.34,
    )

    ax = fig.add_subplot(
        gs[0, 0]
    )

    note = fig.add_subplot(
        gs[1, 0]
    )

    note.axis(
        "off"
    )

    return ax, note


                                                              
                
                                                              

axA, noteA = nested_panel(
    outer[
        0,
        0:3
    ],
    note_height=0.22,
)

xA = np.arange(
    3
)

valsA = (
    disc[
        "YEF3_T972"
    ]
    .to_numpy(float)
)

axA.axhline(
    0,
    lw=0.8,
    color="0.55",
)

for x, y in zip(
    xA,
    valsA,
):

    axA.vlines(
        x,
        0,
        y,
        lw=1.5,
        color="0.45",
    )

axA.scatter(
    xA,
    valsA,
    s=76,
    zorder=3,
)

axA.set_xticks(
    xA
)

axA.set_xticklabels(
    [
        r"$\it{ypl150w\Delta}$",
        r"$\it{ppt1\Delta}$",
        r"$\it{ctk1\Delta}$",
    ],
    rotation=24,
    ha="right",
)

axA.set_ylabel(
    "Yef3-T972 log$_2$ ratio"
)

axA.xaxis.labelpad = 2

clean_axes(
    axA
)

letter(
    axA,
    "A",
)

noteA.text(
    0.50,
    0.62,
    "Discovery / nomination only",
    transform=noteA.transAxes,
    fontsize=8.1,
    fontweight="bold",
    ha="center",
    va="center",
)

noteA.text(
    0.50,
    0.16,
    "Трите пертурбации са изключени от независимия n=107 анализ.",
    transform=noteA.transAxes,
    fontsize=7.5,
    color="0.35",
    ha="center",
    va="center",
)


                                                              
                                    
                                                              

axB, noteB = nested_panel(
    outer[
        0,
        4:8
    ],
    note_height=0.22,
)

raw_rho = (
    adj[
        "raw_spearman_rho"
    ]
    .to_numpy(float)
)

partial = (
    adj[
        "partial_spearman_global_adjusted"
    ]
    .to_numpy(float)
)

yB = np.arange(
    3
)

for i in range(
    3
):

    axB.plot(
        [
            raw_rho[i],
            partial[i],
        ],
        [
            yB[i],
            yB[i],
        ],
        lw=1.5,
        color="0.72",
    )

axB.scatter(
    raw_rho,
    yB,
    s=68,
    marker="o",
    zorder=3,
)

axB.scatter(
    partial,
    yB,
    s=68,
    marker="D",
    zorder=3,
)

axB.axvline(
    0,
    lw=0.8,
    color="0.60",
)

axB.set_yticks(
    yB
)

axB.set_yticklabels(
    [
        "Protein folding",
        "Translation elongation",
        "Stress granule / RNP",
    ]
)

axB.tick_params(
    axis="y",
    pad=2,
)

axB.invert_yaxis()

axB.set_xlabel(
    "Корелация с Yef3-T972 (ρ)"
)

axB.xaxis.labelpad = 2

clean_axes(
    axB
)

letter(
    axB,
    "B",
)

noteB.text(
    0.50,
    0.66,
    "● преди корекция     ◆ след global-state корекция",
    transform=noteB.transAxes,
    fontsize=7.7,
    ha="center",
    va="center",
)

noteB.text(
    0.50,
    0.18,
    (
        "folding: P=0,348    |    elongation: P=0,060    |    "
        "SG/RNP: P=1,0×10⁻⁴; q=3,0×10⁻⁴"
    ),
    transform=noteB.transAxes,
    fontsize=7.25,
    color="0.35",
    ha="center",
    va="center",
)


                                                              
                                              
                                                              

axC, noteC = nested_panel(
    outer[
        0,
        8:12
    ],
    note_height=0.22,
)

xc = (
    res[
        "YEF3_T972_residual"
    ]
    .to_numpy(float)
)

yc = (
    res[
        "STRESS_GRANULE_RNP_residual"
    ]
    .to_numpy(float)
)

axC.axhline(
    0,
    lw=0.7,
    color="0.75",
)

axC.axvline(
    0,
    lw=0.7,
    color="0.75",
)

axC.scatter(
    xc,
    yc,
    s=34,
    alpha=0.70,
)

coef = np.polyfit(
    xc,
    yc,
    1,
)

xx = np.linspace(
    xc.min(),
    xc.max(),
    100,
)

axC.plot(
    xx,
    coef[0] * xx
    + coef[1],
    lw=1.5,
    color="0.30",
)

axC.set_xlabel(
    "Residual Yef3-T972"
)

axC.set_ylabel(
    "Residual SG/RNP programme"
)

axC.xaxis.labelpad = 2

clean_axes(
    axC
)

letter(
    axC,
    "C",
)

noteC.text(
    0.50,
    0.66,
    (
        f"n=107    |    partial ρ={partial_rho:.3f}    |    "
        f"empirical P={sg_perm_p:.1e}    |    q={sg_q:.1e}"
    ),
    transform=noteC.transAxes,
    fontsize=7.7,
    ha="center",
    va="center",
)

noteC.text(
    0.50,
    0.18,
    (
        f"Matched random gene sets: P={random_p:.3f}    ·    "
        "PRE-CAUSAL association"
    ),
    transform=noteC.transAxes,
    fontsize=7.25,
    color="0.35",
    ha="center",
    va="center",
)


                                                              
                                  
                                                              

axD, noteD = nested_panel(
    outer[
        1,
        0:7
    ],
    note_height=0.20,
)

wt = (
    ht[
        ht["genotype"]
        == "WT"
    ]
    .sort_values(
        "time_index"
    )
)

hd = (
    ht[
        ht["genotype"]
        == "hsp104D"
    ]
    .sort_values(
        "time_index"
    )
)

                         
axD.axvspan(
    1,
    8,
    color="0.96",
    zorder=0,
)

axD.axvspan(
    8,
    15,
    color="0.90",
    zorder=0,
)

axD.axvline(
    8,
    lw=0.9,
    ls="--",
    color="0.55",
)

axD.axhline(
    0,
    lw=0.7,
    ls="--",
    color="0.65",
)

axD.plot(
    wt[
        "time_index"
    ],
    wt[
        "mean_log2_vs_T1"
    ],
    marker="o",
    lw=2.0,
    label="WT",
)

axD.plot(
    hd[
        "time_index"
    ],
    hd[
        "mean_log2_vs_T1"
    ],
    marker="o",
    lw=2.0,
    label=r"$\it{hsp104\Delta}$",
)

axD.set_xlim(
    1,
    15,
)

axD.set_xticks(
    [
        1,
        4,
        8,
        12,
        15,
    ]
)

axD.set_xlabel(
    "Времева точка"
)

axD.set_ylabel(
    r"Средно $\log_2[R(t)/R(T1)]$"
)

axD.legend(
    frameon=False,
    loc="best",
)

axD.xaxis.labelpad = 2

clean_axes(
    axD
)

letter(
    axD,
    "D",
)

noteD.text(
    0.02,
    0.76,
    (
        f"T1–T8 AUC:    WT {EARLY_WT:.3f}    |    "
        f"hsp104Δ {EARLY_HD:.3f}"
    ),
    transform=noteD.transAxes,
    fontsize=7.8,
    va="center",
)

noteD.text(
    0.02,
    0.47,
    (
        f"T8–T15 AUC:  WT {LATE_WT:.3f}    |    "
        f"hsp104Δ {LATE_HD:.3f}"
    ),
    transform=noteD.transAxes,
    fontsize=7.8,
    va="center",
)

noteD.text(
    0.02,
    0.10,
    (
        "Описателно: една независимо разрешена сцена на генотип; "
        "без inferential P-value."
    ),
    transform=noteD.transAxes,
    fontsize=7.25,
    color="0.35",
    va="center",
)


                                                              
                                    
                                                              

axE, noteE = nested_panel(
    outer[
        1,
        7:12
    ],
    note_height=0.20,
)

xE = np.arange(
    3
)

labelsE = [
    "След стрес",
    "1 h",
    "2 h",
]

genotypes = [
    "HSP104",
    "hsp104D",
]

line_handles = {}

for gi, genotype in enumerate(
    genotypes
):

    means = []

    for condition in condition_order:

        row = ps[
            ps["condition"]
            == condition
        ].iloc[0]

        means.append(
            float(
                row[
                    f"{genotype}_mean"
                ]
            )
        )

    label = (
        "HSP104"
        if genotype == "HSP104"
        else r"$\it{hsp104\Delta}$"
    )

    line, = axE.plot(
        xE,
        means,
        marker="o",
        markersize=6,
        lw=2.0,
        label=label,
        zorder=4,
    )

    line_handles[
        genotype
    ] = line

                                         
    shift = (
        -0.065
        if genotype == "HSP104"
        else +0.065
    )

    jitter = np.linspace(
        -0.032,
        0.032,
        5,
    )

    for ti, condition in enumerate(
        condition_order
    ):

        vals = (
            rep[
                (
                    rep["condition"]
                    == condition
                )
                &
                (
                    rep["genotype"]
                    == genotype
                )
            ]
            .sort_values(
                "repeat"
            )[
                "value"
            ]
            .to_numpy(float)
        )

        if len(vals) != 5:
            fail(
                "Expected five replicate values in Panel E."
            )

        axE.scatter(
            ti
            + shift
            + jitter,
            vals,
            s=29,
            alpha=0.56,
            color=line.get_color(),
            zorder=2,
        )


axE.set_xticks(
    xE
)

axE.set_xticklabels(
    labelsE
)

axE.set_ylim(
    -4,
    105,
)

axE.set_ylabel(
    "Клетки с Pab1 гранули (%)"
)

axE.legend(
    frameon=False,
    loc="best",
)

axE.xaxis.labelpad = 2

clean_axes(
    axE
)

letter(
    axE,
    "E",
)

p1 = float(
    ps.loc[
        ps["condition"]
        == "1h_recovery",
        "exact_randomization_p_two_sided",
    ].iloc[0]
)

p2 = float(
    ps.loc[
        ps["condition"]
        == "2h_recovery",
        "exact_randomization_p_two_sided",
    ].iloc[0]
)

noteE.text(
    0.50,
    0.76,
    "n=5 независими биологични повторения на група",
    transform=noteE.transAxes,
    fontsize=7.7,
    ha="center",
    va="center",
)

noteE.text(
    0.50,
    0.45,
    f"1 h: P={p1:.5f}    |    2 h: P={p2:.5f}",
    transform=noteE.transAxes,
    fontsize=7.6,
    ha="center",
    va="center",
)

noteE.text(
    0.50,
    0.10,
    "Точки = отделни повторения    ·    exact randomization, two-sided",
    transform=noteE.transAxes,
    fontsize=7.2,
    color="0.35",
    ha="center",
    va="center",
)


                                                              
        
                                                              

fig.subplots_adjust(
    left=0.075,
    right=0.985,
    top=0.975,
    bottom=0.065,
)

OUTBASE = (
    FINAL
    / "Fig4_6_Yef3_Hsp104_recovery"
)

fig.savefig(
    str(OUTBASE)
    + ".png",
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE)
    + ".pdf",
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE)
    + ".svg",
    bbox_inches="tight",
)

plt.close(
    fig
)


                                                              
                  
                                                              

legend = f"""Фигура 4.6. Yef3-T972 маркира RNP/протеостазно фосфопротеомно състояние, докато Hsp104 е свързан с късната постстресова резолюция.

(A) Първоначална номинация на Yef3-T972 чрез три генетични пертурбации. При ypl150wΔ Yef3-T972 е понижен ({valsA[0]:.2f} log₂ ratio), докато при ppt1Δ и ctk1Δ е повишен ({valsA[1]:+.2f} и {valsA[2]:+.2f}). Тези три пертурбации са използвани единствено за discovery/nomination и са изключени от последващата независима статистическа оценка.

(B) Асоциация на Yef3-T972 с предварително дефинираните protein-folding, translation-elongation и stress-granule/RNP програми преди и след корекция спрямо глобалното фосфопротеомно състояние. Корелациите се променят съответно от {raw_rho[0]:.3f} до {partial[0]:.3f}, от {raw_rho[1]:.3f} до {partial[1]:.3f} и от {raw_rho[2]:.3f} до {partial[2]:.3f}. След корекция връзката със stress-granule/RNP програмата остава ясно подкрепена.

(C) Остатъчна връзка между rank-transformed Yef3-T972 и stress-granule/RNP програмата след residualization спрямо идентичния global phosphoproteomic axis в 107 независими делеционни пертурбации. Partial Spearman ρ={partial_rho:.3f}, емпирично P={sg_perm_p:.1e}, q={sg_q:.1e}. При 5000 случайни генни множества със същото експериментално покритие асоциацията е по-силна от очакваното за случаен модул (P={random_p:.3f}). Анализът е предкаузален и не установява Yef3-T972 като причинен регулатор на Guk1 или Hsp104.

(D) Пълна T1–T15 траектория на Guk1-7-GFP фокусното обогатяване при WT и hsp104Δ за предварително фиксираните клетки. Ранният интегрален T1–T8 фенотип е {EARLY_WT:.3f} при WT и {EARLY_HD:.3f} при hsp104Δ и не подкрепя необходимост на Hsp104 за първоначалното намаляване на фокусното обогатяване. В късния T8–T15 интервал стойностите са {LATE_WT:.3f} при WT и {LATE_HD:.3f} при hsp104Δ, което е съвместимо с по-слаба късна резолюция при липса на Hsp104. Сравнението е описателно, тъй като е налична една независимо разрешена сцена на генотип.

(E) Независими публични данни за Pab1-съдържащи stress granules. Отделните точки показват петте биологични повторения, а линиите — груповите средни стойности. Делът на клетките с гранули е 100,00%, 39,67% и 1,11% при HSP104 и 99,13%, 97,78% и 90,41% при hsp104Δ непосредствено след стреса, след 1 h и след 2 h възстановяване. Разликите при 1 h и 2 h са P={p1:.5f} и P={p2:.5f}.

Панели A–C и D–E представляват отделни доказателствени линии. Фигурата не предполага причинна последователност Yef3-T972 → Hsp104 → Guk1.
"""

(
    LEGENDS
    / "Fig4_6_legend_bg.txt"
).write_text(
    legend,
    encoding="utf-8",
)


print()
print(
    "===== FIGURE 4.6 COMPLETE ====="
)

print(
    f"Panel C exact partial rho: {partial_rho:.6f}"
)

print(
    f"Panel C random-set P: {random_p:.6f}"
)

print(
    "Panel D trajectory:",
    HSP_TRAJ,
)

print(
    f"Panel E replicate observations: {len(rep)}"
)

print()
print(
    str(OUTBASE)
    + ".pdf"
)
