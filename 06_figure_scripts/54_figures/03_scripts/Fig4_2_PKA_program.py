#!/usr/bin/env python3

from pathlib import Path
import hashlib
import math

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


                                                              
       
                                                              

ROOT = Path("/")

D32 = ROOT / "32_pka_direct_site_trajectory_corrected"
D17 = ROOT / "17_networkin_temporal_discovery"
D37 = ROOT / "37_final_statistical_freeze"

DIRECT_SUMMARY = D32 / "CORRECTED_DIRECT_A_TRAJECTORY_SUMMARY.tsv"
DIRECT_GROUPS = D32 / "DIRECT_A_10_GROUP_AUCs.tsv"
DIRECT_GENES = D32 / "DIRECT_A_8_GENE_AUCs.tsv"
DIRECT_LOO = D32 / "DIRECT_A_leave_one_gene_out.tsv"

NETWORKIN = D17 / "networkin_global_heat_cold_tests.tsv"
NETWORKIN_SUMMARY = D17 / "NETWORKIN_TEMPORAL_SUMMARY.tsv"

CLAIMS = D37 / "FINAL_CLAIM_MATRIX.tsv"

DIRECT_EMPIRICAL_NULL = (
    ROOT
    / "54_figures"
    / "02_curated_inputs"
    / "Fig4_2"
    / "DIRECT_A_corrected_empirical_null.tsv.gz"
)

OUTROOT = ROOT / "54_figures"

CURATED = OUTROOT / "02_curated_inputs" / "Fig4_2"
FINAL = OUTROOT / "05_final"
LEGENDS = OUTROOT / "06_legends"

for d in [CURATED, FINAL, LEGENDS]:
    d.mkdir(parents=True, exist_ok=True)


                                                              
                  
                                                              

DIRECT_COLOR = "#C84B31"
OPPOSITE_COLOR = "#3478A8"
ORTHO_COLOR = "#765AA5"

NULL_FILL = "#D0D0D0"
NULL_LINE = "#8A8A8A"
OTHER_COLOR = "#AFAFAF"
ZERO_COLOR = "#707070"

SIG_EDGE = "#222222"


                                                              
           
                                                              

def fail(msg):
    raise RuntimeError(msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def normal_pdf(x, mu, sd):
    return (
        np.exp(-0.5 * ((x - mu) / sd) ** 2)
        / (sd * np.sqrt(2.0 * np.pi))
    )


def sci(x, digits=1):
    if x == 0:
        return "0"

    exponent = int(math.floor(math.log10(abs(x))))
    mantissa = x / (10 ** exponent)

    return (
        rf"${mantissa:.{digits}f}"
        rf"\times10^{{{exponent}}}$"
    )


def kinase_label(x):
    x = str(x)

    if x == "TPK2_TPK1_TPK3_group":
        return "TPK1/TPK2/TPK3"

    x = x.replace("_group", "")
    x = x.replace("_", "/")

    return x


def italic_gene(g):
    return rf"$\it{{{g}}}$"


                                                              
           
                                                              

for p in [
    DIRECT_SUMMARY,
    DIRECT_GROUPS,
    DIRECT_GENES,
    DIRECT_LOO,
    NETWORKIN,
    NETWORKIN_SUMMARY,
    CLAIMS,
    DIRECT_EMPIRICAL_NULL,
]:
    if not p.exists():
        fail(f"Missing authoritative input: {p}")


                                                              
      
                                                              

summary = pd.read_csv(DIRECT_SUMMARY, sep="\t")
groups = pd.read_csv(DIRECT_GROUPS, sep="\t")
genes = pd.read_csv(DIRECT_GENES, sep="\t")
loo = pd.read_csv(DIRECT_LOO, sep="\t")

net = pd.read_csv(NETWORKIN, sep="\t")
net_summary = pd.read_csv(NETWORKIN_SUMMARY, sep="\t")

claims = pd.read_csv(CLAIMS, sep="\t")

direct_null_df = pd.read_csv(
    DIRECT_EMPIRICAL_NULL,
    sep="\t",
)

if "null_mean_gene_AUC" not in direct_null_df.columns:
    fail(
        "Exact corrected empirical null file lacks "
        "null_mean_gene_AUC."
    )

direct_null_values = pd.to_numeric(
    direct_null_df["null_mean_gene_AUC"],
    errors="raise",
).to_numpy()

if len(direct_null_values) != 100000:
    fail(
        f"Expected 100000 corrected null permutations; "
        f"found {len(direct_null_values)}"
    )


                                                              
                       
                                                              

if len(summary) != 1:
    fail("Corrected direct-A summary must contain exactly one row.")

S = summary.iloc[0]

required_summary = {
    "n_direct_A_groups",
    "n_direct_A_genes",
    "observed_mean_gene_AUC",
    "observed_median_gene_AUC",
    "negative_genes",
    "null_mean",
    "null_sd",
    "empirical_z",
    "mean_empirical_lower_p",
    "mean_empirical_two_sided_p",
    "all_leave_one_out_means_negative",
    "all_leave_one_out_p_lt_0_05",
    "worst_leave_one_out_p",
}

missing = required_summary - set(summary.columns)

if missing:
    fail(
        "Corrected summary missing columns: "
        + ", ".join(sorted(missing))
    )


n_groups = int(S["n_direct_A_groups"])
n_genes = int(S["n_direct_A_genes"])

obs_direct = float(S["observed_mean_gene_AUC"])
null_direct = float(S["null_mean"])
sd_direct = float(S["null_sd"])
z_direct = float(S["empirical_z"])

p_direct_one = float(S["mean_empirical_lower_p"])
p_direct_two = float(S["mean_empirical_two_sided_p"])

negative_genes = int(S["negative_genes"])


if n_groups != 10:
    fail(f"Expected 10 direct-A groups; found {n_groups}")

if n_genes != 8:
    fail(f"Expected 8 direct-A genes; found {n_genes}")

if len(groups) != 10:
    fail(f"DIRECT_A_10_GROUP_AUCs has {len(groups)} rows")

if genes["gene"].nunique() != 8:
    fail(
        "DIRECT_A_8_GENE_AUCs does not contain exactly 8 genes."
    )

if negative_genes != 7:
    fail(
        f"Expected 7/8 negative genes in corrected analysis; "
        f"found {negative_genes}"
    )

if abs(obs_direct - (-38.5899)) > 0.02:
    fail(
        f"Unexpected corrected direct mean: {obs_direct}"
    )

if abs(null_direct - (-2.7201)) > 0.02:
    fail(
        f"Unexpected corrected null mean: {null_direct}"
    )

if abs(z_direct - (-4.0378)) > 0.02:
    fail(
        f"Unexpected corrected Z: {z_direct}"
    )


                                            
cdc = genes[
    genes["gene"].astype(str).str.upper() == "CDC19"
]

if len(cdc) != 1:
    fail("CDC19 is not uniquely represented.")

if float(cdc.iloc[0]["direct_A_gene_AUC"]) <= 0:
    fail(
        "Corrected CDC19 gene-level AUC is not positive."
    )

negative_actual = (
    pd.to_numeric(
        genes["direct_A_gene_AUC"],
        errors="raise",
    )
    < 0
).sum()

if negative_actual != 7:
    fail(
        f"DIRECT_A_8_GENE_AUCs contains "
        f"{negative_actual} negative genes, expected 7."
    )


                                                              
                               
                                                              

if loo["omitted_gene"].nunique() != 8:
    fail(
        "Corrected leave-one-gene-out table "
        "does not contain exactly 8 genes."
    )

loo["remaining_mean_AUC"] = pd.to_numeric(
    loo["remaining_mean_AUC"],
    errors="raise",
)

loo["empirical_lower_p"] = pd.to_numeric(
    loo["empirical_lower_p"],
    errors="raise",
)

if not (loo["remaining_mean_AUC"] < 0).all():
    fail("Not all corrected LOO module means are negative.")

if not (loo["empirical_lower_p"] < 0.05).all():
    fail("Not all corrected LOO tests have P < 0.05.")


                                                              
                          
                                                              

required_net = {
    "analysis_label",
    "predicted_kinase_group",
    "analysis_type",
    "n_target_sites",
    "n_target_genes",
    "mean_heat_minus_cold",
    "empirical_z",
    "empirical_p",
    "null_mean",
    "null_sd",
    "q_across_standard_groups",
}

missing = required_net - set(net.columns)

if missing:
    fail(
        "NetworKIN table missing columns: "
        + ", ".join(sorted(missing))
    )


orth = net[
    net["analysis_type"]
    .astype(str)
    .eq("orthogonal_positive_control")
].copy()

if len(orth) != 1:
    fail(
        f"Expected exactly one orthogonal positive-control row; "
        f"found {len(orth)}"
    )

O = orth.iloc[0]

orth_groups = int(O["n_target_sites"])
orth_genes = int(O["n_target_genes"])

orth_effect = float(O["mean_heat_minus_cold"])
orth_null = float(O["null_mean"])
orth_sd = float(O["null_sd"])
orth_z = float(O["empirical_z"])
orth_p = float(O["empirical_p"])

if orth_groups != 20:
    fail(f"Expected 20 orthogonal groups; found {orth_groups}")

if orth_genes != 16:
    fail(f"Expected 16 orthogonal genes; found {orth_genes}")

if abs(orth_effect - (-1.279195)) > 0.002:
    fail(f"Unexpected orthogonal PKA effect: {orth_effect}")

if abs(orth_p - 1.4e-4) > 1e-6:
    fail(f"Unexpected orthogonal PKA P: {orth_p}")


                                                              
                              
                                                              

standard = net[
    net["analysis_type"]
    .astype(str)
    .eq("standard")
].copy()

standard["mean_heat_minus_cold"] = pd.to_numeric(
    standard["mean_heat_minus_cold"],
    errors="raise",
)

standard["q_across_standard_groups"] = pd.to_numeric(
    standard["q_across_standard_groups"],
    errors="coerce",
)

standard = standard[
    standard["status"].astype(str).eq("TESTED")
].copy()

if len(standard) != 13:
    fail(
        f"Expected 13 standard NetworKIN groups; "
        f"found {len(standard)}"
    )

n_sig = int(
    (
        standard["q_across_standard_groups"]
        < 0.05
    ).sum()
)

if n_sig != 1:
    fail(
        f"Expected exactly one q<0.05 standard kinase module; "
        f"found {n_sig}"
    )

pka_standard = standard[
    standard["predicted_kinase_group"]
    .astype(str)
    .eq("TPK2_TPK1_TPK3_group")
]

if len(pka_standard) != 1:
    fail(
        "Could not uniquely identify standard NetworKIN PKA family."
    )

if float(
    pka_standard.iloc[0][
        "q_across_standard_groups"
    ]
) >= 0.05:
    fail(
        "Standard NetworKIN PKA family is not FDR significant."
    )


                                                              
            
                                                              

manifest = []

for role, path in [
    ("corrected_direct_summary", DIRECT_SUMMARY),
    ("corrected_direct_10_groups", DIRECT_GROUPS),
    ("corrected_direct_8_genes", DIRECT_GENES),
    ("corrected_direct_leave_one_out", DIRECT_LOO),
    ("networkin_global_heat_cold_tests", NETWORKIN),
    ("networkin_temporal_summary", NETWORKIN_SUMMARY),
    ("final_claim_matrix", CLAIMS),
    ("corrected_direct_empirical_null", DIRECT_EMPIRICAL_NULL),
]:
    manifest.append(
        {
            "role": role,
            "path": str(path),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        }
    )

pd.DataFrame(manifest).to_csv(
    CURATED / "Fig4_2_source_manifest.tsv",
    sep="\t",
    index=False,
)

genes.to_csv(
    CURATED / "Fig4_2_direct_gene_effects.tsv",
    sep="\t",
    index=False,
)

loo.to_csv(
    CURATED / "Fig4_2_leave_one_gene_out.tsv",
    sep="\t",
    index=False,
)

standard.to_csv(
    CURATED / "Fig4_2_networkin_standard_modules.tsv",
    sep="\t",
    index=False,
)


                                                              
                 
                                                              

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.labelsize": 9.4,
        "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3,
        "legend.fontsize": 8.0,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


fig = plt.figure(
    figsize=(18.0, 12.8),
    layout="constrained",
)

outer = fig.add_gridspec(
    2,
    1,
    height_ratios=[1.0, 1.18],
)

top = outer[0].subgridspec(
    1,
    3,
    width_ratios=[1.15, 1.05, 1.10],
    wspace=0.18,
)

bottom = outer[1].subgridspec(
    1,
    2,
    width_ratios=[1.0, 1.45],
    wspace=0.20,
)


                                                              
                                
                                                              

cellA = top[0, 0].subgridspec(
    2, 1,
    height_ratios=[1.0, 0.14],
    hspace=0.04,
)

axA = fig.add_subplot(cellA[0, 0])
noteA = fig.add_subplot(cellA[1, 0])
noteA.axis("off")

A = genes.copy()

A["direct_A_gene_AUC"] = pd.to_numeric(
    A["direct_A_gene_AUC"],
    errors="raise",
)

A = A.sort_values(
    "direct_A_gene_AUC",
    ascending=True,
).reset_index(drop=True)

y = np.arange(len(A))

for i, row in A.iterrows():

    gene = str(row["gene"])
    value = float(row["direct_A_gene_AUC"])

    color = (
        OPPOSITE_COLOR
        if gene.upper() == "CDC19"
        else DIRECT_COLOR
    )

    axA.hlines(
        y=i,
        xmin=min(0, value),
        xmax=max(0, value),
        color=color,
        lw=2.1,
        alpha=0.80,
    )

    axA.scatter(
        value,
        i,
        s=55,
        color=color,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )


labels = []

for g in A["gene"].astype(str):

    if g.upper() == "CDC19":
        labels.append(
            italic_gene(g) + "  (S22)"
        )
    else:
        labels.append(
            italic_gene(g)
        )

axA.set_yticks(y)
axA.set_yticklabels(labels)

axA.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.9,
)

axA.set_xlabel(
    r"$\Delta$AUC = AUC$_{Heat}$ − AUC$_{Cold}$"
)

noteA.text(
    0.02,
    0.55,
    "7/8 гена с отрицателен ефект; Cdc19-S22 е с противоположна посока",
    transform=noteA.transAxes,
    ha="left",
    va="center",
    fontsize=8.2,
)

axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

axA.text(
    -0.15,
    1.04,
    "A",
    transform=axA.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                               
                                                              

cellB = top[0, 1].subgridspec(
    2, 1,
    height_ratios=[1.0, 0.15],
    hspace=0.04,
)

axB = fig.add_subplot(cellB[0, 0])
noteB = fig.add_subplot(cellB[1, 0])
noteB.axis("off")


                                                         
q25, q75 = np.quantile(
    direct_null_values,
    [0.25, 0.75],
)

iqr = q75 - q25

if iqr > 0:
    bin_width = (
        2.0
        * iqr
        / (len(direct_null_values) ** (1.0 / 3.0))
    )

    n_bins = int(
        np.ceil(
            (
                direct_null_values.max()
                - direct_null_values.min()
            )
            / bin_width
        )
    )
else:
    n_bins = 60

n_bins = max(
    40,
    min(n_bins, 100),
)

axB.hist(
    direct_null_values,
    bins=n_bins,
    density=True,
    color=NULL_FILL,
    edgecolor="none",
    alpha=0.90,
)

axB.axvline(
    null_direct,
    color="#555555",
    ls="--",
    lw=1.2,
)

axB.axvline(
    obs_direct,
    color=DIRECT_COLOR,
    lw=2.4,
)

axB.set_xlabel(
    r"Средна генна $\Delta$AUC"
)

axB.set_ylabel(
    "Плътност"
)

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

axB.text(
    -0.15,
    1.04,
    "B",
    transform=axB.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)

noteB.text(
    0.02,
    0.68,
    (
        f"Observed ΔAUC = {obs_direct:.2f}   |   "
        f"8 гена / 10 групи   |   "
        f"Z = {z_direct:.2f}"
    ),
    transform=noteB.transAxes,
    ha="left",
    va="center",
    fontsize=8.1,
)

noteB.text(
    0.02,
    0.20,
    (
        f"двустранно P = {sci(p_direct_two, 1)}   |   "
        f"едностранно P = {sci(p_direct_one, 1)}   |   "
        f"100 000 matched permutations"
    ),
    transform=noteB.transAxes,
    ha="left",
    va="center",
    fontsize=7.8,
    color="#555555",
)


                                                              
                                   
                                                              

cellC = top[0, 2].subgridspec(
    2, 1,
    height_ratios=[1.0, 0.14],
    hspace=0.04,
)

axC = fig.add_subplot(cellC[0, 0])
noteC = fig.add_subplot(cellC[1, 0])
noteC.axis("off")

C = loo.copy()

                                                          
order_map = {
    g: i
    for i, g in enumerate(
        A["gene"].astype(str).tolist()
    )
}

C["_order"] = (
    C["omitted_gene"]
    .astype(str)
    .map(order_map)
)

C = (
    C.sort_values(
        ["_order", "omitted_gene"]
    )
    .reset_index(drop=True)
)

yc = np.arange(len(C))

axC.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.9,
)

axC.axvline(
    obs_direct,
    color=DIRECT_COLOR,
    lw=1.2,
    ls="--",
    alpha=0.75,
)

axC.scatter(
    C["remaining_mean_AUC"],
    yc,
    s=58,
    color="#4A4A4A",
    edgecolor="white",
    linewidth=0.5,
    zorder=3,
)

for i, row in C.iterrows():

    axC.hlines(
        i,
        min(
            obs_direct,
            row["remaining_mean_AUC"],
        ),
        max(
            obs_direct,
            row["remaining_mean_AUC"],
        ),
        color="#B5B5B5",
        lw=1.1,
        zorder=1,
    )


axC.set_yticks(yc)

axC.set_yticklabels(
    [
        italic_gene(g)
        for g in C["omitted_gene"]
        .astype(str)
    ]
)

axC.set_xlabel(
    r"Средна $\Delta$AUC след изключване"
)

axC.set_ylabel(
    "Изключен ген"
)

noteC.text(
    0.02,
    0.55,
    "8/8 leave-one-gene-out модула остават отрицателни; 8/8 P < 0,05",
    transform=noteC.transAxes,
    ha="left",
    va="center",
    fontsize=8.2,
)

axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

axC.text(
    -0.15,
    1.04,
    "C",
    transform=axC.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                                        
                                                              

cellD = bottom[0, 0].subgridspec(
    2, 1,
    height_ratios=[1.0, 0.14],
    hspace=0.04,
)

axD = fig.add_subplot(cellD[0, 0])
noteD = fig.add_subplot(cellD[1, 0])
noteD.axis("off")


                                                       
                                     
axD.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.8,
    alpha=0.65,
)

axD.errorbar(
    x=null_direct if False else orth_null,
    y=0,
    xerr=orth_sd,
    fmt="o",
    markersize=7,
    color=NULL_LINE,
    ecolor=NULL_LINE,
    elinewidth=2.0,
    capsize=5,
    label="Matched null: mean ± SD",
)

axD.scatter(
    orth_effect,
    1,
    s=85,
    marker="D",
    color=ORTHO_COLOR,
    edgecolor="white",
    linewidth=0.6,
    zorder=3,
    label="Observed PKA module",
)

axD.hlines(
    1,
    orth_null,
    orth_effect,
    color=ORTHO_COLOR,
    lw=1.1,
    alpha=0.55,
)

pad = max(
    0.20,
    0.25 * abs(orth_effect - orth_null),
)

xmin = min(
    orth_effect,
    orth_null - orth_sd,
) - pad

xmax = max(
    orth_effect,
    orth_null + orth_sd,
) + pad

axD.set_xlim(
    xmin,
    xmax,
)

axD.set_ylim(
    -0.55,
    1.55,
)

axD.set_yticks(
    [0, 1]
)

axD.set_yticklabels(
    [
        "matched null",
        "observed",
    ]
)

axD.set_xlabel(
    "Среден Heat–Cold ефект"
)

axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)
axD.spines["left"].set_visible(False)

axD.tick_params(
    axis="y",
    length=0,
)

axD.text(
    -0.12,
    1.04,
    "D",
    transform=axD.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)

noteD.text(
    0.02,
    0.62,
    (
        f"Observed = {orth_effect:.2f}   |   "
        f"16 гена / 20 групи   |   "
        f"Z = {orth_z:.2f}   |   "
        f"P = {sci(orth_p, 1)}"
    ),
    transform=noteD.transAxes,
    ha="left",
    va="center",
    fontsize=8.2,
)

noteD.text(
    0.02,
    0.16,
    (
        f"matched-null mean = {orth_null:.2f}; "
        f"SD = {orth_sd:.2f}"
    ),
    transform=noteD.transAxes,
    ha="left",
    va="center",
    fontsize=7.8,
    color="#555555",
)


                                                              
                                   
                                                              

cellE = bottom[0, 1].subgridspec(
    2, 1,
    height_ratios=[1.0, 0.14],
    hspace=0.04,
)

axE = fig.add_subplot(cellE[0, 0])
noteE = fig.add_subplot(cellE[1, 0])
noteE.axis("off")

E = standard.copy()

E["_label"] = (
    E["predicted_kinase_group"]
    .map(kinase_label)
)

E = (
    E.sort_values(
        "mean_heat_minus_cold",
        ascending=True,
    )
    .reset_index(drop=True)
)

ye = np.arange(len(E))

axE.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.9,
)

for i, row in E.iterrows():

    is_pka = (
        row["predicted_kinase_group"]
        == "TPK2_TPK1_TPK3_group"
    )

    is_sig = (
        row["q_across_standard_groups"]
        < 0.05
    )

    color = (
        ORTHO_COLOR
        if is_pka
        else OTHER_COLOR
    )

    marker = (
        "D"
        if is_pka
        else "o"
    )

    size = (
        78
        if is_pka
        else 46
    )

    axE.hlines(
        i,
        0,
        row["mean_heat_minus_cold"],
        color=color,
        lw=1.7 if is_pka else 1.0,
        alpha=0.82,
    )

    axE.scatter(
        row["mean_heat_minus_cold"],
        i,
        s=size,
        marker=marker,
        color=color,
        edgecolor=(
            SIG_EDGE
            if is_sig
            else "white"
        ),
        linewidth=(
            1.1
            if is_sig
            else 0.45
        ),
        zorder=3,
    )


axE.set_yticks(ye)
axE.set_yticklabels(
    E["_label"]
)

axE.set_xlabel(
    "Среден Heat–Cold ефект"
)

pka_e = E[
    E["predicted_kinase_group"]
    == "TPK2_TPK1_TPK3_group"
].iloc[0]

pka_i = int(
    E.index[
        E["predicted_kinase_group"]
        == "TPK2_TPK1_TPK3_group"
    ][0]
)

noteE.text(
    0.02,
    0.55,
    (
        f"PKA-family: q = "
        f"{float(pka_e['q_across_standard_groups']):.3g}"
        f"   |   FDR q < 0,05: {n_sig}/{len(E)} модула"
    ),
    transform=noteE.transAxes,
    ha="left",
    va="center",
    fontsize=8.2,
)

axE.spines["top"].set_visible(False)
axE.spines["right"].set_visible(False)

axE.text(
    -0.09,
    1.04,
    "E",
    transform=axE.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                      
                                                              


                                                              
        
                                                              

OUTBASE = (
    FINAL
    / "Fig4_2_PKA_temperature_program"
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


                                                              
               
                                                              

legend = f"""Фигура 4.2. PKA е най-силно поддържаната киназно-свързана температурна програма.

(A) Генни Heat–Cold ΔAUC стойности за осемте гена, представени сред експериментално потвърдените директни PKA субстрати. Седем от осемте гена показват отрицателен интегрален Heat–Cold ефект, докато Cdc19-S22 е с противоположна посока.

(B) Позиция на наблюдавания среден генен ΔAUC ({obs_direct:.2f}) спрямо съпоставимия matched-null модел (μ={null_direct:.2f}, σ={sd_direct:.2f}). Коригираният директен PKA модул включва 10 фосфопротеомни групи в 8 гена (Z={z_direct:.2f}; двустранно емпирично P={p_direct_two:.2g}; едностранно P={p_direct_one:.2g}). Сивата хистограма представя точното architecture-matched empirical null разпределение от 100 000 пермутации, възстановено детерминистично от коригирания анализ с оригиналния фиксиран random seed.

(C) Leave-one-gene-out анализ на директния PKA модул. Точките показват средния ΔAUC след последователно изключване на всеки от осемте гена; прекъснатата линия показва ефекта на пълния модул. Всички осем редуцирани модула запазват отрицателен ефект и емпирично P<0,05.

(D) Ортогонален NetworKIN PKA-family модул след отстраняване на всички фосфосайтове с предварително налична PKA A/B подкрепа. Модулът включва 20 фосфопротеомни групи в 16 гена и показва среден Heat–Cold ефект {orth_effect:.2f} (Z={orth_z:.2f}; емпирично P={orth_p:.2g}). Matched-null моделът е представен чрез неговите frozen mean ± SD, без да се предполага конкретна форма на null разпределението.

(E) Киназно-широк NetworKIN анализ на всички 13 стандартни модула, преминали предварително зададените критерии за размер и генно покритие. PKA-family модулът е единственият, запазващ статистическа значимост след корекция за множествено тестване (FDR q<0,05); останалите модули са показани в сиво.
"""

legend_path = (
    LEGENDS
    / "Fig4_2_legend_bg.txt"
)

legend_path.write_text(
    legend,
    encoding="utf-8",
)


                                                              
            
                                                              

print()
print("===== FIGURE 4.2 COMPLETE =====")
print()
if abs(
    direct_null_values.mean() - null_direct
) > 1e-6:
    fail(
        "Recovered empirical-null mean does not match "
        f"frozen summary: {direct_null_values.mean()} "
        f"vs {null_direct}"
    )

if abs(
    direct_null_values.std(ddof=1) - sd_direct
) > 1e-6:
    fail(
        "Recovered empirical-null SD does not match "
        f"frozen summary: {direct_null_values.std(ddof=1)} "
        f"vs {sd_direct}"
    )


print(
    f"Direct PKA: {n_groups} groups / {n_genes} genes; "
    f"mean={obs_direct:.4f}; "
    f"null={null_direct:.4f}; "
    f"Z={z_direct:.4f}; "
    f"P2={p_direct_two:.3g}"
)
print(
    f"Orthogonal PKA: {orth_groups} groups / {orth_genes} genes; "
    f"effect={orth_effect:.4f}; "
    f"P={orth_p:.3g}"
)
print(
    f"Standard NetworKIN modules: {len(standard)}; "
    f"q<0.05={n_sig}"
)
print()
print(str(OUTBASE) + ".png")
print(str(OUTBASE) + ".pdf")
print(str(OUTBASE) + ".svg")
print(legend_path)
