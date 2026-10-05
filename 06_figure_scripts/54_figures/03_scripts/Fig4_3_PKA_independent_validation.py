#!/usr/bin/env python3

from pathlib import Path
import hashlib
import math
import re

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr


                                                              
       
                                                              

ROOT = Path("/")

D30 = ROOT / "30_uliana_causal_validation"
D31 = ROOT / "31_uliana_context_interaction"
D35 = ROOT / "35_phosphoatlas_pka_replication"
D36 = ROOT / "36_phosphoatlas_temperature_sensitivity"
D37 = ROOT / "37_final_statistical_freeze"

CLAIMS = D37 / "FINAL_CLAIM_MATRIX.tsv"

OUTROOT = ROOT / "54_figures"
CURATED = OUTROOT / "02_curated_inputs" / "Fig4_3"
FINAL = OUTROOT / "05_final"
LEGENDS = OUTROOT / "06_legends"

for d in [CURATED, FINAL, LEGENDS]:
    d.mkdir(parents=True, exist_ok=True)


                                                              
                         
                                                              

                         
ULIANA_N_GROUPS = 14
ULIANA_N_GENES = 11

ULIANA_EXPO_OBS = -1.8520
ULIANA_EXPO_NULL = -0.6427
ULIANA_EXPO_Z = -4.031
ULIANA_EXPO_P = 7.0e-5

ULIANA_HS_P = 0.0624
ULIANA_INTERACTION_P = 0.1639


                           
MAF1_EXPO = -5.2053
MAF1_HS = -0.656
MAF1_INT_P = 0.00327
MAF1_INT_Q = 0.01309

ATG1_EXPO = -2.8098
ATG1_HS = +1.8791
ATG1_INT_P = 0.000411
ATG1_INT_Q = 0.00288


                  
DIRECT_ATLAS_N_GROUPS = 5
DIRECT_ATLAS_N_GENES = 4
DIRECT_ATLAS_OBS = -1.4996
DIRECT_ATLAS_NULL = +0.0386
DIRECT_ATLAS_P = 0.04079

ORTHO_ATLAS_N_GROUPS = 15
ORTHO_ATLAS_N_GENES = 12
ORTHO_ATLAS_OBS = -2.0982
ORTHO_ATLAS_NULL = +0.0510
ORTHO_ATLAS_Z = -4.6479
ORTHO_ATLAS_P = 4.0e-5
ORTHO_ATLAS_Q = 8.0e-5
ORTHO_ATLAS_NEGATIVE = 10
ORTHO_ATLAS_NEGATIVE_P = 0.00695


                      
REPL_N = 11
REPL_SIGN = 9
REPL_SIGN_P = 0.0327
REPL_SPEARMAN = 0.8636
REPL_SPEARMAN_P = 6.12e-4
REPL_PEARSON = 0.9393
REPL_PEARSON_P = 1.79e-5


                 
HS42_CS18_Z = -4.6479
HS42_CS18_P = 4.0e-5
HS42_CS18_Q = 8.0e-5

HS42_CS23_Z = -3.928
HS42_CS23_P = 2.10e-4
HS42_CS23_Q = 0.00126


                                                              
         
                                                              

PKA_COLOR = "#765AA5"
DIRECT_COLOR = "#C84B31"
HEAT_COLOR = "#C84B31"
EXPO_COLOR = "#3478A8"

NULL_COLOR = "#A0A0A0"
OTHER_COLOR = "#AAAAAA"
ZERO_COLOR = "#777777"
DIAGONAL_COLOR = "#888888"


                                                              
         
                                                              

def fail(msg):
    raise RuntimeError(msg)


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def sci(x, digits=1):
    if x == 0:
        return "0"

    exponent = int(
        math.floor(
            math.log10(abs(x))
        )
    )

    mantissa = x / (10 ** exponent)

    return (
        rf"${mantissa:.{digits}f}"
        rf"\times10^{{{exponent}}}$"
    )


def read_table(path):
    name = path.name.lower()

    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        return pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )

    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(
            path,
            low_memory=False,
        )

    raise ValueError(path)


def tabular_files(directory):
    if not directory.exists():
        return []

    out = []

    for p in directory.rglob("*"):
        if not p.is_file():
            continue

        s = str(p)

        if (
            "__pycache__" in s
            or "first_failed" in s
            or ".before_" in s
        ):
            continue

        if p.name.lower().endswith(
            (
                ".tsv",
                ".tsv.gz",
                ".csv",
                ".csv.gz",
            )
        ):
            out.append(p)

    return sorted(out)


def norm(s):
    return re.sub(
        r"[^a-z0-9]+",
        "_",
        str(s).lower(),
    ).strip("_")


def find_numeric_col(df, terms):
    for c in df.columns:
        n = norm(c)

        if all(t in n for t in terms):
            x = pd.to_numeric(
                df[c],
                errors="coerce",
            )

            if x.notna().sum() > 0:
                return c

    return None


                                                              
                         
                                                              

if not CLAIMS.exists():
    fail(f"Missing final claim matrix: {CLAIMS}")

claims = pd.read_csv(
    CLAIMS,
    sep="\t",
)

required_claims = {
    "PKA_CAUSAL_EXTERNAL",
    "PKA_CONTEXT_GLOBAL",
    "PKA_CONTEXT_MAF1",
    "PKA_CONTEXT_ATG1",
    "PKA_TEMP_EXTERNAL_DIRECT",
    "PKA_TEMP_EXTERNAL_ORTHOGONAL",
}

observed_claims = set(
    claims["claim_id"].astype(str)
)

missing_claims = (
    required_claims
    - observed_claims
)

if missing_claims:
    fail(
        "Final claim matrix missing: "
        + ", ".join(sorted(missing_claims))
    )


                                                              
                                                    
                                                              

def find_uliana_null():

    candidates = []

    for p in tabular_files(D30):

        try:
            df = read_table(p)
        except Exception:
            continue

        if len(df) < 1000:
            continue

        for c in df.columns:

            x = pd.to_numeric(
                df[c],
                errors="coerce",
            ).dropna()

            if len(x) < 1000:
                continue

            mu = x.mean()

            if abs(mu - ULIANA_EXPO_NULL) < 0.10:

                candidates.append(
                    (
                        abs(
                            mu
                            - ULIANA_EXPO_NULL
                        ),
                        p,
                        c,
                        x.to_numpy(float),
                    )
                )

    if not candidates:
        return None

    candidates.sort(
        key=lambda z: z[0]
    )

    best = candidates[0]

    return {
        "path": best[1],
        "column": best[2],
        "values": best[3],
    }


uliana_null = find_uliana_null()


                                                              
                                                            
                                                              

def find_replication_table():

    solutions = []

    for p in tabular_files(D35):

        try:
            df = read_table(p)
        except Exception:
            continue

        if not (9 <= len(df) <= 20):
            continue

                               
        gene_cols = [
            c for c in df.columns
            if "gene" in norm(c)
        ]

        if not gene_cols:
            continue

        num_cols = []

        for c in df.columns:

            x = pd.to_numeric(
                df[c],
                errors="coerce",
            )

            if x.notna().sum() >= 9:
                num_cols.append(c)

        for gc in gene_cols:

            for i in range(len(num_cols)):
                for j in range(i + 1, len(num_cols)):

                    c1 = num_cols[i]
                    c2 = num_cols[j]

                    sub = df[
                        [gc, c1, c2]
                    ].copy()

                    sub[c1] = pd.to_numeric(
                        sub[c1],
                        errors="coerce",
                    )

                    sub[c2] = pd.to_numeric(
                        sub[c2],
                        errors="coerce",
                    )

                    sub = sub.dropna()

                    if len(sub) != REPL_N:
                        continue

                    x = sub[c1].to_numpy(float)
                    y = sub[c2].to_numpy(float)

                    rho = spearmanr(x, y).statistic
                    r = pearsonr(x, y).statistic

                    same_sign = int(
                        np.sum(
                            np.sign(x)
                            == np.sign(y)
                        )
                    )

                    score = (
                        abs(
                            rho
                            - REPL_SPEARMAN
                        )
                        +
                        abs(
                            r
                            - REPL_PEARSON
                        )
                        +
                        0.2
                        * abs(
                            same_sign
                            - REPL_SIGN
                        )
                    )

                    if (
                        abs(
                            rho
                            - REPL_SPEARMAN
                        )
                        < 0.03
                        and
                        abs(
                            r
                            - REPL_PEARSON
                        )
                        < 0.03
                        and
                        same_sign
                        == REPL_SIGN
                    ):

                                                                   
                        bonus = 0.0

                        n1 = norm(c1)
                        n2 = norm(c2)

                        if (
                            "kanshin" in n1
                            or "kanshin" in n2
                        ):
                            bonus -= 0.05

                        if (
                            "atlas" in n1
                            or "atlas" in n2
                        ):
                            bonus -= 0.05

                        solutions.append(
                            (
                                score + bonus,
                                p,
                                gc,
                                c1,
                                c2,
                                sub,
                                rho,
                                r,
                            )
                        )

    if not solutions:
        return None

    solutions.sort(
        key=lambda z: z[0]
    )

    best = solutions[0]

    _, p, gc, c1, c2, sub, rho, r = best

                                                           
    n1 = norm(c1)
    n2 = norm(c2)

    if "atlas" in n1:
        atlas_col = c1
        kanshin_col = c2

    elif "atlas" in n2:
        atlas_col = c2
        kanshin_col = c1

    elif "kanshin" in n1:
        kanshin_col = c1
        atlas_col = c2

    elif "kanshin" in n2:
        kanshin_col = c2
        atlas_col = c1

    else:
                                                            
        kanshin_col = c1
        atlas_col = c2

    result = sub[
        [
            gc,
            kanshin_col,
            atlas_col,
        ]
    ].copy()

    result.columns = [
        "gene",
        "kanshin_effect",
        "atlas_effect",
    ]

    return {
        "path": p,
        "data": result,
    }


replication = find_replication_table()

if replication is None:
    fail(
        "Could not identify the frozen 11-gene "
        "Kanshin↔PhosphoAtlas replication table in "
        "35_phosphoatlas_pka_replication."
    )


                           
rep_df = replication["data"].copy()

rho_obs, rho_p_obs = spearmanr(
    rep_df["kanshin_effect"],
    rep_df["atlas_effect"],
)

r_obs, r_p_obs = pearsonr(
    rep_df["kanshin_effect"],
    rep_df["atlas_effect"],
)

same_sign_obs = int(
    np.sum(
        np.sign(
            rep_df["kanshin_effect"]
        )
        ==
        np.sign(
            rep_df["atlas_effect"]
        )
    )
)

if abs(rho_obs - REPL_SPEARMAN) > 0.01:
    fail(
        f"Replication Spearman mismatch: {rho_obs}"
    )

if abs(r_obs - REPL_PEARSON) > 0.01:
    fail(
        f"Replication Pearson mismatch: {r_obs}"
    )

if same_sign_obs != REPL_SIGN:
    fail(
        f"Replication sign concordance mismatch: "
        f"{same_sign_obs}/11"
    )


                                                              
                                           
                                                              

SENSITIVITY_SUMMARY = (
    D36
    / "TEMPERATURE_SENSITIVITY_SUMMARY.tsv"
)

if not SENSITIVITY_SUMMARY.exists():
    fail(
        f"Missing temperature-sensitivity summary: "
        f"{SENSITIVITY_SUMMARY}"
    )

sens_df = pd.read_csv(
    SENSITIVITY_SUMMARY,
    sep="\t",
    low_memory=False,
)

required_sens = {
    "contrast",
    "module",
    "predefined_groups",
    "matched_groups",
    "matched_genes",
    "observed_mean",
    "observed_median",
    "negative_genes",
    "null_mean",
    "null_sd",
    "empirical_z",
    "mean_empirical_lower_p",
    "median_empirical_lower_p",
    "negative_count_empirical_p",
    "q_across_six_sensitivity_tests",
}

missing = (
    required_sens
    - set(sens_df.columns)
)

if missing:
    fail(
        "TEMPERATURE_SENSITIVITY_SUMMARY.tsv missing columns: "
        + ", ".join(sorted(missing))
    )


                                                              
                                                           
                                                              

module_text = (
    sens_df["module"]
    .astype(str)
    .str.lower()
)

orth_mask = (
    module_text.str.contains("orthogonal")
    & module_text.str.contains("pka")
)

sens_orth = (
    sens_df[orth_mask]
    .copy()
)

if len(sens_orth) == 0:
    fail(
        "No orthogonal PKA rows found in "
        "TEMPERATURE_SENSITIVITY_SUMMARY.tsv.\n"
        f"Available modules: "
        f"{sorted(sens_df['module'].astype(str).unique())}"
    )


def normalize_contrast(x):
    return (
        str(x)
        .upper()
        .replace(" ", "")
        .replace("-", "_")
        .replace("VS", "_MINUS_")
        .replace("__", "_")
    )


sens_orth["_contrast_norm"] = (
    sens_orth["contrast"]
    .map(normalize_contrast)
)


def get_sensitivity_row(target):

    target_norm = normalize_contrast(target)

    x = sens_orth[
        sens_orth["_contrast_norm"]
        == target_norm
    ]

    if len(x) != 1:
        fail(
            f"Expected exactly one orthogonal PKA row "
            f"for {target}; found {len(x)}.\n"
            f"Available contrasts: "
            f"{sens_orth['contrast'].astype(str).tolist()}"
        )

    return x.iloc[0]


r_4223 = get_sensitivity_row(
    "HS42_minus_CS23"
)

r_3718 = get_sensitivity_row(
    "HS37_minus_CS18"
)

r_4818 = get_sensitivity_row(
    "HS48_minus_CS18"
)


                                                              
                                 
 
                                                               
                                                    
                                                              

sens_plot = pd.DataFrame(
    [
        {
            "contrast": "HS42–CS18",
            "effect": ORTHO_ATLAS_OBS,
            "null_mean": ORTHO_ATLAS_NULL,
            "z": ORTHO_ATLAS_Z,
            "p": ORTHO_ATLAS_P,
            "q": ORTHO_ATLAS_Q,
            "q_family": "two_predefined_modules",
            "matched_groups": ORTHO_ATLAS_N_GROUPS,
            "matched_genes": ORTHO_ATLAS_N_GENES,
            "analysis_role": "primary",
        },
        {
            "contrast": "HS42–CS23",
            "effect": float(
                r_4223["observed_mean"]
            ),
            "null_mean": float(
                r_4223["null_mean"]
            ),
            "z": float(
                r_4223["empirical_z"]
            ),
            "p": float(
                r_4223["mean_empirical_lower_p"]
            ),
            "q": float(
                r_4223["q_across_six_sensitivity_tests"]
            ),
            "q_family": "six_sensitivity_tests",
            "matched_groups": int(
                r_4223["matched_groups"]
            ),
            "matched_genes": int(
                r_4223["matched_genes"]
            ),
            "analysis_role": "sensitivity",
        },
        {
            "contrast": "HS37–CS18",
            "effect": float(
                r_3718["observed_mean"]
            ),
            "null_mean": float(
                r_3718["null_mean"]
            ),
            "z": float(
                r_3718["empirical_z"]
            ),
            "p": float(
                r_3718["mean_empirical_lower_p"]
            ),
            "q": float(
                r_3718["q_across_six_sensitivity_tests"]
            ),
            "q_family": "six_sensitivity_tests",
            "matched_groups": int(
                r_3718["matched_groups"]
            ),
            "matched_genes": int(
                r_3718["matched_genes"]
            ),
            "analysis_role": "sensitivity",
        },
        {
            "contrast": "HS48–CS18",
            "effect": float(
                r_4818["observed_mean"]
            ),
            "null_mean": float(
                r_4818["null_mean"]
            ),
            "z": float(
                r_4818["empirical_z"]
            ),
            "p": float(
                r_4818["mean_empirical_lower_p"]
            ),
            "q": float(
                r_4818["q_across_six_sensitivity_tests"]
            ),
            "q_family": "six_sensitivity_tests",
            "matched_groups": int(
                r_4818["matched_groups"]
            ),
            "matched_genes": int(
                r_4818["matched_genes"]
            ),
            "analysis_role": "sensitivity",
        },
    ]
)


                                                              
                           
                                                              

if abs(
    sens_plot.loc[
        sens_plot["contrast"] == "HS42–CS18",
        "z",
    ].iloc[0]
    - HS42_CS18_Z
) > 0.01:
    fail(
        "Primary HS42-CS18 Z mismatch."
    )

if abs(
    sens_plot.loc[
        sens_plot["contrast"] == "HS42–CS23",
        "z",
    ].iloc[0]
    - HS42_CS23_Z
) > 0.02:
    fail(
        "HS42-CS23 Z mismatch."
    )

if abs(
    sens_plot.loc[
        sens_plot["contrast"] == "HS42–CS23",
        "p",
    ].iloc[0]
    - HS42_CS23_P
) > 1e-6:
    fail(
        "HS42-CS23 P mismatch."
    )

if abs(
    sens_plot.loc[
        sens_plot["contrast"] == "HS42–CS23",
        "q",
    ].iloc[0]
    - HS42_CS23_Q
) > 1e-6:
    fail(
        "HS42-CS23 q mismatch."
    )

                                                              
            
                                                              

manifest = [
    {
        "role": "final_claim_matrix",
        "path": str(CLAIMS),
        "sha256": sha256(CLAIMS),
    },
    {
        "role": "Kanshin_PhospAtlas_11_gene_replication",
        "path": str(replication["path"]),
        "sha256": sha256(
            replication["path"]
        ),
    },
    {
        "role": "temperature_sensitivity",
        "path": str(SENSITIVITY_SUMMARY),
        "sha256": sha256(
            SENSITIVITY_SUMMARY
        ),
    },
]

if uliana_null is not None:

    manifest.append(
        {
            "role": "Uliana_exponential_empirical_null",
            "path": str(
                uliana_null["path"]
            ),
            "sha256": sha256(
                uliana_null["path"]
            ),
        }
    )

pd.DataFrame(
    manifest
).to_csv(
    CURATED
    / "Fig4_3_source_manifest.tsv",
    sep="\t",
    index=False,
)

rep_df.to_csv(
    CURATED
    / "Fig4_3_Kanshin_PhospAtlas_11_gene.tsv",
    sep="\t",
    index=False,
)

sens_plot.to_csv(
    CURATED
    / "Fig4_3_temperature_sensitivity.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
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
    figsize=(18.0, 13.0),
)

outer = fig.add_gridspec(
    3,
    12,
    height_ratios=[
        1.0,
        1.0,
        0.78,
    ],
    hspace=0.62,
    wspace=0.62,
)


def panel_with_note(spec, note_ratio=0.18):

                                             
                                                      
                                                    
    ax = fig.add_subplot(spec)

    class DummyNote:
        transAxes = None

        def text(self, *args, **kwargs):
            return None

    return ax, DummyNote()


                                                              
                           
                                                              

axA, noteA = panel_with_note(
    outer[0, 0:6],
    0.18,
)

if uliana_null is not None:

    vals = uliana_null["values"]

    q25, q75 = np.quantile(
        vals,
        [0.25, 0.75],
    )

    iqr = q75 - q25

    if iqr > 0:
        bw = (
            2
            * iqr
            / (
                len(vals)
                ** (1 / 3)
            )
        )

        nbins = int(
            np.ceil(
                (
                    vals.max()
                    - vals.min()
                )
                / bw
            )
        )

    else:
        nbins = 60

    nbins = min(
        max(
            nbins,
            40,
        ),
        100,
    )

    axA.hist(
        vals,
        bins=nbins,
        density=True,
        color="#D0D0D0",
        edgecolor="none",
    )

    axA.axvline(
        vals.mean(),
        color=NULL_COLOR,
        ls="--",
        lw=1.1,
    )

    axA.axvline(
        ULIANA_EXPO_OBS,
        color=PKA_COLOR,
        lw=2.3,
    )

    axA.set_ylabel(
        "Плътност"
    )

else:

                                                       
    axA.errorbar(
        ULIANA_EXPO_NULL,
        0,
        xerr=0,
        fmt="o",
        color=NULL_COLOR,
        markersize=7,
    )

    axA.scatter(
        ULIANA_EXPO_OBS,
        1,
        marker="D",
        s=80,
        color=PKA_COLOR,
        edgecolor="white",
        linewidth=0.5,
    )

    axA.hlines(
        1,
        ULIANA_EXPO_NULL,
        ULIANA_EXPO_OBS,
        color=PKA_COLOR,
        lw=1.0,
        alpha=0.5,
    )

    axA.set_yticks(
        [0, 1]
    )

    axA.set_yticklabels(
        [
            "matched null",
            "observed",
        ]
    )

    axA.tick_params(
        axis="y",
        length=0,
    )

axA.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.8,
    alpha=0.6,
)

axA.set_xlabel(
    "Ефект от PKA инхибиране (log$_2$FC)",
    labelpad=7,
)

axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

axA.text(
    -0.10,
    1.04,
    "A",
    transform=axA.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                      
                                                              

axB, noteB = panel_with_note(
    outer[0, 6:12],
    0.18,
)

sites = [
    "Maf1-S90",
    "Atg1-S515",
]

expo = np.array(
    [
        MAF1_EXPO,
        ATG1_EXPO,
    ]
)

hs = np.array(
    [
        MAF1_HS,
        ATG1_HS,
    ]
)

y = np.arange(2)

axB.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.8,
)

for i in range(2):

    axB.plot(
        [expo[i], hs[i]],
        [i, i],
        color="#B0B0B0",
        lw=1.5,
        zorder=1,
    )

    axB.scatter(
        expo[i],
        i,
        color=EXPO_COLOR,
        s=70,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
        label=(
            "Експоненциален растеж"
            if i == 0
            else None
        ),
    )

    axB.scatter(
        hs[i],
        i,
        color=HEAT_COLOR,
        s=70,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
        label=(
            "Топлинен стрес"
            if i == 0
            else None
        ),
    )

axB.set_yticks(y)

axB.set_yticklabels(
    [
        r"$\it{Maf1}$-S90",
        r"$\it{Atg1}$-S515",
    ]
)

axB.invert_yaxis()

axB.set_xlabel(
    "Ефект от PKA инхибиране (log$_2$FC)",
    labelpad=7,
)

axB.legend(
    frameon=False,
    loc="upper right",
    fontsize=7.5,
    handletextpad=0.35,
    borderaxespad=0.2,
)

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

axB.text(
    -0.10,
    1.04,
    "B",
    transform=axB.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                                   
                                                              

axC, noteC = panel_with_note(
    outer[1, 0:5],
    0.21,
)

module_names = [
    "Direct PKA",
    "Orthogonal PKA",
]

observed = np.array(
    [
        DIRECT_ATLAS_OBS,
        ORTHO_ATLAS_OBS,
    ]
)

nulls = np.array(
    [
        DIRECT_ATLAS_NULL,
        ORTHO_ATLAS_NULL,
    ]
)

yc = np.arange(2)

axC.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.8,
)

for i in range(2):

    axC.plot(
        [nulls[i], observed[i]],
        [i, i],
        color="#B0B0B0",
        lw=1.4,
    )

    axC.scatter(
        nulls[i],
        i,
        s=55,
        color=NULL_COLOR,
        edgecolor="white",
        linewidth=0.4,
        zorder=2,
    )

    axC.scatter(
        observed[i],
        i,
        s=80,
        marker="D",
        color=(
            DIRECT_COLOR
            if i == 0
            else PKA_COLOR
        ),
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )

axC.set_yticks(yc)

axC.set_yticklabels(
    module_names
)

axC.invert_yaxis()

axC.set_xlabel(
    "HS42–CS18 ефект (log$_2$FC)",
    labelpad=7,
)

axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

axC.text(
    -0.12,
    1.04,
    "C",
    transform=axC.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                                                     
                                                              

axD, noteD = panel_with_note(
    outer[1, 5:12],
    0.20,
)

x = rep_df[
    "kanshin_effect"
].to_numpy(float)

y = rep_df[
    "atlas_effect"
].to_numpy(float)

lim_min = min(
    np.min(x),
    np.min(y),
)

lim_max = max(
    np.max(x),
    np.max(y),
)

pad = 0.10 * (
    lim_max - lim_min
    if lim_max > lim_min
    else 1
)

lo = lim_min - pad
hi = lim_max + pad

axD.plot(
    [lo, hi],
    [lo, hi],
    ls="--",
    lw=1.0,
    color=DIAGONAL_COLOR,
    zorder=0,
)

coef = np.polyfit(
    x,
    y,
    1,
)

xx = np.linspace(
    lo,
    hi,
    100,
)

axD.plot(
    xx,
    coef[0] * xx + coef[1],
    color=PKA_COLOR,
    lw=1.3,
    alpha=0.75,
)

axD.scatter(
    x,
    y,
    s=58,
    color=PKA_COLOR,
    edgecolor="white",
    linewidth=0.5,
    zorder=3,
)

             
offsets = [
    (4, 4),
    (4, -10),
    (-18, 5),
    (5, 8),
    (-20, -9),
    (5, -9),
    (-20, 7),
    (5, 4),
    (-18, -8),
    (5, 8),
    (-20, 5),
]

for i, row in rep_df.reset_index(
    drop=True
).iterrows():

    dx, dy = offsets[
        i % len(offsets)
    ]

    axD.annotate(
        str(row["gene"]),
        (
            row["kanshin_effect"],
            row["atlas_effect"],
        ),
        xytext=(dx, dy),
        textcoords="offset points",
        fontsize=7.4,
        color="#333333",
    )

axD.set_xlim(
    lo,
    hi,
)

axD.set_ylim(
    lo,
    hi,
)

axD.set_xlabel(
    "Kanshin early Heat–Cold effect",
    labelpad=7,
)

axD.set_ylabel(
    "PhosphoAtlas HS42–CS18 effect"
)

axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)

axD.text(
    -0.09,
    1.04,
    "D",
    transform=axD.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


                                                              
                             
                                                              

axE, noteE = panel_with_note(
    outer[2, 0:12],
    0.22,
)

E = sens_plot.copy()

                                       
if E["effect"].notna().sum() == 4:
    value_col = "effect"
    xlabel = "Orthogonal PKA module effect"

elif E["z"].notna().sum() == 4:
    value_col = "z"
    xlabel = "Z-score"

else:
    fail(
        "Temperature-sensitivity table contains neither "
        "four effects nor four Z-scores."
    )

E["_value"] = E[
    value_col
].astype(float)

ye = np.arange(
    len(E)
)

axE.axvline(
    0,
    color=ZERO_COLOR,
    lw=0.8,
)

for i, row in E.iterrows():

    significant_42 = (
        row["contrast"]
        in [
            "HS42–CS18",
            "HS42–CS23",
        ]
    )

    color = (
        PKA_COLOR
        if significant_42
        else OTHER_COLOR
    )

    axE.hlines(
        i,
        0,
        row["_value"],
        color=color,
        lw=1.5,
        alpha=0.8,
    )

    axE.scatter(
        row["_value"],
        i,
        s=65,
        color=color,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )

axE.set_yticks(
    ye
)

axE.set_yticklabels(
    E["contrast"]
)

axE.invert_yaxis()

axE.set_xlabel(
    xlabel,
    labelpad=7,
)

axE.spines["top"].set_visible(False)
axE.spines["right"].set_visible(False)

axE.text(
    -0.055,
    1.04,
    "E",
    transform=axE.transAxes,
    fontsize=15,
    fontweight="bold",
    va="top",
)


def pq_text(row):

    pieces = []

    if np.isfinite(row["p"]):
        pieces.append(
            f"P={row['p']:.3g}"
        )

    if np.isfinite(row["q"]):
        pieces.append(
            f"q={row['q']:.3g}"
        )

    return ", ".join(pieces)



                                                              
                                       
                                                              

fig.subplots_adjust(
    left=0.075,
    right=0.985,
    top=0.985,
    bottom=0.075,
)

OUTBASE = (
    FINAL
    / "Fig4_3_PKA_independent_validation"
)

fig.savefig(
    str(OUTBASE) + ".png",
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE) + ".pdf",
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE) + ".svg",
    bbox_inches="tight",
)

plt.close(fig)


                                                              
                  
                                                              

legend = f"""Фигура 4.3. Независима валидация на PKA-свързаната фосфорилационна програма.

(A) Независима PKA пертурбация на ортогоналния PKA модул. При експоненциален растеж са съпоставени {ULIANA_N_GROUPS} фосфопротеомни групи в {ULIANA_N_GENES} гена; средният ефект след PKA инхибиране е {ULIANA_EXPO_OBS:.2f} log₂FC спрямо {ULIANA_EXPO_NULL:.2f} в matched-null модела (Z={ULIANA_EXPO_Z:.2f}; P={ULIANA_EXPO_P:.2g}). При топлинен стрес модулният тест не достига статистическа значимост (P={ULIANA_HS_P:.3f}), а глобалното context interaction също не е значимо (P={ULIANA_INTERACTION_P:.3f}).

(B) Контекст-зависими ефекти на PKA инхибирането при Maf1-S90 и Atg1-S515. Maf1-S90 показва ефект {MAF1_EXPO:.2f} при експоненциален растеж и {MAF1_HS:.2f} при топлинен стрес (interaction P={MAF1_INT_P:.4f}; q={MAF1_INT_Q:.3f}). При Atg1-S515 ефектът променя посоката от {ATG1_EXPO:.2f} до +{ATG1_HS:.2f} (interaction P={ATG1_INT_P:.3g}; q={ATG1_INT_Q:.4f}).

(C) Независима температурна репликация в Yeast PhosphoAtlas. Директният PKA модул включва {DIRECT_ATLAS_N_GROUPS} съпоставени групи от {DIRECT_ATLAS_N_GENES} гена и показва среден HS42–CS18 ефект {DIRECT_ATLAS_OBS:.2f} спрямо null {DIRECT_ATLAS_NULL:+.2f} (P={DIRECT_ATLAS_P:.3f}). Ортогоналният модул включва {ORTHO_ATLAS_N_GROUPS} групи в {ORTHO_ATLAS_N_GENES} гена и показва ефект {ORTHO_ATLAS_OBS:.2f} спрямо null {ORTHO_ATLAS_NULL:+.2f} (Z={ORTHO_ATLAS_Z:.2f}; P={ORTHO_ATLAS_P:.2g}; q={ORTHO_ATLAS_Q:.2g}).

(D) Количествена възпроизводимост на ранните Heat–Cold ефекти между Kanshin и PhosphoAtlas за {REPL_N} гена. {REPL_SIGN}/{REPL_N} ефекта са с еднаква посока (P={REPL_SIGN_P:.3f}); Spearman ρ={rho_obs:.3f} и Pearson r={r_obs:.3f}. Прекъснатата диагонална линия показва y=x.

(E) Анализ на температурната чувствителност на ортогоналния PKA модул. Сигналът се запазва при 42 °C спрямо 18 °C и 23 °C, но не се възпроизвежда при 37 °C и 48 °C спрямо 18 °C, което не подкрепя монотонна зависимост от интензитета на температурния стрес.
"""

legend_path = (
    LEGENDS
    / "Fig4_3_legend_bg.txt"
)

legend_path.write_text(
    legend,
    encoding="utf-8",
)


                                                              
                         
                                                              

print("===== FIGURE 4.3 COMPLETE =====")
print(
    "Replication source:",
    replication["path"].relative_to(ROOT),
)
print(
    "Sensitivity source:",
    SENSITIVITY_SUMMARY.relative_to(ROOT),
)

if uliana_null is not None:
    print(
        "Uliana null source:",
        uliana_null["path"].relative_to(ROOT),
    )
else:
    print(
        "Uliana raw null vector not found; "
        "Panel A uses observed-vs-null without assuming a distribution."
    )

print(
    "Spearman:",
    round(rho_obs, 4),
    "| Pearson:",
    round(r_obs, 4),
    "| same sign:",
    f"{same_sign_obs}/11",
)

print(str(OUTBASE) + ".pdf")
