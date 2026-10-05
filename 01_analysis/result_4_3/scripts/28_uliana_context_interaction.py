#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import t as tdist


ROOT = Path("/")

ULIANA = (
    ROOT
    / "29_external_validation"
    / "PXD052971"
    / "raw"
    / "Uliana_2026_Table_S2.xlsx"
)

S1FILE = (
    ROOT
    / "00_raw"
    / "Kanshin_2015_Table_S1.xlsx"
)

DIRECT_A = (
    ROOT
    / "07_regulator_validation"
    / "FINAL_DIRECT_A_kinase_edges.tsv"
)

KANSHIN_AUC = (
    ROOT
    / "26_pka_global_trajectory"
    / "ALL_S1_heat_cold_AUC.tsv.gz"
)

OUT = (
    ROOT
    / "31_uliana_context_interaction"
)
OUT.mkdir(parents=True, exist_ok=True)


                                                              
         
                                                              

def norm_id(x):
    if pd.isna(x):
        return ""

    s = str(x).strip()

    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except Exception:
        pass

    return s


def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def clean_gene_series(x):
    return (
        x.where(x.notna(), "")
        .astype(str)
        .str.strip()
        .replace({
            "nan": "",
            "NaN": "",
            "None": "",
        })
    )


def split_sites(x):
    s = clean(x)

    if not s:
        return []

    return [
        z.strip()
        for z in s.split("/")
        if z.strip()
    ]


def bh(p):
    p = np.asarray(
        p,
        dtype=float,
    )

    out = np.full(
        len(p),
        np.nan,
    )

    good = np.isfinite(p)

    if not good.any():
        return out

    pv = p[good]

    order = np.argsort(pv)
    ranked = pv[order]

    m = len(ranked)

    q = (
        ranked
        * m
        / np.arange(1, m + 1)
    )

    q = np.minimum.accumulate(
        q[::-1]
    )[::-1]

    q = np.minimum(
        q,
        1.0,
    )

    restored = np.empty_like(q)
    restored[order] = q

    out[good] = restored

    return out


def parse_uliana_id(x):
    s = clean(x)

    m = re.search(
        r"_([STY])_(\d+)$",
        s,
    )

    if not m:
        return "", np.nan

    return (
        m.group(1),
        int(m.group(2)),
    )


def parse_targeted_id(x):
    s = clean(x)

    m = re.match(
        r"^(.+)_([STY])(\d+)(?:_.*)?$",
        s,
    )

    if not m:
        return "", "", np.nan

    return (
        m.group(1),
        m.group(2),
        int(m.group(3)),
    )


def interaction_model(
    expo_ctrl,
    expo_pkai,
    hs_ctrl,
    hs_pkai,
):
    """
    y = b0 + b1*HS + b2*PKAi + b3*(HS x PKAi)

    b3 =
      (HS_PKAi - HS_ctrl)
      -
      (Expo_PKAi - Expo_ctrl)
    """

    groups = [
        (expo_ctrl, 0, 0),
        (expo_pkai, 0, 1),
        (hs_ctrl,   1, 0),
        (hs_pkai,   1, 1),
    ]

    yy = []
    xx = []

    for values, hs, pkai in groups:

        for v in values:

            if not np.isfinite(v):
                continue

            yy.append(v)

            xx.append([
                1.0,
                float(hs),
                float(pkai),
                float(hs * pkai),
            ])

    y = np.asarray(
        yy,
        dtype=float,
    )

    X = np.asarray(
        xx,
        dtype=float,
    )

    if len(y) <= X.shape[1]:
        return {
            "interaction_beta": np.nan,
            "interaction_se": np.nan,
            "interaction_t": np.nan,
            "interaction_p": np.nan,
            "model_df_resid": np.nan,
        }

    if np.linalg.matrix_rank(X) < X.shape[1]:
        return {
            "interaction_beta": np.nan,
            "interaction_se": np.nan,
            "interaction_t": np.nan,
            "interaction_p": np.nan,
            "model_df_resid": np.nan,
        }

    beta = np.linalg.lstsq(
        X,
        y,
        rcond=None,
    )[0]

    resid = y - X @ beta

    df = len(y) - X.shape[1]

    s2 = float(
        np.sum(resid ** 2)
        / df
    )

    cov = (
        s2
        * np.linalg.inv(
            X.T @ X
        )
    )

    se = float(
        np.sqrt(
            cov[3, 3]
        )
    )

    if se == 0:
        tval = np.nan
        pval = np.nan
    else:
        tval = float(
            beta[3] / se
        )

        pval = float(
            2.0
            * tdist.sf(
                abs(tval),
                df,
            )
        )

    return {
        "interaction_beta":
            float(beta[3]),

        "interaction_se":
            se,

        "interaction_t":
            tval,

        "interaction_p":
            pval,

        "model_df_resid":
            df,
    }


def target_restricted_unique(
    df,
    targets,
    label,
):
    """
    Remove unmappable rows, restrict to keys that could actually
    match our target set, then require uniqueness.
    """

    d = df.copy()

    d = d[
        d["gene"].notna()
        & (d["gene"] != "")
        & d["aa"].isin(["S", "T", "Y"])
        & pd.to_numeric(
            d["position"],
            errors="coerce",
        ).notna()
    ].copy()

    d["position"] = pd.to_numeric(
        d["position"],
        errors="raise",
    ).astype(int)

    target_keys = (
        targets[
            [
                "gene",
                "aa",
                "position",
            ]
        ]
        .drop_duplicates()
    )

    d = d.merge(
        target_keys,
        on=[
            "gene",
            "aa",
            "position",
        ],
        how="inner",
    )

    dup = d[
        d.duplicated(
            [
                "gene",
                "aa",
                "position",
            ],
            keep=False,
        )
    ].sort_values(
        [
            "gene",
            "aa",
            "position",
        ]
    )

    if len(dup):
        print()
        print(
            f"DUPLICATE TARGETABLE ROWS IN {label}"
        )
        print(
            dup.to_string(
                index=False
            )
        )

        raise RuntimeError(
            f"{label}: duplicate exact-site keys "
            "among our target sites"
        )

    return d


                                                              
                                   
                                                              

edges = pd.read_csv(
    DIRECT_A,
    sep="\t",
)

s1 = pd.read_excel(
    S1FILE,
    sheet_name="Phosphopeptides",
)

auc = pd.read_csv(
    KANSHIN_AUC,
    sep="\t",
    compression="gzip",
)


edges["s1_id"] = (
    edges["s1_id"]
    .map(norm_id)
)

s1["s1_id"] = (
    s1["id"]
    .map(norm_id)
)

auc["s1_id"] = (
    auc["s1_id"]
    .map(norm_id)
)


pka_edges = (
    edges[
        edges["regulator"] == "TPK1"
    ][
        ["s1_id"]
    ]
    .drop_duplicates()
)


pka_groups = pka_edges.merge(
    s1,
    on="s1_id",
    how="left",
    validate="1:1",
)


gene_col = (
    "Gene"
    if "Gene" in pka_groups.columns
    else "Standard Name"
)


target_rows = []

for _, r in pka_groups.iterrows():

    gene = clean(
        r[gene_col]
    )

    for site in split_sites(
        r["pSites"]
    ):

        m = re.fullmatch(
            r"([STY])(\d+)",
            site,
        )

        if not m:
            raise RuntimeError(
                f"Cannot parse {gene} {site}"
            )

        target_rows.append({
            "s1_id":
                r["s1_id"],

            "gene":
                gene,

            "site":
                site,

            "aa":
                m.group(1),

            "position":
                int(
                    m.group(2)
                ),
        })


targets = pd.DataFrame(
    target_rows
)


if len(targets) != 10:
    raise RuntimeError(
        f"Expected 10 direct-A exact sites; "
        f"got {len(targets)}"
    )


if targets[
    [
        "gene",
        "site",
    ]
].duplicated().any():

    raise RuntimeError(
        "Duplicate exact direct-A sites"
    )


                                                              
                                                           
                                                              

auc_group = auc[
    [
        "s1_id",
        "gene",
        "pSites",
        "heat_minus_cold_auc",
    ]
].copy()


if auc_group[
    "s1_id"
].duplicated().any():

    dup = auc_group[
        auc_group[
            "s1_id"
        ].duplicated(
            keep=False
        )
    ]

    print(
        dup.to_string(
            index=False
        )
    )

    raise RuntimeError(
        "Unexpected duplicate s1_id in Kanshin AUC table"
    )


auc_group = auc_group.rename(
    columns={
        "gene":
            "Kanshin_gene",

        "pSites":
            "Kanshin_group_pSites",

        "heat_minus_cold_auc":
            "Kanshin_heat_minus_cold_AUC",
    }
)


targets = targets.merge(
    auc_group,
    on="s1_id",
    how="left",
    validate="m:1",
)


targets.to_csv(
    OUT
    / "DIRECT_A_exact_targets.tsv",
    sep="\t",
    index=False,
)


                                                              
                                        
                                                              

stat = pd.read_excel(
    ULIANA,
    sheet_name="Statistical analysis phospho",
)


parsed = stat[
    "ID"
].map(
    parse_uliana_id
)


stat["aa"] = [
    x[0]
    for x in parsed
]

stat["position"] = [
    x[1]
    for x in parsed
]

stat["gene"] = clean_gene_series(
    stat[
        "Gene name"
    ]
)


st = stat[
    [
        "gene",
        "aa",
        "position",
        "ID",

        "Expo_PKAi-Expo_ctrl",
        "HS_PKAi-HS_ctrl",

        "pval_adj_BH_Expo_PKAi-Expo_ctrl",
        "pval_adj_BH_HS_PKAi-HS_ctrl",
    ]
].rename(
    columns={
        "ID":
            "Uliana_ID",

        "Expo_PKAi-Expo_ctrl":
            "Uliana_Expo_effect",

        "HS_PKAi-HS_ctrl":
            "Uliana_HS_effect",

        "pval_adj_BH_Expo_PKAi-Expo_ctrl":
            "Uliana_Expo_q",

        "pval_adj_BH_HS_PKAi-HS_ctrl":
            "Uliana_HS_q",
    }
)


st = target_restricted_unique(
    st,
    targets,
    "Uliana statistical analysis phospho",
)


print(
    "Untargeted exact target rows:",
    len(st),
)


res = targets.merge(
    st,
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="left",
    validate="1:1",
)


res[
    "observed_untargeted"
] = res[
    "Uliana_ID"
].notna()


                                                              
                                                 
                                                              

ints = pd.read_excel(
    ULIANA,
    sheet_name="Phosphosite intensities",
)


parsed = ints[
    "ID"
].map(
    parse_uliana_id
)


ints["aa"] = [
    x[0]
    for x in parsed
]

ints["position"] = [
    x[1]
    for x in parsed
]

ints["gene"] = clean_gene_series(
    ints[
        "gene_name"
    ]
)


ints = target_restricted_unique(
    ints,
    targets,
    "Uliana phosphosite intensities",
)


print(
    "Untargeted intensity exact target rows:",
    len(ints),
)


interaction_rows = []


for _, r in ints.iterrows():

    expo_ctrl = np.asarray(
        [
            r["Expo_ctrl_1"],
            r["Expo_ctrl_2"],
            r["Expo_ctrl_3"],
        ],
        dtype=float,
    )

    expo_pkai = np.asarray(
        [
            r["Expo_PKAi_1"],
            r["Expo_PKAi_2"],
            r["Expo_PKAi_3"],
        ],
        dtype=float,
    )

    hs_ctrl = np.asarray(
        [
            r["HS_ctrl_1"],
            r["HS_ctrl_2"],
            r["HS_ctrl_3"],
        ],
        dtype=float,
    )

    hs_pkai = np.asarray(
        [
            r["HS_PKAi_1"],
            r["HS_PKAi_2"],
            r["HS_PKAi_3"],
        ],
        dtype=float,
    )


    model = interaction_model(
        expo_ctrl,
        expo_pkai,
        hs_ctrl,
        hs_pkai,
    )


    expo_effect = float(
        np.nanmean(
            expo_pkai
        )
        - np.nanmean(
            expo_ctrl
        )
    )

    hs_effect = float(
        np.nanmean(
            hs_pkai
        )
        - np.nanmean(
            hs_ctrl
        )
    )


    interaction_rows.append({
        "gene":
            r["gene"],

        "aa":
            r["aa"],

        "position":
            int(
                r["position"]
            ),

        "Uliana_intensity_ID":
            r["ID"],

        "Expo_ctrl_mean":
            float(
                np.nanmean(
                    expo_ctrl
                )
            ),

        "Expo_PKAi_mean":
            float(
                np.nanmean(
                    expo_pkai
                )
            ),

        "HS_ctrl_mean":
            float(
                np.nanmean(
                    hs_ctrl
                )
            ),

        "HS_PKAi_mean":
            float(
                np.nanmean(
                    hs_pkai
                )
            ),

        "Expo_PKAi_minus_ctrl_replicates":
            expo_effect,

        "HS_PKAi_minus_ctrl_replicates":
            hs_effect,

        "HSctrl_minus_Expoctrl":
            float(
                np.nanmean(
                    hs_ctrl
                )
                - np.nanmean(
                    expo_ctrl
                )
            ),

        "abs_PKAi_effect_attenuation_in_HS":
            (
                abs(
                    expo_effect
                )
                - abs(
                    hs_effect
                )
            ),

        **model,
    })


inter = pd.DataFrame(
    interaction_rows
)


res = res.merge(
    inter,
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="left",
    validate="1:1",
)


observed_mask = (
    res[
        "observed_untargeted"
    ]
    & res[
        "interaction_p"
    ].notna()
)


res[
    "interaction_q_directA_overlap"
] = np.nan


res.loc[
    observed_mask,
    "interaction_q_directA_overlap",
] = bh(
    res.loc[
        observed_mask,
        "interaction_p",
    ]
)


                                                              
                                   
                                                              

tar = pd.read_excel(
    ULIANA,
    sheet_name="Targeted phosphosite int norm",
)


parsed = tar[
    "ID"
].map(
    parse_targeted_id
)


tar["gene"] = [
    x[0]
    for x in parsed
]

tar["aa"] = [
    x[1]
    for x in parsed
]

tar["position"] = [
    x[2]
    for x in parsed
]


tar_rows = []


for _, r in tar.iterrows():

    if (
        not r["gene"]
        or r["aa"] not in [
            "S",
            "T",
            "Y",
        ]
        or pd.isna(
            r["position"]
        )
    ):
        continue


    expo_ctrl = np.asarray(
        [
            r["Expo_ctrl_1"],
            r["Expo_ctrl_2"],
            r["Expo_ctrl_3"],
        ],
        dtype=float,
    )

    expo_pkai = np.asarray(
        [
            r["Expo_PKAi_1"],
            r["Expo_PKAi_2"],
            r["Expo_PKAi_3"],
        ],
        dtype=float,
    )

    hs_ctrl = np.asarray(
        [
            r["HS_ctrl_1"],
            r["HS_ctrl_2"],
            r["HS_ctrl_3"],
        ],
        dtype=float,
    )

    hs_pkai = np.asarray(
        [
            r["HS_PKAi_1"],
            r["HS_PKAi_2"],
            r["HS_PKAi_3"],
        ],
        dtype=float,
    )


    model = interaction_model(
        expo_ctrl,
        expo_pkai,
        hs_ctrl,
        hs_pkai,
    )


    expo_effect = float(
        np.nanmean(
            expo_pkai
        )
        - np.nanmean(
            expo_ctrl
        )
    )

    hs_effect = float(
        np.nanmean(
            hs_pkai
        )
        - np.nanmean(
            hs_ctrl
        )
    )


    tar_rows.append({
        "gene":
            r["gene"],

        "aa":
            r["aa"],

        "position":
            int(
                r["position"]
            ),

        "Targeted_ID":
            r["ID"],

        "targeted_Expo_PKAi_minus_ctrl":
            expo_effect,

        "targeted_HS_PKAi_minus_ctrl":
            hs_effect,

        "targeted_HSctrl_minus_Expoctrl":
            float(
                np.nanmean(
                    hs_ctrl
                )
                - np.nanmean(
                    expo_ctrl
                )
            ),

        "targeted_abs_PKAi_effect_attenuation_in_HS":
            (
                abs(
                    expo_effect
                )
                - abs(
                    hs_effect
                )
            ),

        "targeted_interaction_beta":
            model[
                "interaction_beta"
            ],

        "targeted_interaction_se":
            model[
                "interaction_se"
            ],

        "targeted_interaction_t":
            model[
                "interaction_t"
            ],

        "targeted_interaction_p":
            model[
                "interaction_p"
            ],
    })


tar_inter = pd.DataFrame(
    tar_rows
)


tar_inter[
    "targeted_interaction_q_all7"
] = bh(
    tar_inter[
        "targeted_interaction_p"
    ]
)


tar_inter.to_csv(
    OUT
    / "ALL7_targeted_PRM_context_interactions.tsv",
    sep="\t",
    index=False,
)


                                                               
                                             
tar_for_merge = tar_inter.merge(
    targets[
        [
            "gene",
            "aa",
            "position",
        ]
    ],
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="inner",
)


if tar_for_merge.duplicated(
    [
        "gene",
        "aa",
        "position",
    ]
).any():

    raise RuntimeError(
        "Duplicate targeted PRM exact-site key"
    )


res = res.merge(
    tar_for_merge,
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="left",
    validate="1:1",
)


res[
    "observed_targeted_PRM"
] = res[
    "Targeted_ID"
].notna()


                                                              
             
                                                              

if len(res) != 10:
    raise RuntimeError(
        f"Expected 10 final direct-A exact sites; "
        f"got {len(res)}"
    )


if res[
    [
        "gene",
        "site",
    ]
].duplicated().any():

    raise RuntimeError(
        "Final table contains duplicated exact sites"
    )


res = res.sort_values(
    [
        "observed_untargeted",
        "gene",
        "position",
    ],
    ascending=[
        False,
        True,
        True,
    ],
)


res.to_csv(
    OUT
    / "DIRECT_A_PKA_context_validation.tsv",
    sep="\t",
    index=False,
)


                                                              
            
                                                              

summary = pd.DataFrame(
    [
        {
            "direct_A_exact_sites":
                len(res),

            "untargeted_exact_overlap":
                int(
                    res[
                        "observed_untargeted"
                    ].sum()
                ),

            "untargeted_HS_effect_negative":
                int(
                    (
                        res.loc[
                            res[
                                "observed_untargeted"
                            ],
                            "Uliana_HS_effect",
                        ]
                        < 0
                    ).sum()
                ),

            "untargeted_HS_q05":
                int(
                    (
                        res[
                            "Uliana_HS_q"
                        ]
                        < 0.05
                    ).sum()
                ),

            "untargeted_Expo_q05":
                int(
                    (
                        res[
                            "Uliana_Expo_q"
                        ]
                        < 0.05
                    ).sum()
                ),

            "untargeted_interaction_p05":
                int(
                    (
                        res[
                            "interaction_p"
                        ]
                        < 0.05
                    ).sum()
                ),

            "untargeted_interaction_q05":
                int(
                    (
                        res[
                            "interaction_q_directA_overlap"
                        ]
                        < 0.05
                    ).sum()
                ),

            "targeted_PRM_exact_overlap":
                int(
                    res[
                        "observed_targeted_PRM"
                    ].sum()
                ),
        }
    ]
)


summary.to_csv(
    OUT
    / "DIRECT_A_CONTEXT_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

print()
print("=" * 110)
print("ULIANA CONTEXT-DEPENDENCE AUDIT")
print("=" * 110)

print()
print(
    summary.to_string(
        index=False
    )
)


print()
print("DIRECT-A EXACT SITES")
print("-" * 110)

cols = [
    "gene",
    "site",
    "s1_id",
    "Kanshin_heat_minus_cold_AUC",

    "observed_untargeted",

    "Uliana_Expo_effect",
    "Uliana_Expo_q",

    "Uliana_HS_effect",
    "Uliana_HS_q",

    "HSctrl_minus_Expoctrl",

    "interaction_beta",
    "interaction_p",
    "interaction_q_directA_overlap",

    "observed_targeted_PRM",

    "targeted_Expo_PKAi_minus_ctrl",
    "targeted_HS_PKAi_minus_ctrl",

    "targeted_interaction_beta",
    "targeted_interaction_p",
]


print(
    res[
        cols
    ].to_string(
        index=False
    )
)


print()
print("ALL TARGETED PRM CONTEXT INTERACTIONS")
print("-" * 110)

print(
    tar_inter.to_string(
        index=False
    )
)


print()
print("Output:", OUT)
