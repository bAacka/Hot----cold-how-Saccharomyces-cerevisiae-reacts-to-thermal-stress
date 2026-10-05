#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(
    "/"
)

MEMBERSHIP = (
    ROOT
    / "17_networkin_temporal_discovery"
    / "networkin_temporal_candidate_membership.tsv"
)

ULIANA = (
    ROOT
    / "29_external_validation"
    / "PXD052971"
    / "raw"
    / "Uliana_2026_Table_S2.xlsx"
)

OUT = (
    ROOT
    / "33_uliana_networkin_pka_validation"
)
OUT.mkdir(
    parents=True,
    exist_ok=True,
)


TARGET_LABEL = (
    "TPK2_TPK1_TPK3_group"
    "__ORTHOGONAL_NO_CURATED_AB"
)

PKA_GROUP = (
    "TPK2_TPK1_TPK3_group"
)

N_PERM = 100000
BASE_SEED = 20260915


                                                              
         
                                                              

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


def extract_sites(
    gene,
    site_id,
):
    """
    Examples:
        AVT1_T181/S187
        NUP60_S10
        YLR257W_S123
    """

    gene = clean(gene)
    sid = clean(site_id)

    prefix = gene + "_"

    if sid.startswith(prefix):
        suffix = sid[len(prefix):]
    else:
        suffix = sid

    sites = re.findall(
        r"[STY]\d+",
        suffix,
    )

    if not sites:
        raise RuntimeError(
            f"Cannot parse sites from "
            f"gene={gene!r}, site_id={sid!r}"
        )

    return sites


def empirical_lower(
    observed,
    null,
):
    return float(
        (
            1
            + np.sum(
                null <= observed
            )
        )
        / (
            len(null)
            + 1
        )
    )


def empirical_upper(
    observed,
    null,
):
    return float(
        (
            1
            + np.sum(
                null >= observed
            )
        )
        / (
            len(null)
            + 1
        )
    )


def empirical_two_sided(
    observed,
    null,
):
    center = float(
        np.mean(
            null
        )
    )

    distance = abs(
        observed
        - center
    )

    return float(
        (
            1
            + np.sum(
                np.abs(
                    null
                    - center
                )
                >= distance
            )
        )
        / (
            len(null)
            + 1
        )
    )


def atomize(
    df,
    role,
):
    rows = []

    for _, r in df.iterrows():

        gene = clean(
            r["gene"]
        )

        site_id = clean(
            r["site_id"]
        )

        sites = extract_sites(
            gene,
            site_id,
        )

        for site in sites:

            m = re.fullmatch(
                r"([STY])(\d+)",
                site,
            )

            rows.append({
                "role":
                    role,

                "gene":
                    gene,

                "site_id":
                    site_id,

                "atom_site":
                    site,

                "aa":
                    m.group(1),

                "position":
                    int(
                        m.group(2)
                    ),
            })

    return pd.DataFrame(
        rows
    )


def matched_group_table(
    atom,
    external,
):
    d = atom.merge(
        external,
        on=[
            "gene",
            "aa",
            "position",
        ],
        how="left",
        validate="m:1",
    )

    d[
        "external_match"
    ] = d[
        "Uliana_exact_key"
    ].notna()

    d.to_csv(
        OUT
        / (
            f"{d['role'].iloc[0]}"
            "_atom_level_matches.tsv"
        ),
        sep="\t",
        index=False,
    )

    rows = []

    for (
        role,
        gene,
        site_id,
    ), g in d.groupby(
        [
            "role",
            "gene",
            "site_id",
        ],
        sort=False,
    ):

        gm = g[
            g[
                "external_match"
            ]
        ].copy()

        if len(gm) == 0:
            continue

        rows.append({
            "role":
                role,

            "gene":
                gene,

            "site_id":
                site_id,

            "n_atoms_in_group":
                len(g),

            "n_exact_atoms_matched":
                len(gm),

            "matched_atoms":
                "/".join(
                    gm[
                        "atom_site"
                    ].astype(str)
                ),

            "Expo_PKAi_minus_ctrl":
                float(
                    np.nanmedian(
                        gm[
                            "Expo_PKAi_minus_ctrl"
                        ]
                    )
                ),

            "HS_PKAi_minus_ctrl":
                float(
                    np.nanmedian(
                        gm[
                            "HS_PKAi_minus_ctrl"
                        ]
                    )
                ),

            "context_interaction":
                float(
                    np.nanmedian(
                        gm[
                            "context_interaction"
                        ]
                    )
                ),

            "Expo_q_min":
                float(
                    np.nanmin(
                        gm[
                            "Expo_q"
                        ]
                    )
                )
                if gm[
                    "Expo_q"
                ].notna().any()
                else np.nan,

            "HS_q_min":
                float(
                    np.nanmin(
                        gm[
                            "HS_q"
                        ]
                    )
                )
                if gm[
                    "HS_q"
                ].notna().any()
                else np.nan,
        })

    return pd.DataFrame(
        rows
    )


def architecture_test(
    target_groups,
    background_groups,
    value_col,
    alternative,
    seed,
):
    """
    Preserve the number of matched groups per target gene.

    For each null module:
      - choose the same number of distinct background genes
      - require each chosen gene to contain at least the
        corresponding number of matched phosphogroups
      - sample that number of groups within each gene
      - collapse within gene by median
      - calculate the mean across genes
    """

    tg = target_groups[
        np.isfinite(
            target_groups[
                value_col
            ]
        )
    ].copy()

    bg = background_groups[
        np.isfinite(
            background_groups[
                value_col
            ]
        )
    ].copy()

    target_gene = (
        tg
        .groupby(
            "gene",
            as_index=False,
        )
        .agg(
            n_groups=(
                "site_id",
                "nunique",
            ),
            value=(
                value_col,
                "median",
            ),
        )
    )

    if len(target_gene) < 3:
        raise RuntimeError(
            f"Only {len(target_gene)} target genes "
            f"available for {value_col}"
        )

    multiplicities = sorted(
        target_gene[
            "n_groups"
        ].astype(int).tolist(),
        reverse=True,
    )

    observed = float(
        target_gene[
            "value"
        ].mean()
    )

    bg_by_gene = {
        gene:
            g[
                value_col
            ].to_numpy(
                dtype=float
            )

        for gene, g in bg.groupby(
            "gene"
        )
    }

    eligible = {
        m: np.array(
            [
                gene
                for gene, vals
                in bg_by_gene.items()
                if len(vals) >= m
            ],
            dtype=object,
        )

        for m in set(
            multiplicities
        )
    }

    for m in sorted(
        eligible
    ):
        if len(
            eligible[m]
        ) < len(
            multiplicities
        ):
            raise RuntimeError(
                f"Too few background genes "
                f"with >= {m} groups: "
                f"{len(eligible[m])}"
            )

    rng = np.random.default_rng(
        seed
    )

    null = np.empty(
        N_PERM,
        dtype=float,
    )

    for i in range(
        N_PERM
    ):

        chosen = set()
        scores = []

        for m in multiplicities:

            pool = eligible[m]

            while True:
                gene = pool[
                    rng.integers(
                        0,
                        len(pool),
                    )
                ]

                if gene not in chosen:
                    break

            chosen.add(
                gene
            )

            vals = bg_by_gene[
                gene
            ]

            if m == 1:
                picked = np.array([
                    vals[
                        rng.integers(
                            0,
                            len(vals),
                        )
                    ]
                ])
            else:
                idx = rng.choice(
                    len(vals),
                    size=m,
                    replace=False,
                )

                picked = vals[
                    idx
                ]

            scores.append(
                float(
                    np.median(
                        picked
                    )
                )
            )

        null[i] = float(
            np.mean(
                scores
            )
        )

    null_mean = float(
        np.mean(
            null
        )
    )

    null_sd = float(
        np.std(
            null,
            ddof=1,
        )
    )

    z = float(
        (
            observed
            - null_mean
        )
        / null_sd
    )

    if alternative == "lower":
        p = empirical_lower(
            observed,
            null,
        )

    elif alternative == "upper":
        p = empirical_upper(
            observed,
            null,
        )

    elif alternative == "two-sided":
        p = empirical_two_sided(
            observed,
            null,
        )

    else:
        raise ValueError(
            alternative
        )

    return {
        "metric":
            value_col,

        "alternative":
            alternative,

        "n_target_groups":
            len(tg),

        "n_target_genes":
            len(target_gene),

        "multiplicity_pattern":
            ",".join(
                map(
                    str,
                    multiplicities,
                )
            ),

        "observed_mean_gene_effect":
            observed,

        "observed_median_gene_effect":
            float(
                target_gene[
                    "value"
                ].median()
            ),

        "negative_genes":
            int(
                (
                    target_gene[
                        "value"
                    ] < 0
                ).sum()
            ),

        "positive_genes":
            int(
                (
                    target_gene[
                        "value"
                    ] > 0
                ).sum()
            ),

        "null_mean":
            null_mean,

        "null_sd":
            null_sd,

        "empirical_z":
            z,

        "empirical_p":
            p,
    }, target_gene


                                                              
                         
                                                              

mem = pd.read_csv(
    MEMBERSHIP,
    sep="\t",
)


target = mem[
    mem[
        "analysis_label"
    ] == TARGET_LABEL
].copy()


target = target.drop_duplicates(
    [
        "gene",
        "site_id",
    ]
)


print(
    "Orthogonal PKA rows:",
    len(target),
)


if len(target) != 20:
    print(
        target[
            [
                "gene",
                "site_id",
                "analysis_type",
                "predicted_kinase_group",
                "has_curated_PKA_AB",
            ]
        ].to_string(
            index=False
        )
    )

    raise RuntimeError(
        f"Expected 20 orthogonal PKA "
        f"measurement groups; got {len(target)}"
    )


if (
    target[
        "has_curated_PKA_AB"
    ]
    .astype(str)
    .str.lower()
    .isin(
        [
            "true",
            "1",
        ]
    )
    .any()
):
    raise RuntimeError(
        "Orthogonal target contains curated PKA A/B overlap"
    )


target_genes = set(
    target[
        "gene"
    ]
)


                              
standard = mem[
    mem[
        "analysis_type"
    ] == "standard"
].copy()


                                                          
standard_pka = standard[
    standard[
        "predicted_kinase_group"
    ] == PKA_GROUP
].copy()


standard_pka_keys = set(
    zip(
        standard_pka[
            "gene"
        ],
        standard_pka[
            "site_id"
        ],
    )
)


                         
                                                        
                                
                           
background_rows = []

for _, r in standard.iterrows():

    key = (
        r["gene"],
        r["site_id"],
    )

    if key in standard_pka_keys:
        continue

    if r["gene"] in target_genes:
        continue

    background_rows.append(
        r
    )


background = pd.DataFrame(
    background_rows
)


background = background.drop_duplicates(
    [
        "gene",
        "site_id",
    ]
)


print(
    "Competitive background groups:",
    len(background),
)


                                                              
                                       
                                                              

target_atom = atomize(
    target,
    "orthogonal_PKA",
)

background_atom = atomize(
    background,
    "networkin_non_PKA_background",
)


target_atom.to_csv(
    OUT
    / "ORTHOGONAL_PKA_atomized_candidates.tsv",
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


stat = stat[
    (stat["gene"] != "")
    & stat["aa"].isin(
        [
            "S",
            "T",
            "Y",
        ]
    )
    & pd.to_numeric(
        stat["position"],
        errors="coerce",
    ).notna()
].copy()


stat["position"] = (
    pd.to_numeric(
        stat["position"],
        errors="raise",
    )
    .astype(int)
)


                                                              
                           
stat_exact = (
    stat
    .groupby(
        [
            "gene",
            "aa",
            "position",
        ],
        as_index=False,
    )
    .agg(
        Uliana_stat_rows=(
            "ID",
            "size",
        ),

        Expo_PKAi_minus_ctrl=(
            "Expo_PKAi-Expo_ctrl",
            "median",
        ),

        HS_PKAi_minus_ctrl=(
            "HS_PKAi-HS_ctrl",
            "median",
        ),

        Expo_q=(
            "pval_adj_BH_Expo_PKAi-Expo_ctrl",
            "min",
        ),

        HS_q=(
            "pval_adj_BH_HS_PKAi-HS_ctrl",
            "min",
        ),
    )
)


                                                              
                                               
                                                              

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


ints = ints[
    (ints["gene"] != "")
    & ints["aa"].isin(
        [
            "S",
            "T",
            "Y",
        ]
    )
    & pd.to_numeric(
        ints["position"],
        errors="coerce",
    ).notna()
].copy()


ints["position"] = (
    pd.to_numeric(
        ints["position"],
        errors="raise",
    )
    .astype(int)
)


ints[
    "Expo_effect_rep"
] = (
    ints[
        [
            "Expo_PKAi_1",
            "Expo_PKAi_2",
            "Expo_PKAi_3",
        ]
    ].mean(
        axis=1
    )
    -
    ints[
        [
            "Expo_ctrl_1",
            "Expo_ctrl_2",
            "Expo_ctrl_3",
        ]
    ].mean(
        axis=1
    )
)


ints[
    "HS_effect_rep"
] = (
    ints[
        [
            "HS_PKAi_1",
            "HS_PKAi_2",
            "HS_PKAi_3",
        ]
    ].mean(
        axis=1
    )
    -
    ints[
        [
            "HS_ctrl_1",
            "HS_ctrl_2",
            "HS_ctrl_3",
        ]
    ].mean(
        axis=1
    )
)


ints[
    "context_interaction"
] = (
    ints[
        "HS_effect_rep"
    ]
    -
    ints[
        "Expo_effect_rep"
    ]
)


int_exact = (
    ints
    .groupby(
        [
            "gene",
            "aa",
            "position",
        ],
        as_index=False,
    )
    .agg(
        Uliana_intensity_rows=(
            "ID",
            "size",
        ),

        context_interaction=(
            "context_interaction",
            "median",
        ),
    )
)


                                 
external = stat_exact.merge(
    int_exact,
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="inner",
    validate="1:1",
)


external[
    "Uliana_exact_key"
] = (
    external["gene"]
    + "_"
    + external["aa"]
    + external["position"].astype(str)
)


                                                              
                                                          
                       
                                                              

target_groups = matched_group_table(
    target_atom,
    external,
)


background_groups = matched_group_table(
    background_atom,
    external,
)


target_groups.to_csv(
    OUT
    / "ORTHOGONAL_PKA_Uliana_group_matches.tsv",
    sep="\t",
    index=False,
)


background_groups.to_csv(
    OUT
    / "BACKGROUND_Uliana_group_matches.tsv",
    sep="\t",
    index=False,
)


print()
print(
    "Orthogonal PKA exact-matched measurement groups:",
    len(target_groups),
)

print(
    "Orthogonal PKA exact-matched genes:",
    target_groups[
        "gene"
    ].nunique(),
)

print(
    "Background exact-matched groups:",
    len(background_groups),
)

print(
    "Background exact-matched genes:",
    background_groups[
        "gene"
    ].nunique(),
)


if target_groups[
    "gene"
].nunique() < 3:
    raise RuntimeError(
        "Too few orthogonal PKA genes overlap Uliana "
        "for a module-level test."
    )


                                                              
                                   
 
                      
                                                   
                          
 
                    
                                 
 
                         
                                                  
                                                     
                                          
                                                              

tests = []
gene_tables = []


r, g = architecture_test(
    target_groups,
    background_groups,
    "Expo_PKAi_minus_ctrl",
    "lower",
    BASE_SEED + 1,
)

r["test_name"] = (
    "orthogonal_PKA_causal_support_exponential"
)

tests.append(
    r
)

g["test_name"] = r[
    "test_name"
]

gene_tables.append(
    g
)


r, g = architecture_test(
    target_groups,
    background_groups,
    "HS_PKAi_minus_ctrl",
    "lower",
    BASE_SEED + 2,
)

r["test_name"] = (
    "orthogonal_PKA_causal_support_heat"
)

tests.append(
    r
)

g["test_name"] = r[
    "test_name"
]

gene_tables.append(
    g
)


r, g = architecture_test(
    target_groups,
    background_groups,
    "context_interaction",
    "upper",
    BASE_SEED + 3,
)

r["test_name"] = (
    "orthogonal_PKA_positive_context_shift_in_heat"
)

tests.append(
    r
)

g["test_name"] = r[
    "test_name"
]

gene_tables.append(
    g
)


tests = pd.DataFrame(
    tests
)


tests = tests[
    [
        "test_name",
        "metric",
        "alternative",
        "n_target_groups",
        "n_target_genes",
        "multiplicity_pattern",
        "observed_mean_gene_effect",
        "observed_median_gene_effect",
        "negative_genes",
        "positive_genes",
        "null_mean",
        "null_sd",
        "empirical_z",
        "empirical_p",
    ]
]


tests.to_csv(
    OUT
    / "ORTHOGONAL_PKA_EXTERNAL_MODULE_TESTS.tsv",
    sep="\t",
    index=False,
)


pd.concat(
    gene_tables,
    ignore_index=True,
).to_csv(
    OUT
    / "ORTHOGONAL_PKA_EXTERNAL_GENE_EFFECTS.tsv",
    sep="\t",
    index=False,
)


                                                              
          
                                                              

print()
print("=" * 110)
print("ORTHOGONAL NETWORKIN PKA -> ULIANA EXTERNAL VALIDATION")
print("=" * 110)


print()
print("MODULE TESTS")
print("-" * 110)

print(
    tests.to_string(
        index=False
    )
)


print()
print("MATCHED ORTHOGONAL PKA GROUPS")
print("-" * 110)

print(
    target_groups[
        [
            "gene",
            "site_id",
            "n_atoms_in_group",
            "n_exact_atoms_matched",
            "matched_atoms",
            "Expo_PKAi_minus_ctrl",
            "Expo_q_min",
            "HS_PKAi_minus_ctrl",
            "HS_q_min",
            "context_interaction",
        ]
    ].to_string(
        index=False
    )
)


print()
print("Output:", OUT)
