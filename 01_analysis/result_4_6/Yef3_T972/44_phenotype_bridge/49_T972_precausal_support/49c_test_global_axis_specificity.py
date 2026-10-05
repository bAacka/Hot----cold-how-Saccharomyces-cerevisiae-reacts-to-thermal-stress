#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import re

import numpy as np
import pandas as pd

from scipy.stats import (
    rankdata,
    spearmanr,
    pearsonr,
)


                                                              
       
                                                              

ROOT = Path(
    "/"
)

REG = (
    ROOT /
    "19_structural_phosphosite_context" /
    "07_YEF3_T972_regulatory_evidence"
)

EXT = (
    REG /
    "external_kinase_perturbation"
)

S2 = (
    EXT /
    "NIHMS1605463-supplement-Table_S2.xlsx"
)

OUT = (
    ROOT /
    "44_phenotype_bridge" /
    "49_T972_precausal_support"
)

MODULE_FILE = (
    OUT /
    "49b_GO_MODULE_GENES.tsv"
)

TARGET = "YEF3_pT972"

DISCOVERY = {
    "ypl150w",
    "ppt1",
    "ctk1",
}

LEAKAGE_GENES = {
    "YEF3",
    "YPL150W",
    "PPT1",
    "CTK1",
}

MODULES = [
    "PROTEIN_FOLDING",
    "TRANSLATION_ELONGATION",
    "STRESS_GRANULE_RNP",
]

MIN_MEASURED_GENES = 5

N_PERM = 20_000
N_RANDOM = 5_000

SEED = 20260922

rng = np.random.default_rng(
    SEED
)


                                                              
         
                                                              

def clean_site(x):

    return re.sub(
        r"\s+",
        "",
        str(x),
    )


def parse_del_col(c):

    m = re.fullmatch(
        r"[∆Δ]([^_]+)_Rep(\d+)",
        str(c).strip(),
        flags=re.I,
    )

    if not m:
        return None

    return (
        m.group(1).lower(),
        int(m.group(2)),
    )


def parse_site_gene(site):

    m = re.match(
        r"^(.+?)_p",
        clean_site(site),
        flags=re.I,
    )

    if not m:
        return ""

    return (
        m.group(1)
        .upper()
    )


def bh_fdr(pvals):

    p = np.asarray(
        pvals,
        dtype=float,
    )

    result = np.full(
        len(p),
        np.nan,
    )

    good = np.isfinite(
        p
    )

    if not good.any():
        return result

    idx = np.where(
        good
    )[0]

    x = p[
        idx
    ]

    order = np.argsort(
        x
    )

    sx = x[
        order
    ]

    m = len(
        sx
    )

    q = (
        sx
        *
        m
        /
        np.arange(
            1,
            m + 1,
        )
    )

    q = np.minimum.accumulate(
        q[::-1]
    )[::-1]

    q = np.minimum(
        q,
        1.0,
    )

    tmp = np.empty(
        m
    )

    tmp[
        order
    ] = q

    result[
        idx
    ] = tmp

    return result


def rank_z(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    r = rankdata(
        x,
        method="average",
    )

    sd = r.std(
        ddof=0,
    )

    if sd == 0:
        raise RuntimeError(
            "Zero-variance rank vector"
        )

    return (
        r - r.mean()
    ) / sd


def residualize_on_global(
    x,
    global_score,
):

    x = rank_z(
        x
    )

    g = rank_z(
        global_score
    )

    X = np.column_stack([
        np.ones(
            len(g)
        ),
        g,
    ])

    beta = np.linalg.lstsq(
        X,
        x,
        rcond=None,
    )[0]

    return (
        x
        -
        X @ beta
    )


def residual_corr(
    a,
    b,
    g,
):

    ra = residualize_on_global(
        a,
        g,
    )

    rb = residualize_on_global(
        b,
        g,
    )

    return float(
        pearsonr(
            ra,
            rb,
        ).statistic
    )


def permutation_p(
    target_residual,
    programme_residual,
):

    obs = float(
        pearsonr(
            target_residual,
            programme_residual,
        ).statistic
    )

    exceed = 0

    for _ in range(
        N_PERM
    ):

        rp = rng.permutation(
            target_residual
        )

        stat = float(
            pearsonr(
                rp,
                programme_residual,
            ).statistic
        )

        if (
            abs(stat)
            >=
            abs(obs) - 1e-15
        ):

            exceed += 1

    return (
        obs,
        (
            1 + exceed
        )
        /
        (
            N_PERM + 1
        ),
    )


                                                              
                          
                                                              

print(
    "=" * 110
)

print(
    "LOAD FROZEN 49b MODULE MEMBERSHIP"
)

print(
    "=" * 110
)

m = pd.read_csv(
    MODULE_FILE,
    sep="\t",
)

required = {
    "module",
    "gene",
    "present_in_phosphoproteome",
}

if not required.issubset(
    m.columns
):

    raise RuntimeError(
        "Unexpected 49b module file schema"
    )


present = (
    m[
        "present_in_phosphoproteome"
    ]
    .astype(str)
    .str.lower()
    .isin(
        {
            "true",
            "1",
            "yes",
        }
    )
)

module_genes = {}

for module in MODULES:

    genes = sorted(
        set(
            m.loc[
                m[
                    "module"
                ].eq(
                    module
                )
                &
                present,
                "gene",
            ]
            .astype(str)
            .str.upper()
        )
    )

    module_genes[
        module
    ] = genes

    print(
        module,
        len(
            genes
        ),
    )


                                                              
         
                                                              

print()
print(
    "=" * 110
)
print(
    "LOAD PROTEIN-NORMALIZED DELETION PHOSPHOPROTEOME"
)
print(
    "=" * 110
)

x = pd.read_excel(
    S2,
    sheet_name="Phosphosite wProNormalization",
)

x = x.dropna(
    axis=1,
    how="all",
).copy()


SITE_COL = (
    "Gene Symbol pSTY position"
)

if SITE_COL not in x.columns:

    raise RuntimeError(
        f"Missing {SITE_COL}"
    )


x["_site"] = (
    x[
        SITE_COL
    ]
    .map(
        clean_site
    )
)

x["_gene"] = (
    x["_site"]
    .map(
        parse_site_gene
    )
)


                                                              
                  
                                                              

deletion_cols = defaultdict(
    list
)

for c in x.columns:

    z = parse_del_col(
        c
    )

    if z is None:
        continue

    gene, rep = z

    deletion_cols[
        gene
    ].append(
        c
    )


print(
    "deletion strains:",
    len(
        deletion_cols
    ),
)

if len(
    deletion_cols
) != 110:

    raise RuntimeError(
        "Expected 110 deletion strains; "
        f"found {len(deletion_cols)}"
    )


                                                              
                        
                                                              

site_cols = {}

for deletion, cols in sorted(
    deletion_cols.items()
):

    q = (
        x[
            cols
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1,
            skipna=True,
        )
        .to_numpy(
            dtype=float
        )
    )

    site_cols[
        deletion
    ] = q


site_matrix = pd.DataFrame(
    site_cols
)

site_matrix.index = pd.Index(
    x[
        "_site"
    ].to_numpy(),
    name="_site",
)


finite = int(
    np.isfinite(
        site_matrix
        .to_numpy(
            float
        )
    ).sum()
)

print(
    "finite site x deletion:",
    f"{finite:,}",
)

if finite < 100_000:

    raise RuntimeError(
        "Unexpectedly sparse matrix"
    )


                                                              
      
                                                              

target_idx = np.where(
    x[
        "_site"
    ]
    .eq(
        TARGET
    )
    .to_numpy()
)[0]

if len(
    target_idx
) != 1:

    raise RuntimeError(
        f"Expected one {TARGET}; "
        f"found {len(target_idx)}"
    )


t972 = (
    site_matrix
    .iloc[
        int(
            target_idx[0]
        )
    ]
)

print(
    "T972 measured deletions:",
    int(
        t972.notna().sum()
    ),
)

if int(
    t972.notna().sum()
) != 110:

    raise RuntimeError(
        "Expected T972 in all 110 deletions"
    )


                                                              
                               
                                                              

tmp = (
    site_matrix
    .copy()
)

tmp["_gene"] = (
    x[
        "_gene"
    ]
    .to_numpy()
)

tmp = tmp[
    tmp[
        "_gene"
    ].ne("")
].copy()


gene_matrix = (
    tmp
    .groupby(
        "_gene",
        sort=True,
    )
    .median(
        numeric_only=True
    )
)


print(
    "gene-level matrix:",
    gene_matrix.shape,
)


                                                              
                               
                                                              

global_genes = sorted(
    set(
        gene_matrix.index
    )
    -
    LEAKAGE_GENES
)

global_matrix = (
    gene_matrix
    .loc[
        global_genes
    ]
)

global_score = (
    global_matrix
    .median(
        axis=0,
        skipna=True,
    )
)

global_n = (
    global_matrix
    .notna()
    .sum(
        axis=0
    )
)


print(
    "global genes:",
    len(
        global_genes
    ),
)

print(
    "global measured genes / deletion:",
    int(
        global_n.min()
    ),
    "to",
    int(
        global_n.max()
    ),
)


                                                              
                  
                                                              

programme_scores = {}

programme_counts = {}

for module, genes in module_genes.items():

    gm = gene_matrix.loc[
        genes
    ]

    score = gm.median(
        axis=0,
        skipna=True,
    )

    count = gm.notna().sum(
        axis=0
    )

    score[
        count
        <
        MIN_MEASURED_GENES
    ] = np.nan

    programme_scores[
        module
    ] = score

    programme_counts[
        module
    ] = count


                                                              
                                     
                                                              

primary_deletions = [
    d
    for d in sorted(
        deletion_cols
    )
    if d
    not in
    DISCOVERY
]


if len(
    primary_deletions
) != 107:

    raise RuntimeError(
        f"Expected 107 primary deletions; "
        f"found {len(primary_deletions)}"
    )


                                                              
                         
                                                              

diag_rows = []

yy = np.array(
    [
        t972[
            d
        ]
        for d
        in
        primary_deletions
    ],
    dtype=float,
)

gg = np.array(
    [
        global_score[
            d
        ]
        for d
        in
        primary_deletions
    ],
    dtype=float,
)


sg = spearmanr(
    yy,
    gg,
)

diag_rows.append({
    "variable":
        "YEF3_T972",

    "spearman_vs_global":
        float(
            sg.statistic
        ),

    "p":
        float(
            sg.pvalue
        ),
})


for module in MODULES:

    xx = np.array(
        [
            programme_scores[
                module
            ][
                d
            ]
            for d
            in
            primary_deletions
        ],
        dtype=float,
    )

    ok = (
        np.isfinite(
            xx
        )
        &
        np.isfinite(
            gg
        )
    )

    sr = spearmanr(
        xx[
            ok
        ],
        gg[
            ok
        ],
    )

    diag_rows.append({
        "variable":
            module,

        "spearman_vs_global":
            float(
                sr.statistic
            ),

        "p":
            float(
                sr.pvalue
            ),
    })


diagnostics = pd.DataFrame(
    diag_rows
)


                                                              
                                
                                                              

result_rows = []

for module in MODULES:

    xx = np.array(
        [
            programme_scores[
                module
            ][
                d
            ]
            for d
            in
            primary_deletions
        ],
        dtype=float,
    )

    yy = np.array(
        [
            t972[
                d
            ]
            for d
            in
            primary_deletions
        ],
        dtype=float,
    )

    gg = np.array(
        [
            global_score[
                d
            ]
            for d
            in
            primary_deletions
        ],
        dtype=float,
    )

    ok = (
        np.isfinite(
            xx
        )
        &
        np.isfinite(
            yy
        )
        &
        np.isfinite(
            gg
        )
    )

    xx = xx[
        ok
    ]

    yy = yy[
        ok
    ]

    gg = gg[
        ok
    ]

    raw = spearmanr(
        yy,
        xx,
    )

    target_resid = (
        residualize_on_global(
            yy,
            gg,
        )
    )

    module_resid = (
        residualize_on_global(
            xx,
            gg,
        )
    )

    partial_rho, perm_p = (
        permutation_p(
            target_resid,
            module_resid,
        )
    )

    result_rows.append({
        "module":
            module,

        "n_deletions":
            len(
                xx
            ),

        "raw_spearman_rho":
            float(
                raw.statistic
            ),

        "raw_spearman_p":
            float(
                raw.pvalue
            ),

        "partial_spearman_global_adjusted":
            float(
                partial_rho
            ),

        "permutation_p_two_sided":
            float(
                perm_p
            ),
    })


results = pd.DataFrame(
    result_rows
)

results[
    "BH_q_across_3_modules"
] = bh_fdr(
    results[
        "permutation_p_two_sided"
    ].to_numpy(
        float
    )
)


                                                              
                                                     
                                                              

print()
print(
    "=" * 110
)
print(
    "GLOBAL-ADJUSTED RANDOM GENE-SET NULL"
)
print(
    "=" * 110
)


eligible_genes = sorted(
    set(
        gene_matrix.index
    )
    -
    LEAKAGE_GENES
)


                                         

                                                     
target_vec = np.array(
    [
        t972[d]
        for d in primary_deletions
    ],
    dtype=float,
)

global_vec = np.array(
    [
        global_score[d]
        for d in primary_deletions
    ],
    dtype=float,
)

target_resid_fixed = (
    residualize_on_global(
        target_vec,
        global_vec,
    )
)


random_rows = []

for module in MODULES:

    size = len(
        module_genes[
            module
        ]
    )

    observed = float(
        results.loc[
            results[
                "module"
            ].eq(
                module
            ),
            "partial_spearman_global_adjusted",
        ].iloc[0]
    )

    null = []

    for _ in range(
        N_RANDOM
    ):

        genes = rng.choice(
            eligible_genes,
            size=size,
            replace=False,
        )

        gm = gene_matrix.loc[
            genes,
            primary_deletions,
        ]

        score = gm.median(
            axis=0,
            skipna=True,
        )

        count = gm.notna().sum(
            axis=0
        )

        score[
            count
            <
            MIN_MEASURED_GENES
        ] = np.nan

        random_vec = score.to_numpy(
            dtype=float
        )

        ok = (
            np.isfinite(
                random_vec
            )
            &
            np.isfinite(
                target_vec
            )
            &
            np.isfinite(
                global_vec
            )
        )

        if ok.sum() < 10:
            continue

        tr = residualize_on_global(
            target_vec[
                ok
            ],
            global_vec[
                ok
            ],
        )

        rr = residualize_on_global(
            random_vec[
                ok
            ],
            global_vec[
                ok
            ],
        )

        stat = float(
            pearsonr(
                tr,
                rr,
            ).statistic
        )

        if np.isfinite(
            stat
        ):
            null.append(
                stat
            )


    null = np.asarray(
        null,
        dtype=float,
    )

    empirical = (
        1
        +
        int(
            np.sum(
                np.abs(
                    null
                )
                >=
                abs(
                    observed
                )
            )
        )
    ) / (
        len(
            null
        )
        + 1
    )


    random_rows.append({
        "module":
            module,

        "module_size":
            size,

        "n_random_sets":
            len(
                null
            ),

        "observed_partial_rho":
            observed,

        "random_set_empirical_p_two_sided":
            empirical,

        "null_abs_rho_median":
            float(
                np.median(
                    np.abs(
                        null
                    )
                )
            ),

        "null_abs_rho_q95":
            float(
                np.quantile(
                    np.abs(
                        null
                    ),
                    0.95,
                )
            ),
    })


random_summary = pd.DataFrame(
    random_rows
)


                                                              
                                                           
                  
                                                              

disc_rows = []

for deletion in [
    "ypl150w",
    "ppt1",
    "ctk1",
]:

    row = {
        "deleted_gene":
            deletion,

        "YEF3_T972":
            t972[
                deletion
            ],

        "global_score":
            global_score[
                deletion
            ],
    }

    for module in MODULES:

        row[
            f"{module}_score"
        ] = (
            programme_scores[
                module
            ][
                deletion
            ]
        )

    disc_rows.append(
        row
    )


disc = pd.DataFrame(
    disc_rows
)


                                                              
       
                                                              

diagnostics.to_csv(
    OUT /
    "49c_GLOBAL_AXIS_DIAGNOSTICS.tsv",
    sep="\t",
    index=False,
)

results.to_csv(
    OUT /
    "49c_PRIMARY_GLOBAL_ADJUSTED_RESULTS.tsv",
    sep="\t",
    index=False,
)

random_summary.to_csv(
    OUT /
    "49c_GLOBAL_ADJUSTED_RANDOM_GENESET_NULL.tsv",
    sep="\t",
    index=False,
)

disc.to_csv(
    OUT /
    "49c_DISCOVERY_PERTURBATIONS_GLOBAL_CONTEXT.tsv",
    sep="\t",
    index=False,
)


print()
print(
    "=" * 110
)
print(
    "GLOBAL AXIS DIAGNOSTICS"
)
print(
    "=" * 110
)

print(
    diagnostics.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6g}",
    )
)


print()
print(
    "=" * 110
)
print(
    "PRIMARY GLOBAL-ADJUSTED RESULTS"
)
print(
    "=" * 110
)

print(
    results.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6g}",
    )
)


print()
print(
    "=" * 110
)
print(
    "GLOBAL-ADJUSTED RANDOM-GENE-SET SPECIFICITY"
)
print(
    "=" * 110
)

print(
    random_summary.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6g}",
    )
)


print()
print(
    "=" * 110
)
print(
    "DISCOVERY PERTURBATIONS — DESCRIPTIVE ONLY"
)
print(
    "=" * 110
)

print(
    disc.to_string(
        index=False,
        float_format=lambda z:
            f"{z:+.6g}",
    )
)


print()
print(
    "=" * 110
)
print(
    "WRITTEN"
)
print(
    "=" * 110
)

for f in [
    "49c_GLOBAL_AXIS_DIAGNOSTICS.tsv",
    "49c_PRIMARY_GLOBAL_ADJUSTED_RESULTS.tsv",
    "49c_GLOBAL_ADJUSTED_RANDOM_GENESET_NULL.tsv",
    "49c_DISCOVERY_PERTURBATIONS_GLOBAL_CONTEXT.tsv",
]:

    print(
        OUT / f
    )
