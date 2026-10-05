#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import re

import numpy as np
import pandas as pd

from scipy.stats import (
    rankdata,
    pearsonr,
    spearmanr,
)


                                                              
       
                                                              

ROOT = Path(
    "/"
)

STR = (
    ROOT /
    "19_structural_phosphosite_context"
)

REG = (
    STR /
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

GAF = (
    ROOT /
    "18_functional_programs" /
    "resources" /
    "sgd.gaf"
)

OBO = (
    ROOT /
    "18_functional_programs" /
    "resources" /
    "go-basic.obo"
)

OUT = (
    ROOT /
    "44_phenotype_bridge" /
    "49_T972_precausal_support"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                  
                                                              

TARGET = "YEF3_pT972"

DISCOVERY_DELETIONS = {
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

MODULE_ROOTS = {
    "PROTEIN_FOLDING": [
        "GO:0006457",
    ],

    "TRANSLATION_ELONGATION": [
        "GO:0006414",
    ],

    "STRESS_GRANULE_RNP": [
        "GO:0034063",
        "GO:0010494",
    ],
}

MIN_MEASURED_GENES = 5

N_PERM = 20_000
N_RANDOM_GENESETS = 5_000

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

    s = clean_site(site)

    m = re.match(
        r"^(.+?)_p",
        s,
        flags=re.I,
    )

    if not m:
        return ""

    return m.group(1).upper()


def bh_fdr(pvals):

    p = np.asarray(
        pvals,
        dtype=float,
    )

    out = np.full(
        len(p),
        np.nan,
    )

    ok = np.isfinite(p)

    if not ok.any():
        return out

    idx = np.where(ok)[0]
    v = p[idx]

    order = np.argsort(v)
    ranked = v[order]

    m = len(ranked)

    q = (
        ranked
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
        m,
        dtype=float,
    )

    tmp[order] = q

    out[idx] = tmp

    return out


def standardized_rank(x):

    r = rankdata(
        np.asarray(
            x,
            dtype=float,
        ),
        method="average",
    )

    sd = r.std(
        ddof=0,
    )

    if sd == 0:
        return None

    return (
        r - r.mean()
    ) / sd


def permutation_spearman_p(
    x,
    y,
    n_perm,
):

    x = np.asarray(
        x,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )

    rx = standardized_rank(x)
    ry = standardized_rank(y)

    if (
        rx is None
        or
        ry is None
    ):
        return np.nan

    obs = float(
        np.mean(
            rx * ry
        )
    )

    exceed = 0

    for _ in range(
        n_perm
    ):

        yp = rng.permutation(
            ry
        )

        stat = float(
            np.mean(
                rx * yp
            )
        )

        if (
            abs(stat)
            >=
            abs(obs) - 1e-15
        ):
            exceed += 1

    return (
        1.0 + exceed
    ) / (
        n_perm + 1.0
    )


                                                              
                     
                                                              

print(
    "=" * 110
)

print(
    "LOAD GO HIERARCHY"
)

print(
    "=" * 110
)

if not OBO.is_file():
    raise RuntimeError(
        f"Missing OBO: {OBO}"
    )

parents = defaultdict(set)

current_id = None
current_obsolete = False

with OBO.open(
    "r",
    encoding="utf-8",
    errors="replace",
) as fh:

    for raw in fh:

        line = raw.rstrip(
            "\n"
        )

        if line == "[Term]":

            current_id = None
            current_obsolete = False
            continue

        if line.startswith(
            "["
        ):

            current_id = None
            current_obsolete = False
            continue

        if line.startswith(
            "id: GO:"
        ):

            current_id = (
                line.split(
                    "id:",
                    1,
                )[1]
                .strip()
            )

            continue

        if (
            current_id
            and
            line.startswith(
                "is_obsolete: true"
            )
        ):

            current_obsolete = True
            continue

        if (
            current_id
            and
            not current_obsolete
            and
            line.startswith(
                "is_a: GO:"
            )
        ):

            parent = (
                line.split(
                    "is_a:",
                    1,
                )[1]
                .split()[0]
                .strip()
            )

            parents[
                current_id
            ].add(
                parent
            )


children = defaultdict(set)

for child, ps in parents.items():

    for parent in ps:

        children[
            parent
        ].add(
            child
        )


def descendants_including_self(
    roots,
):

    seen = set(
        roots
    )

    stack = list(
        roots
    )

    while stack:

        x = stack.pop()

        for c in children.get(
            x,
            set(),
        ):

            if c not in seen:

                seen.add(
                    c
                )

                stack.append(
                    c
                )

    return seen


module_terms = {}

for module, roots in MODULE_ROOTS.items():

    terms = descendants_including_self(
        roots
    )

    module_terms[
        module
    ] = terms

    print(
        module,
        "GO terms:",
        len(terms),
    )


                                                              
               
                                                              

print()
print(
    "=" * 110
)
print(
    "LOAD SGD GAF"
)
print(
    "=" * 110
)

if not GAF.is_file():
    raise RuntimeError(
        f"Missing GAF: {GAF}"
    )

go_to_genes = defaultdict(set)

with GAF.open(
    "r",
    encoding="utf-8",
    errors="replace",
) as fh:

    for line in fh:

        if (
            not line
            or
            line.startswith("!")
        ):
            continue

        f = line.rstrip(
            "\n"
        ).split(
            "\t"
        )

        if len(f) < 5:
            continue

        symbol = (
            f[2]
            .strip()
            .upper()
        )

        qualifier = (
            f[3]
            .strip()
        )

        go_id = (
            f[4]
            .strip()
        )

        if (
            not symbol
            or
            not go_id
        ):
            continue

        qualifiers = {
            q.strip().upper()
            for q in qualifier.split("|")
            if q.strip()
        }

        if "NOT" in qualifiers:
            continue

        go_to_genes[
            go_id
        ].add(
            symbol
        )


module_genes_raw = {}

for module, terms in module_terms.items():

    genes = set()

    for term in terms:

        genes.update(
            go_to_genes.get(
                term,
                set(),
            )
        )

    genes -= LEAKAGE_GENES

    module_genes_raw[
        module
    ] = genes

    print(
        module,
        "GO genes after leakage exclusions:",
        len(genes),
    )


                                                              
                                                   
                                                              

print()
print(
    "=" * 110
)
print(
    "LOAD TABLE S2 — PROTEIN-NORMALIZED PHOSPHOSITES"
)
print(
    "=" * 110
)

if not S2.is_file():
    raise RuntimeError(
        f"Missing S2: {S2}"
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
        f"Missing site column: {SITE_COL}"
    )


x["_site"] = (
    x[SITE_COL]
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


                                                              
                                     
                                                              

deletion_cols = defaultdict(list)

for c in x.columns:

    parsed = parse_del_col(
        c
    )

    if parsed is None:
        continue

    gene, rep = parsed

    deletion_cols[
        gene
    ].append(
        c
    )


print(
    "deletion strains:",
    len(deletion_cols),
)

if len(
    deletion_cols
) < 100:

    raise RuntimeError(
        "Expected approximately 110 deletion strains; "
        f"found {len(deletion_cols)}"
    )


                                                              
                                                      
                                                              

site_matrix_cols = {}

for deletion, cols in sorted(
    deletion_cols.items()
):

    z = x[
        cols
    ].apply(
        pd.to_numeric,
        errors="coerce",
    )

                
                                              
     
                                                                 
                                                                   
                                                                  
                     
     
                                                                
    site_matrix_cols[
        deletion
    ] = (
        z.mean(
            axis=1,
            skipna=True,
        )
        .to_numpy(
            dtype=float
        )
    )


site_matrix = pd.DataFrame(
    site_matrix_cols
)

site_matrix.index = pd.Index(
    x["_site"].to_numpy(),
    name="_site",
)


                                                              
                                                      
                                                              

n_finite = int(
    np.isfinite(
        site_matrix.to_numpy(
            dtype=float
        )
    ).sum()
)

print(
    "finite site x deletion measurements:",
    f"{n_finite:,}",
)

if n_finite < 100_000:
    raise RuntimeError(
        "Phosphosite-by-deletion matrix contains unexpectedly "
        f"few finite measurements: {n_finite:,}"
    )


target_mask_qc = (
    x["_site"]
    .eq(
        TARGET
    )
    .to_numpy()
)

if int(
    target_mask_qc.sum()
) != 1:
    raise RuntimeError(
        f"Expected exactly one {TARGET} row before target extraction; "
        f"found {int(target_mask_qc.sum())}"
    )

target_qc = (
    site_matrix
    .iloc[
        np.where(
            target_mask_qc
        )[0][0]
    ]
)

n_target_qc = int(
    target_qc.notna().sum()
)

print(
    "T972 finite measurements after matrix construction:",
    n_target_qc,
)

if n_target_qc == 0:
    raise RuntimeError(
        "T972 is still entirely missing after matrix construction."
    )


                                                              
                    
                                                              

target_rows = np.where(
    x["_site"].eq(
        TARGET
    ).to_numpy()
)[0]

if len(
    target_rows
) != 1:

    raise RuntimeError(
        f"Expected exactly one {TARGET} row; "
        f"found {len(target_rows)}"
    )

target_i = int(
    target_rows[0]
)

t972 = (
    site_matrix
    .iloc[
        target_i
    ]
    .copy()
)

t972.name = (
    "YEF3_T972_log2ratio"
)

print(
    "T972 measured deletions:",
    int(
        t972.notna().sum()
    ),
)


                                                              
                                           
 
                                                      
                                                              

tmp = (
    site_matrix
    .copy()
)

tmp["_gene"] = (
    x["_gene"]
    .to_numpy()
)

tmp = tmp[
    tmp["_gene"].ne("")
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
    "genes with phosphosite measurements:",
    len(
        gene_matrix
    ),
)


                                                              
                                          
                                                              

module_genes = {}

module_gene_rows = []

for module, genes in module_genes_raw.items():

    observed = sorted(
        genes.intersection(
            set(
                gene_matrix.index
            )
        )
    )

    module_genes[
        module
    ] = observed

    for gene in sorted(
        genes
    ):

        module_gene_rows.append({
            "module":
                module,

            "gene":
                gene,

            "present_in_phosphoproteome":
                gene in gene_matrix.index,
        })

    print()
    print(
        module,
        "assay-covered genes:",
        len(observed),
    )

    print(
        ", ".join(
            observed[:40]
        )
        +
        (
            " ..."
            if len(observed) > 40
            else ""
        )
    )


pd.DataFrame(
    module_gene_rows
).to_csv(
    OUT /
    "49b_GO_MODULE_GENES.tsv",
    sep="\t",
    index=False,
)


                                                              
                           
                                                              

score_rows = []

module_scores = {}
module_counts = {}

for module, genes in module_genes.items():

    if len(
        genes
    ) < MIN_MEASURED_GENES:

        print(
            f"WARNING: {module} has only "
            f"{len(genes)} assay-covered genes"
        )

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

    module_scores[
        module
    ] = score

    module_counts[
        module
    ] = count


for deletion in sorted(
    deletion_cols
):

    row = {
        "deleted_gene":
            deletion,

        "YEF3_T972_log2ratio":
            t972.get(
                deletion,
                np.nan,
            ),

        "is_discovery_deletion":
            deletion
            in
            DISCOVERY_DELETIONS,
    }

    for module in MODULE_ROOTS:

        row[
            f"{module}_score"
        ] = (
            module_scores[
                module
            ].get(
                deletion,
                np.nan,
            )
        )

        row[
            f"{module}_n_genes"
        ] = int(
            module_counts[
                module
            ].get(
                deletion,
                0,
            )
        )

    score_rows.append(
        row
    )


scores = pd.DataFrame(
    score_rows
)

scores.to_csv(
    OUT /
    "49b_T972_MODULE_SCORES_BY_DELETION.tsv",
    sep="\t",
    index=False,
)


                                                              
                    
                                                              

def evaluate(
    module,
    deletion_filter,
    label,
):

    score_col = (
        f"{module}_score"
    )

    z = scores[
        deletion_filter(
            scores,
            module,
        )
    ][
        [
            "deleted_gene",
            "YEF3_T972_log2ratio",
            score_col,
        ]
    ].dropna()

    if len(z) < 10:

        return {
            "analysis":
                label,

            "module":
                module,

            "n_deletions":
                len(z),

            "spearman_rho":
                np.nan,

            "pearson_r":
                np.nan,

            "permutation_p_two_sided":
                np.nan,
        }

    a = (
        z[
            "YEF3_T972_log2ratio"
        ]
        .to_numpy(
            float
        )
    )

    b = (
        z[
            score_col
        ]
        .to_numpy(
            float
        )
    )

    sr = spearmanr(
        a,
        b,
    )

    pr = pearsonr(
        a,
        b,
    )

    pp = permutation_spearman_p(
        a,
        b,
        N_PERM,
    )

    return {
        "analysis":
            label,

        "module":
            module,

        "n_deletions":
            len(z),

        "spearman_rho":
            float(
                sr.statistic
            ),

        "spearman_asymptotic_p":
            float(
                sr.pvalue
            ),

        "pearson_r":
            float(
                pr.statistic
            ),

        "pearson_p":
            float(
                pr.pvalue
            ),

        "permutation_p_two_sided":
            float(
                pp
            ),
    }


def filt_primary(
    df,
    module,
):

    return (
        ~df[
            "is_discovery_deletion"
        ]
        &
        df[
            "YEF3_T972_log2ratio"
        ].notna()
        &
        df[
            f"{module}_score"
        ].notna()
    )


def filt_all(
    df,
    module,
):

    return (
        df[
            "YEF3_T972_log2ratio"
        ].notna()
        &
        df[
            f"{module}_score"
        ].notna()
    )


def filt_primary_no_self(
    df,
    module,
):

    members = {
        x.lower()
        for x in module_genes[
            module
        ]
    }

    return (
        filt_primary(
            df,
            module,
        )
        &
        ~df[
            "deleted_gene"
        ].isin(
            members
        )
    )


                                                              
              
                                                              

primary_rows = []

for module in MODULE_ROOTS:

    primary_rows.append(
        evaluate(
            module,
            filt_primary,
            "PRIMARY_EXCLUDE_YPL150W_PPT1_CTK1",
        )
    )


primary = pd.DataFrame(
    primary_rows
)

primary[
    "BH_q_across_3_modules"
] = bh_fdr(
    primary[
        "permutation_p_two_sided"
    ].to_numpy(
        float
    )
)


                                                              
                   
                                                              

sensitivity_rows = []

for module in MODULE_ROOTS:

    sensitivity_rows.append(
        evaluate(
            module,
            filt_all,
            "SENSITIVITY_INCLUDE_ALL_DELETIONS",
        )
    )

    sensitivity_rows.append(
        evaluate(
            module,
            filt_primary_no_self,
            "SENSITIVITY_PRIMARY_PLUS_EXCLUDE_SELF_DELETIONS",
        )
    )


sensitivity = pd.DataFrame(
    sensitivity_rows
)


                                                              
                                   
 
                             
                                                              

print()
print(
    "=" * 110
)
print(
    "MATCHED-SIZE RANDOM GENE-SET NULL"
)
print(
    "=" * 110
)

eligible_pool = sorted(
    set(
        gene_matrix.index
    )
    -
    LEAKAGE_GENES
)

random_rows = []

primary_deletions = [
    d
    for d in scores[
        "deleted_gene"
    ].tolist()
    if (
        d
        not in
        DISCOVERY_DELETIONS
        and
        np.isfinite(
            t972.get(
                d,
                np.nan,
            )
        )
    )
]


for module, genes in module_genes.items():

    m = len(
        genes
    )

    if (
        m < MIN_MEASURED_GENES
        or
        m > len(
            eligible_pool
        )
    ):

        random_rows.append({
            "module":
                module,

            "module_size":
                m,

            "n_valid_random_sets":
                0,

            "observed_primary_spearman_rho":
                primary.loc[
                    primary[
                        "module"
                    ].eq(
                        module
                    ),
                    "spearman_rho",
                ].iloc[0],

            "random_gene_set_empirical_p_two_sided":
                np.nan,

            "random_abs_rho_median":
                np.nan,

            "random_abs_rho_q95":
                np.nan,
        })

        continue

    observed = float(
        primary.loc[
            primary[
                "module"
            ].eq(
                module
            ),
            "spearman_rho",
        ].iloc[0]
    )

    null = []

    for _ in range(
        N_RANDOM_GENESETS
    ):

        chosen = rng.choice(
            eligible_pool,
            size=m,
            replace=False,
        )

        gm = gene_matrix.loc[
            chosen,
            primary_deletions,
        ]

        rs = gm.median(
            axis=0,
            skipna=True,
        )

        rc = gm.notna().sum(
            axis=0
        )

        rs[
            rc
            <
            MIN_MEASURED_GENES
        ] = np.nan

        yy = np.array(
            [
                t972.get(
                    d,
                    np.nan,
                )
                for d
                in primary_deletions
            ],
            dtype=float,
        )

        xx = rs.to_numpy(
            float
        )

        ok = (
            np.isfinite(
                xx
            )
            &
            np.isfinite(
                yy
            )
        )

        if ok.sum() < 10:
            continue

        rho = spearmanr(
            yy[
                ok
            ],
            xx[
                ok
            ],
        ).statistic

        if np.isfinite(
            rho
        ):
            null.append(
                float(
                    rho
                )
            )


    null = np.asarray(
        null,
        dtype=float,
    )

    if len(
        null
    ) == 0:

        emp = np.nan
        med = np.nan
        q95 = np.nan

    else:

        emp = (
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

        med = float(
            np.median(
                np.abs(
                    null
                )
            )
        )

        q95 = float(
            np.quantile(
                np.abs(
                    null
                ),
                0.95,
            )
        )

    random_rows.append({
        "module":
            module,

        "module_size":
            m,

        "n_valid_random_sets":
            len(
                null
            ),

        "observed_primary_spearman_rho":
            observed,

        "random_gene_set_empirical_p_two_sided":
            emp,

        "random_abs_rho_median":
            med,

        "random_abs_rho_q95":
            q95,
    })


random_summary = pd.DataFrame(
    random_rows
)


                                                              
                                          
                                                              

disc = scores[
    scores[
        "deleted_gene"
    ].isin(
        DISCOVERY_DELETIONS
    )
].copy()

disc.to_csv(
    OUT /
    "49b_DISCOVERY_REGULATOR_PROJECTION.tsv",
    sep="\t",
    index=False,
)


                                                              
               
                                                              

primary.to_csv(
    OUT /
    "49b_PRIMARY_T972_PROGRAM_CORRELATIONS.tsv",
    sep="\t",
    index=False,
)

sensitivity.to_csv(
    OUT /
    "49b_SENSITIVITY_T972_PROGRAM_CORRELATIONS.tsv",
    sep="\t",
    index=False,
)

random_summary.to_csv(
    OUT /
    "49b_RANDOM_GENESET_NULL_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

print()
print(
    "=" * 110
)
print(
    "PRIMARY: T972 vs PREDEFINED PROGRAMMES"
)
print(
    "=" * 110
)

print(
    primary.to_string(
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
    "RANDOM GENE-SET SPECIFICITY"
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
    "DISCOVERY PERTURBATIONS — NOT INCLUDED IN PRIMARY CORRELATION"
)
print(
    "=" * 110
)

show = [
    "deleted_gene",
    "YEF3_T972_log2ratio",
]

for module in MODULE_ROOTS:
    show += [
        f"{module}_score",
        f"{module}_n_genes",
    ]

print(
    disc[
        show
    ].to_string(
        index=False,
        float_format=lambda z:
            f"{z:+.5f}",
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

for p in [
    "49b_GO_MODULE_GENES.tsv",
    "49b_T972_MODULE_SCORES_BY_DELETION.tsv",
    "49b_PRIMARY_T972_PROGRAM_CORRELATIONS.tsv",
    "49b_SENSITIVITY_T972_PROGRAM_CORRELATIONS.tsv",
    "49b_RANDOM_GENESET_NULL_SUMMARY.tsv",
    "49b_DISCOVERY_REGULATOR_PROJECTION.tsv",
]:

    print(
        OUT / p
    )

