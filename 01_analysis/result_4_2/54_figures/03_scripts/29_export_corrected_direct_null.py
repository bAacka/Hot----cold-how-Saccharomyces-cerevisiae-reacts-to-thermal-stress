#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path("/")

ALL = (
    ROOT
    / "26_pka_global_trajectory"
    / "ALL_S1_heat_cold_AUC.tsv.gz"
)

TARGETS = (
    ROOT
    / "31_uliana_context_interaction"
    / "DIRECT_A_exact_targets.tsv"
)

OUT = (
    ROOT
    / "32_pka_direct_site_trajectory_corrected"
)
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260915)

N_PERM = 100000
N_LOO = 50000


def empirical_lower(obs, null):
    null = np.asarray(null, dtype=float)
    return float(
        (1 + np.sum(null <= obs))
        / (len(null) + 1)
    )


def empirical_two_sided(obs, null):
    null = np.asarray(null, dtype=float)

    center = np.mean(null)

    dist_obs = abs(obs - center)
    dist_null = np.abs(null - center)

    return float(
        (1 + np.sum(dist_null >= dist_obs))
        / (len(null) + 1)
    )


                                                              
                               
                                                              

all_auc = pd.read_csv(
    ALL,
    sep="\t",
    compression="gzip",
)

all_auc["s1_id"] = (
    all_auc["s1_id"]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
)

all_auc["gene"] = (
    all_auc["gene"]
    .where(all_auc["gene"].notna(), "")
    .astype(str)
    .str.strip()
)

all_auc["heat_minus_cold_auc"] = pd.to_numeric(
    all_auc["heat_minus_cold_auc"],
    errors="coerce",
)

all_auc = all_auc[
    (all_auc["gene"] != "")
    & np.isfinite(
        all_auc["heat_minus_cold_auc"]
    )
].copy()


if all_auc["s1_id"].duplicated().any():
    raise RuntimeError(
        "ALL_S1 AUC table contains duplicate s1_id"
    )


                                                              
                           
                                                              

targets = pd.read_csv(
    TARGETS,
    sep="\t",
    dtype={"s1_id": str},
)

targets["s1_id"] = (
    targets["s1_id"]
    .str.replace(r"\.0$", "", regex=True)
)

direct_group_ids = (
    targets["s1_id"]
    .drop_duplicates()
    .tolist()
)

if len(direct_group_ids) != 10:
    raise RuntimeError(
        f"Expected 10 direct-A measurement groups; "
        f"found {len(direct_group_ids)}"
    )


direct = all_auc[
    all_auc["s1_id"].isin(
        direct_group_ids
    )
].copy()


if len(direct) != 10:
    raise RuntimeError(
        f"Expected 10 direct-A AUC rows; "
        f"found {len(direct)}"
    )


missing = (
    set(direct_group_ids)
    - set(direct["s1_id"])
)

if missing:
    raise RuntimeError(
        f"Missing direct-A s1_ids: {sorted(missing)}"
    )


                                                              
                                                    
                                                              

direct_gene = (
    direct
    .groupby("gene", as_index=False)
    .agg(
        direct_A_n_groups=(
            "s1_id",
            "size",
        ),
        direct_A_gene_AUC=(
            "heat_minus_cold_auc",
            "median",
        ),
    )
)


if len(direct_gene) != 8:
    raise RuntimeError(
        f"Expected 8 direct-A genes; got {len(direct_gene)}"
    )


multiplicities = sorted(
    direct_gene[
        "direct_A_n_groups"
    ].astype(int).tolist(),
    reverse=True,
)


print("DIRECT-A MULTIPLICITY PATTERN:")
print(multiplicities)


                                                              
            
 
                                                 
                                                    
                                                        
                                                  
 
                   
                    
           
                          
                       
                                                              

pka_genes = set(
    direct_gene["gene"]
)

bg = all_auc[
    ~all_auc["gene"].isin(
        pka_genes
    )
].copy()


bg_by_gene = {
    gene: grp[
        "heat_minus_cold_auc"
    ].to_numpy(dtype=float)
    for gene, grp in bg.groupby("gene")
}


eligible = {
    m: np.array(
        [
            gene
            for gene, vals in bg_by_gene.items()
            if len(vals) >= m
        ],
        dtype=object,
    )
    for m in set(multiplicities)
}


for m in sorted(eligible):
    print(
        f"Background genes with >= {m} groups:",
        len(eligible[m]),
    )


def draw_pseudomodule(mults):
    chosen = set()
    gene_scores = []

    for m in mults:

        pool = eligible[m]

        while True:
            gene = pool[
                RNG.integers(
                    0,
                    len(pool),
                )
            ]

            if gene not in chosen:
                break

        chosen.add(gene)

        vals = bg_by_gene[gene]

        if m == 1:
            picked = np.array(
                [
                    vals[
                        RNG.integers(
                            0,
                            len(vals),
                        )
                    ]
                ]
            )
        else:
            idx = RNG.choice(
                len(vals),
                size=m,
                replace=False,
            )

            picked = vals[idx]

        gene_scores.append(
            float(
                np.median(
                    picked
                )
            )
        )

    gene_scores = np.asarray(
        gene_scores,
        dtype=float,
    )

    return (
        float(
            np.mean(
                gene_scores
            )
        ),
        float(
            np.median(
                gene_scores
            )
        ),
        int(
            np.sum(
                gene_scores < 0
            )
        ),
    )


                                                              
                            
                                                              

observed_gene_scores = direct_gene[
    "direct_A_gene_AUC"
].to_numpy(dtype=float)

obs_mean = float(
    np.mean(
        observed_gene_scores
    )
)

obs_median = float(
    np.median(
        observed_gene_scores
    )
)

obs_negative = int(
    np.sum(
        observed_gene_scores < 0
    )
)


                                                              
                                   
                                                              

null_mean = np.empty(
    N_PERM,
    dtype=float,
)

null_median = np.empty(
    N_PERM,
    dtype=float,
)

null_negative = np.empty(
    N_PERM,
    dtype=int,
)


for i in range(N_PERM):
    (
        null_mean[i],
        null_median[i],
        null_negative[i],
    ) = draw_pseudomodule(
        multiplicities
    )



                                                              
                                                      
                                                          
                    
                                                              

FIG_NULL_OUT = Path(
    "/"
    "54_figures/02_curated_inputs/Fig4_2/"
    "DIRECT_A_corrected_empirical_null.tsv.gz"
)

FIG_NULL_OUT.parent.mkdir(parents=True, exist_ok=True)

pd.DataFrame(
    {
        "null_mean_gene_AUC": null_mean,
    }
).to_csv(
    FIG_NULL_OUT,
    sep="\t",
    index=False,
    compression="gzip",
)

print(
    "EXPORTED EXACT CORRECTED NULL:",
    FIG_NULL_OUT,
    "n=",
    len(null_mean),
)

raise SystemExit(0)


mean_p_lower = empirical_lower(
    obs_mean,
    null_mean,
)

mean_p_two = empirical_two_sided(
    obs_mean,
    null_mean,
)

median_p_lower = empirical_lower(
    obs_median,
    null_median,
)

negative_p_upper = float(
    (
        1
        + np.sum(
            null_negative
            >= obs_negative
        )
    )
    / (
        N_PERM
        + 1
    )
)


mean_z = float(
    (
        obs_mean
        - np.mean(
            null_mean
        )
    )
    / np.std(
        null_mean,
        ddof=1,
    )
)


                                                              
                               
                                                              

loo_rows = []


for omitted in direct_gene["gene"]:

    remain = direct_gene[
        direct_gene["gene"] != omitted
    ].copy()

    obs = float(
        remain[
            "direct_A_gene_AUC"
        ].mean()
    )

    mults = sorted(
        remain[
            "direct_A_n_groups"
        ].astype(int).tolist(),
        reverse=True,
    )

    null = np.empty(
        N_LOO,
        dtype=float,
    )

    for i in range(N_LOO):

        scores = []
        chosen = set()

        for m in mults:

            pool = eligible[m]

            while True:
                gene = pool[
                    RNG.integers(
                        0,
                        len(pool),
                    )
                ]

                if gene not in chosen:
                    break

            chosen.add(gene)

            vals = bg_by_gene[gene]

            if m == 1:
                picked = [
                    vals[
                        RNG.integers(
                            0,
                            len(vals),
                        )
                    ]
                ]
            else:
                picked = vals[
                    RNG.choice(
                        len(vals),
                        size=m,
                        replace=False,
                    )
                ]

            scores.append(
                np.median(
                    picked
                )
            )

        null[i] = np.mean(
            scores
        )

    loo_rows.append({
        "omitted_gene":
            omitted,

        "n_remaining_genes":
            len(remain),

        "remaining_mean_AUC":
            obs,

        "empirical_lower_p":
            empirical_lower(
                obs,
                null,
            ),
    })


loo = pd.DataFrame(
    loo_rows
).sort_values(
    "empirical_lower_p",
    ascending=False,
)


                                                              
                    
                                                              

abs_total = float(
    np.sum(
        np.abs(
            observed_gene_scores
        )
    )
)


direct_gene[
    "negative"
] = (
    direct_gene[
        "direct_A_gene_AUC"
    ] < 0
)


direct_gene[
    "absolute_contribution"
] = (
    np.abs(
        direct_gene[
            "direct_A_gene_AUC"
        ]
    )
    / abs_total
)


direct_gene = direct_gene.sort_values(
    "direct_A_gene_AUC"
)


                                                              
        
                                                              

direct[
    [
        "s1_id",
        "gene",
        "pSites",
        "heat_minus_cold_auc",
    ]
].sort_values(
    [
        "gene",
        "s1_id",
    ]
).to_csv(
    OUT
    / "DIRECT_A_10_GROUP_AUCs.tsv",
    sep="\t",
    index=False,
)


direct_gene.to_csv(
    OUT
    / "DIRECT_A_8_GENE_AUCs.tsv",
    sep="\t",
    index=False,
)


loo.to_csv(
    OUT
    / "DIRECT_A_leave_one_gene_out.tsv",
    sep="\t",
    index=False,
)


summary = pd.DataFrame(
    [
        {
            "n_direct_A_groups":
                len(direct),

            "n_direct_A_genes":
                len(direct_gene),

            "observed_mean_gene_AUC":
                obs_mean,

            "observed_median_gene_AUC":
                obs_median,

            "negative_genes":
                obs_negative,

            "null_mean":
                float(
                    np.mean(
                        null_mean
                    )
                ),

            "null_sd":
                float(
                    np.std(
                        null_mean,
                        ddof=1,
                    )
                ),

            "empirical_z":
                mean_z,

            "mean_empirical_lower_p":
                mean_p_lower,

            "mean_empirical_two_sided_p":
                mean_p_two,

            "median_empirical_lower_p":
                median_p_lower,

            "negative_count_empirical_p":
                negative_p_upper,

            "all_leave_one_out_means_negative":
                bool(
                    (
                        loo[
                            "remaining_mean_AUC"
                        ] < 0
                    ).all()
                ),

            "all_leave_one_out_p_lt_0_05":
                bool(
                    (
                        loo[
                            "empirical_lower_p"
                        ] < 0.05
                    ).all()
                ),

            "worst_leave_one_out_p":
                float(
                    loo[
                        "empirical_lower_p"
                    ].max()
                ),
        }
    ]
)


summary.to_csv(
    OUT
    / "CORRECTED_DIRECT_A_TRAJECTORY_SUMMARY.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 100)
print("CORRECTED DIRECT-A PKA GLOBAL TRAJECTORY TEST")
print("=" * 100)

print()
print("SUMMARY")
print("-" * 100)
print(
    summary.to_string(
        index=False
    )
)

print()
print("DIRECT-A GENE SCORES")
print("-" * 100)
print(
    direct_gene.to_string(
        index=False
    )
)

print()
print("LEAVE-ONE-GENE-OUT")
print("-" * 100)
print(
    loo.to_string(
        index=False
    )
)

print()
print("Output:", OUT)
