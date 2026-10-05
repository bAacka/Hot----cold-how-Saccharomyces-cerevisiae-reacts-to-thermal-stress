#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/")

RESP = (
    ROOT
    / "14_response_programs"
    / "response_classification.tsv"
)

TRAJ = (
    ROOT
    / "14_response_programs"
    / "response_class_trajectories_long.tsv.gz"
)

NET = (
    ROOT
    / "16_networkin_response_discovery"
    / "networkin_predictions_COLLAPSED.tsv"
)

AB = (
    ROOT
    / "07_regulator_validation"
    / "FINAL_EXPERIMENTAL_AB_kinase_edges.tsv"
)

OUT = (
    ROOT
    / "17_networkin_temporal_discovery"
)
OUT.mkdir(parents=True, exist_ok=True)


PKA_GROUP = "TPK2_TPK1_TPK3_group"

PKA_REGULATORS = {
    "TPK1",
    "TPK2",
    "TPK3",
}

MIN_PREDICTED_GROUPS = 5
MIN_TARGET_GENES = 3

N_PERM = 50000

RNG = np.random.default_rng(
    20260914
)


                                                              
         
                                                              

def norm_id(x):

    if pd.isna(x):
        return ""

    if isinstance(x, float) and x.is_integer():
        return str(int(x))

    s = str(x).strip()

    if s.endswith(".0"):
        try:
            return str(
                int(float(s))
            )
        except Exception:
            pass

    return s


def bh(pvalues):

    p = np.asarray(
        pvalues,
        dtype=float,
    )

    q = np.full(
        len(p),
        np.nan,
        dtype=float,
    )

    good = np.isfinite(p)

    if not good.any():
        return q

    pv = p[good]

    order = np.argsort(pv)

    ranked = pv[order]

    m = len(ranked)

    adj = (
        ranked
        * m
        / np.arange(
            1,
            m + 1,
        )
    )

    adj = np.minimum.accumulate(
        adj[::-1]
    )[::-1]

    adj = np.minimum(
        adj,
        1.0,
    )

    restored = np.empty_like(
        adj
    )

    restored[order] = adj

    q[good] = restored

    return q


def permutation_mean_test(
    observed,
    background,
    n_target,
):

    background = np.asarray(
        background,
        dtype=float,
    )

    background = background[
        np.isfinite(
            background
        )
    ]

    if n_target > len(background):
        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )

    null = np.empty(
        N_PERM,
        dtype=float,
    )

    for i in range(N_PERM):

        ix = RNG.choice(
            len(background),
            size=n_target,
            replace=False,
        )

        null[i] = (
            background[ix]
            .mean()
        )

    center = float(
        null.mean()
    )

    sd = float(
        null.std(ddof=1)
    )

    z = (
        (observed - center)
        / sd
        if sd > 0
        else np.nan
    )

    p = (
        1
        + np.sum(
            np.abs(
                null - center
            )
            >= abs(
                observed - center
            )
        )
    ) / (
        N_PERM + 1
    )

    return (
        float(p),
        float(z),
        center,
        sd,
    )


                                                              
      
                                                              

resp = pd.read_csv(
    RESP,
    sep="\t",
)

traj = pd.read_csv(
    TRAJ,
    sep="\t",
    compression="gzip",
)

pred = pd.read_csv(
    NET,
    sep="\t",
)

ab = pd.read_csv(
    AB,
    sep="\t",
)


if len(resp) != 347:
    raise RuntimeError(
        f"Expected 347 dynamic groups; got {len(resp)}"
    )


resp["s1_id"] = (
    resp["id"]
    .map(norm_id)
)

ab["s1_id"] = (
    ab["s1_id"]
    .map(norm_id)
)


                                                              
                                          
                                                              

pka_ab = ab[
    ab["regulator"]
    .isin(
        PKA_REGULATORS
    )
].copy()


pka_ab_sites = set(
    pka_ab.merge(
        resp[
            [
                "s1_id",
                "site_id",
            ]
        ],
        on="s1_id",
        how="inner",
    )[
        "site_id"
    ]
)


                                                              
                           
 
                    
                                                         
 
                              
                                                           
                                         
                                                              

candidate_sets = []


for kinase_group, g in pred.groupby(
    "predicted_kinase_group"
):

    ids = set(
        g["dynamic_group_id"]
    )

    if len(ids) < MIN_PREDICTED_GROUPS:
        continue

    candidate_sets.append({
        "analysis_label":
            kinase_group,

        "predicted_kinase_group":
            kinase_group,

        "analysis_type":
            "standard",

        "known_positive_control":
            kinase_group == PKA_GROUP,

        "site_ids":
            ids,
    })


pka_all = set(
    pred.loc[
        pred[
            "predicted_kinase_group"
        ] == PKA_GROUP,
        "dynamic_group_id",
    ]
)


pka_orthogonal = (
    pka_all
    - pka_ab_sites
)


candidate_sets.append({
    "analysis_label":
        (
            PKA_GROUP
            + "__ORTHOGONAL_NO_CURATED_AB"
        ),

    "predicted_kinase_group":
        PKA_GROUP,

    "analysis_type":
        "orthogonal_positive_control",

    "known_positive_control":
        True,

    "site_ids":
        pka_orthogonal,
})


                                                              
                  
                                                              

membership_rows = []


for candidate in candidate_sets:

    label = candidate[
        "analysis_label"
    ]

    for sid in sorted(
        candidate["site_ids"]
    ):

        r = resp[
            resp["site_id"] == sid
        ]

        if len(r) != 1:
            raise RuntimeError(
                f"Bad response mapping for {sid}"
            )

        r = r.iloc[0]

        membership_rows.append({
            "analysis_label":
                label,

            "analysis_type":
                candidate[
                    "analysis_type"
                ],

            "predicted_kinase_group":
                candidate[
                    "predicted_kinase_group"
                ],

            "site_id":
                sid,

            "gene":
                r["gene_final"],

            "response_class":
                r[
                    "response_class_final"
                ],

            "response_class_origin":
                r[
                    "response_class_origin"
                ],

            "mean_heat_minus_cold":
                r[
                    "mean_heat_minus_cold"
                ],

            "has_curated_PKA_AB":
                sid in pka_ab_sites,
        })


membership = pd.DataFrame(
    membership_rows
)


membership.to_csv(
    OUT
    / "networkin_temporal_candidate_membership.tsv",
    sep="\t",
    index=False,
)


                                                              
                                    
 
                                                              
                                                  
 
                                                        
 
                                                      
                                                              

site_scores = resp[
    [
        "site_id",
        "gene_final",
        "mean_heat_minus_cold",
    ]
].copy()


site_scores[
    "mean_heat_minus_cold"
] = pd.to_numeric(
    site_scores[
        "mean_heat_minus_cold"
    ],
    errors="coerce",
)


site_scores = site_scores.dropna(
    subset=[
        "gene_final",
        "mean_heat_minus_cold",
    ]
)


                            
                                                                
background_gene = (
    site_scores.groupby(
        "gene_final"
    )[
        "mean_heat_minus_cold"
    ]
    .median()
)


                                                              
                      
                                                              

rows = []


for candidate in candidate_sets:

    ids = candidate[
        "site_ids"
    ]

    target_sites = site_scores[
        site_scores[
            "site_id"
        ].isin(ids)
    ].copy()


                                                              
    target_gene = (
        target_sites.groupby(
            "gene_final"
        )[
            "mean_heat_minus_cold"
        ]
        .median()
    )


    n_sites = (
        target_sites[
            "site_id"
        ].nunique()
    )

    n_genes = len(
        target_gene
    )


    if n_genes < MIN_TARGET_GENES:

        rows.append({
            "analysis_label":
                candidate[
                    "analysis_label"
                ],

            "predicted_kinase_group":
                candidate[
                    "predicted_kinase_group"
                ],

            "analysis_type":
                candidate[
                    "analysis_type"
                ],

            "known_positive_control":
                candidate[
                    "known_positive_control"
                ],

            "n_target_sites":
                n_sites,

            "n_target_genes":
                n_genes,

            "status":
                "TOO_FEW_GENES",
        })

        continue


    vals = target_gene.to_numpy(
        dtype=float,
    )


    observed_mean = float(
        vals.mean()
    )

    observed_median = float(
        np.median(vals)
    )


    target_genes = set(
        target_gene.index
    )


                       
                                                
    bg = background_gene[
        ~background_gene.index.isin(
            target_genes
        )
    ]


    p, z, null_mean, null_sd = (
        permutation_mean_test(
            observed_mean,
            bg.to_numpy(
                dtype=float
            ),
            n_genes,
        )
    )


                                              
    loo_means = []

    if n_genes >= 2:

        for i in range(
            n_genes
        ):

            leave = np.delete(
                vals,
                i,
            )

            loo_means.append(
                float(
                    leave.mean()
                )
            )

    same_direction = (
        all(
            np.sign(v)
            == np.sign(
                observed_mean
            )
            for v in loo_means
        )
        if loo_means
        else False
    )


    denominator = float(
        np.abs(vals).sum()
    )

    max_contribution = (
        float(
            np.abs(vals).max()
            / denominator
        )
        if denominator > 0
        else np.nan
    )


    rows.append({
        "analysis_label":
            candidate[
                "analysis_label"
            ],

        "predicted_kinase_group":
            candidate[
                "predicted_kinase_group"
            ],

        "analysis_type":
            candidate[
                "analysis_type"
            ],

        "known_positive_control":
            candidate[
                "known_positive_control"
            ],

        "n_target_sites":
            n_sites,

        "n_target_genes":
            n_genes,

        "mean_heat_minus_cold":
            observed_mean,

        "median_heat_minus_cold":
            observed_median,

        "fraction_heat_lower_than_cold":
            float(
                np.mean(
                    vals < 0
                )
            ),

        "background_gene_mean":
            float(
                bg.mean()
            ),

        "background_gene_median":
            float(
                bg.median()
            ),

        "empirical_z":
            z,

        "empirical_p":
            p,

        "null_mean":
            null_mean,

        "null_sd":
            null_sd,

        "all_leave_one_out_same_direction":
            same_direction,

        "loo_mean_closest_to_zero":
            (
                float(
                    loo_means[
                        np.argmin(
                            np.abs(
                                loo_means
                            )
                        )
                    ]
                )
                if loo_means
                else np.nan
            ),

        "max_single_gene_abs_contribution":
            max_contribution,

        "status":
            "TESTED",
    })


results = pd.DataFrame(
    rows
)


                                                              
                                           
 
                                                          
                                                         
                                                              

results[
    "q_across_standard_groups"
] = np.nan


standard_idx = results.index[
    (
        results[
            "analysis_type"
        ] == "standard"
    )
    & (
        results[
            "status"
        ] == "TESTED"
    )
]


results.loc[
    standard_idx,
    "q_across_standard_groups",
] = bh(
    results.loc[
        standard_idx,
        "empirical_p",
    ]
)


results = results.sort_values(
    [
        "analysis_type",
        "q_across_standard_groups",
        "empirical_p",
        "analysis_label",
    ],
    na_position="last",
)


results.to_csv(
    OUT
    / "networkin_global_heat_cold_tests.tsv",
    sep="\t",
    index=False,
)


                                                              
                
                                                              

discoveries = results[
    (
        results[
            "analysis_type"
        ] == "standard"
    )
    & (
        results[
            "q_across_standard_groups"
        ] < 0.05
    )
    & (
        ~results[
            "known_positive_control"
        ]
    )
].copy()


discoveries.to_csv(
    OUT
    / "DISCOVERY_networkin_temporal_hits.tsv",
    sep="\t",
    index=False,
)


                                                              
                                    
 
                                                  
                                                               
                                                              

h = traj[
    traj["condition"]
    == "heat"
][
    [
        "site_id",
        "gene",
        "time_min",
        "log2FC",
    ]
].rename(
    columns={
        "log2FC":
            "heat_log2FC",
    }
)


c = traj[
    traj["condition"]
    == "cold"
][
    [
        "site_id",
        "gene",
        "time_min",
        "log2FC",
    ]
].rename(
    columns={
        "log2FC":
            "cold_log2FC",
    }
)


delta = h.merge(
    c,
    on=[
        "site_id",
        "gene",
        "time_min",
    ],
    how="inner",
    validate="one_to_one",
)


delta["heat_minus_cold"] = (
    delta["heat_log2FC"]
    - delta["cold_log2FC"]
)


profile_rows = []


for candidate in candidate_sets:

    ids = candidate[
        "site_ids"
    ]


    for time_min, g in delta[
        delta["site_id"]
        .isin(ids)
    ].groupby(
        "time_min"
    ):

        gene_values = (
            g.dropna(
                subset=[
                    "heat_minus_cold"
                ]
            )
            .groupby(
                "gene"
            )[
                "heat_minus_cold"
            ]
            .median()
        )


        if len(gene_values) == 0:
            continue


        vals = gene_values.to_numpy(
            dtype=float,
        )


        profile_rows.append({
            "analysis_label":
                candidate[
                    "analysis_label"
                ],

            "analysis_type":
                candidate[
                    "analysis_type"
                ],

            "predicted_kinase_group":
                candidate[
                    "predicted_kinase_group"
                ],

            "time_min":
                int(
                    time_min
                ),

            "n_genes":
                len(vals),

            "mean_heat_minus_cold":
                float(
                    vals.mean()
                ),

            "median_heat_minus_cold":
                float(
                    np.median(vals)
                ),

            "fraction_heat_lower_than_cold":
                float(
                    np.mean(
                        vals < 0
                    )
                ),
        })


profiles = pd.DataFrame(
    profile_rows
)


profiles.to_csv(
    OUT
    / "networkin_temporal_profiles.tsv",
    sep="\t",
    index=False,
)


                                                              
                                  
                                                              

pka_audit = pd.DataFrame(
    [
        {
            "metric":
                "NetworKIN_PKA_dynamic_groups",

            "value":
                len(pka_all),
        },

        {
            "metric":
                "dynamic_groups_with_curated_TPK_AB",

            "value":
                len(
                    pka_ab_sites
                ),
        },

        {
            "metric":
                "NetworKIN_PKA_groups_overlapping_curated_TPK_AB",

            "value":
                len(
                    pka_all
                    & pka_ab_sites
                ),
        },

        {
            "metric":
                "orthogonal_NetworKIN_PKA_groups_remaining",

            "value":
                len(
                    pka_orthogonal
                ),
        },
    ]
)


pka_audit.to_csv(
    OUT
    / "PKA_ORTHOGONALITY_AUDIT.tsv",
    sep="\t",
    index=False,
)


                                                              
         
                                                              

standard = results[
    (
        results[
            "analysis_type"
        ] == "standard"
    )
    & (
        results[
            "status"
        ] == "TESTED"
    )
].copy()


summary = pd.DataFrame(
    [
        {
            "n_standard_prediction_groups_tested":
                len(standard),

            "n_standard_global_q05":
                int(
                    (
                        standard[
                            "q_across_standard_groups"
                        ] < 0.05
                    ).sum()
                ),

            "n_novel_discovery_hits":
                len(
                    discoveries
                ),
        }
    ]
)


summary.to_csv(
    OUT
    / "NETWORKIN_TEMPORAL_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

print()
print("=" * 80)
print("NETWOR-KIN GLOBAL TEMPORAL DIFFERENTIAL COMPLETE")
print("=" * 80)

print()
print("PKA ORTHOGONALITY")
print("-" * 80)

print(
    pka_audit.to_string(
        index=False
    )
)

print()
print("SUMMARY")
print("-" * 80)

print(
    summary.to_string(
        index=False
    )
)

print()
print("STANDARD TESTS")
print("-" * 80)

print(
    standard[
        [
            "predicted_kinase_group",
            "n_target_sites",
            "n_target_genes",
            "mean_heat_minus_cold",
            "fraction_heat_lower_than_cold",
            "empirical_z",
            "empirical_p",
            "q_across_standard_groups",
            "all_leave_one_out_same_direction",
            "max_single_gene_abs_contribution",
            "known_positive_control",
        ]
    ]
    .sort_values(
        [
            "q_across_standard_groups",
            "empirical_p",
        ]
    )
    .to_string(
        index=False
    )
)

print()
print("ORTHOGONAL PKA CONTROL")
print("-" * 80)

orth = results[
    results[
        "analysis_type"
    ]
    == "orthogonal_positive_control"
]

print(
    orth.to_string(
        index=False
    )
)

print()
print("NOVEL DISCOVERY HITS")
print("-" * 80)

if len(discoveries):

    print(
        discoveries.to_string(
            index=False
        )
    )

else:
    print("NONE")

print()
print("Output:", OUT)
