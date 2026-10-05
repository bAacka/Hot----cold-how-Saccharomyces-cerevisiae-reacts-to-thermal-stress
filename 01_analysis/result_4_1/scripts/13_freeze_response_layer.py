#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd


ROOT = Path("/")

S1FILE = ROOT / "00_raw" / "Kanshin_2015_Table_S1.xlsx"
S2FILE = ROOT / "00_raw" / "Kanshin_2015_Table_S2.xlsx"

OUT = ROOT / "14_response_programs"
OUT.mkdir(parents=True, exist_ok=True)

TIMES = np.arange(0, 30, 2)


                                                              
                   
 
            
                                        
                                                           
                                                              
                                                    
                                                              

CURATED_HEAT = {
    "HOG1_Y176",
    "MSN2_S288",
    "MSN2_S633",
    "MSN4_S263",
    "MSN4_S316",
    "MSN4_S558",
}


def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def norm(x):
    return re.sub(
        r"\s+",
        "",
        clean(x).upper(),
    )


def psite_only(x):
    """
    S2 pSite:
        ABP1_T165/S169

    S1 pSites:
        T165/S169
    """
    s = clean(x)
    return re.sub(
        r"^[^_]+_",
        "",
        s,
    )


def site_gene(x):
    s = clean(x)

    if "_" not in s:
        return ""

    return s.split("_", 1)[0]


def direction(x, tol=1e-12):

    if not np.isfinite(x):
        return "missing"

    if x > tol:
        return "increase"

    if x < -tol:
        return "decrease"

    return "zero"


def temporal_features(row, prefix):

    vals = np.asarray(
        [
            pd.to_numeric(
                row.get(
                    f"{prefix}.T{t:02d}",
                    np.nan,
                ),
                errors="coerce",
            )
            for t in TIMES
        ],
        dtype=float,
    )

    good = np.isfinite(vals)

    rec = {
        "n_observed":
            int(good.sum()),
    }

    if not good.any():

        rec.update({
            "peak_time_min": np.nan,
            "peak_log2FC": np.nan,
            "mean_log2FC": np.nan,
            "median_log2FC": np.nan,
            "signed_auc": np.nan,
            "endpoint_log2FC": np.nan,
        })

        return rec

    gtimes = TIMES[good]
    gvals = vals[good]

    peak_i = np.argmax(
        np.abs(gvals)
    )

    if len(gvals) >= 2:
        auc = float(
            np.trapezoid(
                gvals,
                gtimes,
            )
        )
    else:
        auc = np.nan

    rec.update({
        "peak_time_min":
            int(
                gtimes[peak_i]
            ),

        "peak_log2FC":
            float(
                gvals[peak_i]
            ),

        "mean_log2FC":
            float(
                np.mean(gvals)
            ),

        "median_log2FC":
            float(
                np.median(gvals)
            ),

        "signed_auc":
            auc,

        "endpoint_log2FC":
            float(
                gvals[-1]
            ),
    })

    return rec


                                                              
                       
                                                              

s1 = pd.read_excel(
    S1FILE,
    sheet_name="Phosphopeptides",
)

s2 = pd.read_excel(
    S2FILE,
    sheet_name="Dynamic_pSites",
)

assert len(s1) == 2777
assert len(s2) == 347


                                                              
                        
                                                              

s1["join_key"] = (
    s1["Uniprot"].map(norm)
    + "|"
    + s1["pSites"].map(norm)
)

s2["psite_only"] = (
    s2["pSite"]
    .map(psite_only)
)

s2["gene_S2"] = (
    s2["pSite"]
    .map(site_gene)
)

s2["join_key"] = (
    s2["Uniprot"].map(norm)
    + "|"
    + s2["psite_only"].map(norm)
)


if s1["join_key"].duplicated().any():

    dup = s1[
        s1["join_key"].duplicated(
            keep=False
        )
    ][
        [
            "join_key",
            "Gene",
            "pSites",
        ]
    ]

    print(dup.to_string(index=False))

    raise RuntimeError(
        "S1 join key is not unique"
    )


                         
keep = [
    "join_key",
    "id",
    "Gene",
    "Standard Name",
    "Behaviour",
    "FitClust_Heat",
    "FitClust_Cold",
]

for prefix in ["Heat", "Cold"]:
    for t in TIMES:
        keep.append(
            f"{prefix}.T{t:02d}"
        )


x = s2.merge(
    s1[keep],
    on="join_key",
    how="left",
    validate="one_to_one",
    indicator=True,
)


if not (
    x["_merge"] == "both"
).all():

    bad = x[
        x["_merge"] != "both"
    ]

    print(
        bad[
            [
                "Uniprot",
                "pSite",
                "join_key",
            ]
        ].to_string(index=False)
    )

    raise RuntimeError(
        "S2 -> S1 mapping incomplete"
    )


x = x.drop(
    columns="_merge"
)


                                                              
                 
                                                              

x["site_id"] = (
    x["gene_S2"]
    + "_"
    + x["psite_only"]
)

x["gene_final"] = np.where(
    x["Gene"].notna()
    & x["Gene"]
        .astype(str)
        .str.strip()
        .ne(""),
    x["Gene"].astype(str),
    x["gene_S2"],
)


                                                              
                                    
                                                              

x["behaviour_source"] = (
    x["Behaviour"]
    .fillna("")
    .astype(str)
    .str.strip()
)


VALID_SOURCE = {
    "Heat",
    "Cold",
    "Bidirectional",
    "Temp_Independent",
    "",
}

unexpected = sorted(
    set(
        x["behaviour_source"]
    )
    - VALID_SOURCE
)

if unexpected:
    raise RuntimeError(
        "Unexpected Behaviour values: "
        + repr(unexpected)
    )


                                                              
                      
 
                             
                                                        
                                                              

def final_class(row):

    b = row["behaviour_source"]

    if b:
        return b

    if row["site_id"] in CURATED_HEAT:
        return "Heat"

    return "UNRESOLVED"


x["response_class_final"] = (
    x.apply(
        final_class,
        axis=1,
    )
)


x["response_class_origin"] = np.where(
    x["behaviour_source"] != "",
    "source_Behaviour",
    np.where(
        x["site_id"].isin(
            CURATED_HEAT
        ),
        "manual_heat_candidate",
        "unresolved",
    ),
)


unresolved = x[
    x["response_class_final"]
    == "UNRESOLVED"
].copy()

unresolved.to_csv(
    OUT / "UNRESOLVED_response_classes.tsv",
    sep="\t",
    index=False,
)

if len(unresolved):
    raise RuntimeError(
        f"{len(unresolved)} response classes remain unresolved"
    )


                                                              
                      
                                                              

SEMANTICS = {
    "Heat":
        "heat_specific",

    "Cold":
        "cold_specific",

    "Bidirectional":
        "shared_opposite_direction",

    "Temp_Independent":
        "shared_same_direction",
}

x["response_semantics"] = (
    x["response_class_final"]
    .map(SEMANTICS)
)


x["heat_response_final"] = (
    x["response_class_final"]
    .isin(
        [
            "Heat",
            "Bidirectional",
            "Temp_Independent",
        ]
    )
)

x["cold_response_final"] = (
    x["response_class_final"]
    .isin(
        [
            "Cold",
            "Bidirectional",
            "Temp_Independent",
        ]
    )
)


                                                              
              
                                                           
                                                              

x["heat_kinetic_cluster"] = (
    pd.to_numeric(
        x["ClusterID_Heat"],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
)

x["cold_kinetic_cluster"] = (
    pd.to_numeric(
        x["ClusterID_Cold"],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
)

x["heat_kinetic_cluster_member"] = (
    x["heat_kinetic_cluster"] > 0
)

x["cold_kinetic_cluster_member"] = (
    x["cold_kinetic_cluster"] > 0
)


                                                              
                                
                                                              

feature_rows = []

for _, row in x.iterrows():

    hf = temporal_features(
        row,
        "Heat",
    )

    cf = temporal_features(
        row,
        "Cold",
    )

    paired_h = []
    paired_c = []

    for t in TIMES:

        h = pd.to_numeric(
            row.get(
                f"Heat.T{t:02d}",
                np.nan,
            ),
            errors="coerce",
        )

        c = pd.to_numeric(
            row.get(
                f"Cold.T{t:02d}",
                np.nan,
            ),
            errors="coerce",
        )

        if np.isfinite(h) and np.isfinite(c):

            paired_h.append(h)
            paired_c.append(c)


    if len(paired_h) >= 3:

        paired_h = np.asarray(
            paired_h,
            dtype=float,
        )

        paired_c = np.asarray(
            paired_c,
            dtype=float,
        )

        if (
            np.std(paired_h) > 0
            and np.std(paired_c) > 0
        ):

            corr = float(
                np.corrcoef(
                    paired_h,
                    paired_c,
                )[0, 1]
            )

        else:
            corr = np.nan

        mean_delta = float(
            np.mean(
                paired_h
                - paired_c
            )
        )

    else:

        corr = np.nan
        mean_delta = np.nan


    feature_rows.append({
        "site_id":
            row["site_id"],

        "heat_n_observed":
            hf["n_observed"],

        "heat_peak_time_min":
            hf["peak_time_min"],

        "heat_peak_log2FC":
            hf["peak_log2FC"],

        "heat_mean_log2FC":
            hf["mean_log2FC"],

        "heat_median_log2FC":
            hf["median_log2FC"],

        "heat_signed_auc":
            hf["signed_auc"],

        "heat_endpoint_log2FC":
            hf["endpoint_log2FC"],

        "cold_n_observed":
            cf["n_observed"],

        "cold_peak_time_min":
            cf["peak_time_min"],

        "cold_peak_log2FC":
            cf["peak_log2FC"],

        "cold_mean_log2FC":
            cf["mean_log2FC"],

        "cold_median_log2FC":
            cf["median_log2FC"],

        "cold_signed_auc":
            cf["signed_auc"],

        "cold_endpoint_log2FC":
            cf["endpoint_log2FC"],

        "heat_cold_paired_correlation":
            corr,

        "mean_heat_minus_cold":
            mean_delta,
    })


features = pd.DataFrame(
    feature_rows
)

x = x.merge(
    features,
    on="site_id",
    how="left",
    validate="one_to_one",
)


                                                              
                                 
 
                                                            
                                                              

x["heat_auc_direction"] = (
    x["heat_signed_auc"]
    .map(direction)
)

x["cold_auc_direction"] = (
    x["cold_signed_auc"]
    .map(direction)
)


def shared_direction_qc(row):

    cls = row[
        "response_class_final"
    ]

    h = row[
        "heat_auc_direction"
    ]

    c = row[
        "cold_auc_direction"
    ]

    if (
        h in {"missing", "zero"}
        or c in {"missing", "zero"}
    ):
        return "not_testable"

    same = h == c

    if cls == "Bidirectional":
        return (
            "consistent"
            if not same
            else "discordant"
        )

    if cls == "Temp_Independent":
        return (
            "consistent"
            if same
            else "discordant"
        )

    return "not_applicable"


x["shared_direction_qc"] = (
    x.apply(
        shared_direction_qc,
        axis=1,
    )
)


                                                              
                                           
                                                              

final_cols = [
    "site_id",
    "gene_final",
    "Uniprot",
    "psite_only",
    "id",

    "behaviour_source",
    "response_class_final",
    "response_class_origin",
    "response_semantics",

    "heat_response_final",
    "cold_response_final",

    "heat_kinetic_cluster",
    "cold_kinetic_cluster",
    "heat_kinetic_cluster_member",
    "cold_kinetic_cluster_member",

    "FitClust_Heat",
    "FitClust_Cold",

    "heat_n_observed",
    "cold_n_observed",

    "heat_peak_time_min",
    "heat_peak_log2FC",
    "heat_signed_auc",
    "heat_auc_direction",

    "cold_peak_time_min",
    "cold_peak_log2FC",
    "cold_signed_auc",
    "cold_auc_direction",

    "heat_cold_paired_correlation",
    "mean_heat_minus_cold",

    "shared_direction_qc",
]


classification = x[
    final_cols
].copy()


classification.to_csv(
    OUT / "response_classification.tsv",
    sep="\t",
    index=False,
)


                                                              
                       
                                                              

long_rows = []

for _, row in x.iterrows():

    for condition, prefix in [
        ("heat", "Heat"),
        ("cold", "Cold"),
    ]:

        for t in TIMES:

            value = pd.to_numeric(
                row.get(
                    f"{prefix}.T{t:02d}",
                    np.nan,
                ),
                errors="coerce",
            )

            long_rows.append({
                "site_id":
                    row["site_id"],

                "gene":
                    row["gene_final"],

                "psite_group":
                    row["psite_only"],

                "response_class":
                    row[
                        "response_class_final"
                    ],

                "response_semantics":
                    row[
                        "response_semantics"
                    ],

                "condition":
                    condition,

                "time_min":
                    int(t),

                "log2FC":
                    value,
            })


long = pd.DataFrame(
    long_rows
)

long.to_csv(
    OUT
    / "response_class_trajectories_long.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


                                                              
              
                                                              

counts = (
    classification.groupby(
        [
            "response_class_final",
            "response_semantics",
            "response_class_origin",
        ]
    )
    .size()
    .rename("n_groups")
    .reset_index()
    .sort_values(
        [
            "response_class_final",
            "response_class_origin",
        ]
    )
)

counts.to_csv(
    OUT / "response_class_counts.tsv",
    sep="\t",
    index=False,
)


                                                              
                          
                                                              

class_summary = (
    classification.groupby(
        [
            "response_class_final",
            "response_semantics",
        ]
    )
    .agg(
        n_groups=(
            "site_id",
            "size",
        ),

        n_genes=(
            "gene_final",
            "nunique",
        ),

        median_heat_peak_log2FC=(
            "heat_peak_log2FC",
            "median",
        ),

        median_cold_peak_log2FC=(
            "cold_peak_log2FC",
            "median",
        ),

        median_heat_auc=(
            "heat_signed_auc",
            "median",
        ),

        median_cold_auc=(
            "cold_signed_auc",
            "median",
        ),

        median_heat_cold_corr=(
            "heat_cold_paired_correlation",
            "median",
        ),

        median_heat_minus_cold=(
            "mean_heat_minus_cold",
            "median",
        ),
    )
    .reset_index()
)

class_summary.to_csv(
    OUT / "response_class_summary.tsv",
    sep="\t",
    index=False,
)


                                                              
                                             
                                                              

shared = classification[
    classification[
        "response_class_final"
    ].isin(
        [
            "Bidirectional",
            "Temp_Independent",
        ]
    )
].copy()


direction_qc = (
    shared.groupby(
        [
            "response_class_final",
            "shared_direction_qc",
        ]
    )
    .size()
    .rename("n_groups")
    .reset_index()
)

direction_qc.to_csv(
    OUT / "shared_direction_consistency.tsv",
    sep="\t",
    index=False,
)


discordant = shared[
    shared["shared_direction_qc"]
    == "discordant"
].copy()

discordant.to_csv(
    OUT / "shared_direction_DISCORDANT.tsv",
    sep="\t",
    index=False,
)


                                                              
                                           
                                                              

with open(
    OUT / "DEPRECATED_INTERPRETATIONS.txt",
    "w",
) as fh:

    fh.write(
        "DO NOT USE:\n"
        "-----------\n"
        "1. ClusterID_Heat > 0 as a biological "
        "'heat-response' flag.\n"
        "2. ClusterID_Cold > 0 as a biological "
        "'cold-response' flag.\n"
        "3. The previous 282 heat / 80 cold counts as "
        "response counts.\n"
        "4. The previous 15-site heat/cold overlap as the "
        "set of sites responding to both temperatures.\n"
        "5. The previous 6 same / 9 opposite result derived "
        "from that 15-site intersection.\n\n"

        "USE INSTEAD:\n"
        "------------\n"
        "S1 Behaviour = biological temperature-response class.\n"
        "S2 ClusterID_Heat/Cold = kinetic fuzzy-cluster "
        "membership only.\n"
        "Raw trajectories = quantitative amplitude/timing layer.\n"
    )


                                                              
                  
                                                              

expected_source = {
    "Heat": 157,
    "Cold": 4,
    "Bidirectional": 124,
    "Temp_Independent": 56,
    "": 6,
}

observed_source = (
    x["behaviour_source"]
    .value_counts()
    .to_dict()
)

for key, val in expected_source.items():

    if observed_source.get(
        key,
        0,
    ) != val:

        raise RuntimeError(
            f"Unexpected source Behaviour count "
            f"{key!r}: "
            f"{observed_source.get(key,0)} "
            f"!= {val}"
        )


expected_final = {
    "Heat": 163,
    "Cold": 4,
    "Bidirectional": 124,
    "Temp_Independent": 56,
}

observed_final = (
    x["response_class_final"]
    .value_counts()
    .to_dict()
)

if observed_final != expected_final:
    raise RuntimeError(
        "Unexpected final response counts:\n"
        + repr(observed_final)
    )


                                                              
       
                                                              

print()
print("=" * 80)
print("AUTHORITATIVE RESPONSE LAYER FROZEN")
print("=" * 80)

print()
print("SOURCE BEHAVIOUR")
print("-" * 80)

print(
    pd.Series(
        observed_source
    )
    .rename("n")
    .to_string()
)

print()
print("FINAL RESPONSE CLASSES")
print("-" * 80)

print(
    pd.Series(
        observed_final
    )
    .rename("n")
    .to_string()
)

print()
print("CLASS SUMMARY")
print("-" * 80)

print(
    class_summary.to_string(
        index=False
    )
)

print()
print("SHARED-DIRECTION QC")
print("-" * 80)

print(
    direction_qc.to_string(
        index=False
    )
)

print()
print(
    "Shared-class discordant by signed AUC:",
    len(discordant),
)

print()
print("Output:", OUT)
