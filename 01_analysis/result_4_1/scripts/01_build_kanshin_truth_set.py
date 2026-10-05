#!/usr/bin/env python3

from pathlib import Path
import hashlib
import re

import numpy as np
import pandas as pd


ROOT = Path("/")
RAW = ROOT / "00_raw"
OUT = ROOT / "01_kanshin"
QC = ROOT / "02_qc"

OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

S1 = RAW / "Kanshin_2015_Table_S1.xlsx"
S2 = RAW / "Kanshin_2015_Table_S2.xlsx"


                                                              
         
                                                              

def clean(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in {"nan", "none", "na"}:
        return ""
    return s


def norm(x):
    return re.sub(r"\s+", "", clean(x).upper())


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def numeric(x):
    return pd.to_numeric(x, errors="coerce")


def split_s2_psite(x):
    """
    Example:
        ABP1_T165/S169 -> gene=ABP1, residue_group=T165/S169
        HOG1_Y176      -> gene=HOG1, residue_group=Y176
    """
    x = clean(x)

    if "_" not in x:
        return "", x

    return x.split("_", 1)


def extract_atomic_sites(residue_group):
    """
    Parse explicit S/T/Y residue positions.

    Important:
    This expands notation mechanically but DOES NOT claim that every
    expanded residue represents an independently quantified phosphosite.
    """
    return [
        x.upper()
        for x in re.findall(r"(?i)[STY]\d+", clean(residue_group))
    ]


                                                              
                           
                                                              

source_manifest = pd.DataFrame([
    {
        "file": S1.name,
        "sha256": sha256(S1),
        "role": "Kanshin Supplementary Table S1: quantified phosphopeptides/time courses",
    },
    {
        "file": S2.name,
        "sha256": sha256(S2),
        "role": "Kanshin Supplementary Table S2: author-defined dynamic phosphosite groups",
    },
])

source_manifest.to_csv(
    OUT / "source_manifest.tsv",
    sep="\t",
    index=False,
)


                                                              
                                
                                                              

s1 = pd.read_excel(
    S1,
    sheet_name="Phosphopeptides",
    header=0,
)

s1 = s1.dropna(axis=1, how="all").copy()

expected_s1 = [
    "id",
    "Uniprot",
    "SGD",
    "Standard Name",
    "Gene",
    "Behaviour",
    "pSites",
    "r^2_Heat",
    "r^2_Cold",
    "FitClust_Heat",
    "FitClust_Cold",
]

missing = [x for x in expected_s1 if x not in s1.columns]

if missing:
    raise RuntimeError(
        "Missing expected S1 columns: " + ", ".join(missing)
    )


                     
heat_cols = [
    f"Heat.T{x:02d}"
    for x in range(0, 30, 2)
]

cold_cols = [
    f"Cold.T{x:02d}"
    for x in range(0, 30, 2)
]

for col in heat_cols + cold_cols:
    if col not in s1.columns:
        raise RuntimeError(f"Missing time-course column: {col}")


                          
s1["s1_row"] = np.arange(2, len(s1) + 2)

s1["uniprot"] = s1["Uniprot"].map(clean)
s1["sgd"] = s1["SGD"].map(clean)
s1["orf"] = s1["Standard Name"].map(clean)
s1["gene"] = s1["Gene"].map(clean)
s1["psite_group"] = s1["pSites"].map(clean)
s1["behaviour"] = s1["Behaviour"].map(clean)

s1["cluster_heat"] = numeric(s1["FitClust_Heat"]).fillna(0).astype(int)
s1["cluster_cold"] = numeric(s1["FitClust_Cold"]).fillna(0).astype(int)

s1["join_key"] = (
    s1["uniprot"].map(norm)
    + "|"
    + s1["psite_group"].map(norm)
)


                                                              
                                            
                                                              

s2 = pd.read_excel(
    S2,
    sheet_name="Dynamic_pSites",
    header=0,
)

s2 = s2.dropna(axis=1, how="all").copy()

expected_s2 = [
    "Name",
    "Description",
    "Uniprot",
    "pSite",
    "ClusterID_Heat",
    "ClusterID_Cold",
]

missing = [x for x in expected_s2 if x not in s2.columns]

if missing:
    raise RuntimeError(
        "Missing expected S2 columns: " + ", ".join(missing)
    )

s2["s2_row"] = np.arange(2, len(s2) + 2)

tmp = s2["pSite"].map(split_s2_psite)

s2["gene"] = [x[0] for x in tmp]
s2["psite_group"] = [x[1] for x in tmp]

s2["uniprot"] = s2["Uniprot"].map(clean)

s2["cluster_heat"] = (
    numeric(s2["ClusterID_Heat"])
    .fillna(0)
    .astype(int)
)

s2["cluster_cold"] = (
    numeric(s2["ClusterID_Cold"])
    .fillna(0)
    .astype(int)
)

s2["join_key"] = (
    s2["uniprot"].map(norm)
    + "|"
    + s2["psite_group"].map(norm)
)


                                                              
                                  
                                                 
                                                              

s1_duplicate_keys = (
    s1.groupby("join_key")
      .size()
      .rename("n_s1_rows")
      .reset_index()
)

s1_duplicate_keys = s1_duplicate_keys[
    s1_duplicate_keys["n_s1_rows"] > 1
].copy()

s1_duplicate_keys.to_csv(
    QC / "S1_duplicate_join_keys.tsv",
    sep="\t",
    index=False,
)


join_columns = [
    "join_key",
    "s1_row",
    "id",
    "sgd",
    "orf",
    "gene",
    "behaviour",
    "r^2_Heat",
    "r^2_Cold",
    "cluster_heat",
    "cluster_cold",
] + heat_cols + cold_cols

joined = s2.merge(
    s1[join_columns],
    on="join_key",
    how="left",
    suffixes=("_S2", "_S1"),
    indicator=True,
)

joined.to_csv(
    OUT / "kanshin_dynamic_groups_full.tsv",
    sep="\t",
    index=False,
)


                                                              
                                          
                                                              

unmatched = joined[joined["_merge"] != "both"].copy()

unmatched.to_csv(
    QC / "S2_unmatched_to_S1.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                           
                                                              

cluster_discordant = joined[
    joined["_merge"].eq("both")
    & (
        joined["cluster_heat_S2"].ne(joined["cluster_heat_S1"])
        | joined["cluster_cold_S2"].ne(joined["cluster_cold_S1"])
    )
].copy()

cluster_discordant.to_csv(
    QC / "S1_S2_cluster_discordance.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                         
                                                              

atomic_rows = []

for _, row in joined.iterrows():

    residues = extract_atomic_sites(row["psite_group"])

    for residue in residues:

        atomic_rows.append({
            "s2_row": row["s2_row"],
            "s1_row": row.get("s1_row", ""),
            "gene": row["gene_S2"],
            "uniprot": row["uniprot"],
            "psite_group": row["psite_group"],
            "atomic_residue": residue,
            "atomic_site_id":
                f"{row['gene_S2']}_{residue}",
            "cluster_heat": row["cluster_heat_S2"],
            "cluster_cold": row["cluster_cold_S2"],
            "behaviour_S1": row.get("behaviour", ""),
        })

atomic = pd.DataFrame(atomic_rows)

atomic.to_csv(
    OUT / "kanshin_atomic_site_expansion.tsv",
    sep="\t",
    index=False,
)


                                                              
                                          
                                                              

atomic_counts = (
    atomic.groupby("atomic_site_id")
          .size()
          .rename("n_occurrences")
          .reset_index()
)

duplicate_atomic = atomic_counts[
    atomic_counts["n_occurrences"] > 1
].copy()

duplicate_atomic.to_csv(
    QC / "duplicate_atomic_site_ids.tsv",
    sep="\t",
    index=False,
)

if len(duplicate_atomic):

    duplicated_rows = atomic[
        atomic["atomic_site_id"].isin(
            duplicate_atomic["atomic_site_id"]
        )
    ].sort_values(
        ["atomic_site_id", "s2_row"]
    )

else:
    duplicated_rows = atomic.iloc[0:0].copy()

duplicated_rows.to_csv(
    QC / "duplicate_atomic_site_rows.tsv",
    sep="\t",
    index=False,
)


                                                              
                      
                                                              

multisite_rows = []

for _, row in joined.iterrows():

    residues = extract_atomic_sites(row["psite_group"])

    if len(residues) > 1:
        multisite_rows.append({
            "s2_row": row["s2_row"],
            "gene": row["gene_S2"],
            "uniprot": row["uniprot"],
            "psite_group": row["psite_group"],
            "n_explicit_residues": len(residues),
            "residues": "|".join(residues),
            "cluster_heat": row["cluster_heat_S2"],
            "cluster_cold": row["cluster_cold_S2"],
        })

multisite = pd.DataFrame(multisite_rows)

multisite.to_csv(
    QC / "multisite_groups.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                            
                                                              

behaviour_counts = (
    s1["behaviour"]
    .replace("", pd.NA)
    .value_counts(dropna=False)
    .rename_axis("behaviour")
    .reset_index(name="n")
)

behaviour_counts.to_csv(
    QC / "S1_behaviour_counts.tsv",
    sep="\t",
    index=False,
)


s2_without_behaviour = joined[
    joined["behaviour"].map(clean).eq("")
].copy()

s2_without_behaviour.to_csv(
    QC / "S2_dynamic_without_S1_behaviour.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                                
 
            
                                           
                                                
                                                              

def peak_from_row(row, cols):

    vals = pd.to_numeric(
        row[cols],
        errors="coerce",
    )

    vals = vals.dropna()

    if vals.empty:
        return pd.Series({
            "peak_time": np.nan,
            "peak_log2FC": np.nan,
            "peak_abs_log2FC": np.nan,
            "peak_direction": "missing",
        })

    abs_vals = vals.abs()

    max_abs = abs_vals.max()

                                 
                                      
    candidates = abs_vals[
        np.isclose(abs_vals, max_abs)
    ]

    chosen_col = candidates.index[0]
    signed_value = float(vals[chosen_col])

    time_match = re.search(
        r"T(\d+)$",
        chosen_col
    )

    time = (
        int(time_match.group(1))
        if time_match
        else np.nan
    )

    if signed_value > 0:
        direction = "increase"
    elif signed_value < 0:
        direction = "decrease"
    else:
        direction = "zero"

    return pd.Series({
        "peak_time": time,
        "peak_log2FC": signed_value,
        "peak_abs_log2FC": abs(signed_value),
        "peak_direction": direction,
    })


heat_peak = joined.apply(
    lambda r: peak_from_row(r, heat_cols),
    axis=1,
).add_prefix("heat_")

cold_peak = joined.apply(
    lambda r: peak_from_row(r, cold_cols),
    axis=1,
).add_prefix("cold_")

truth = pd.concat(
    [
        joined.reset_index(drop=True),
        heat_peak.reset_index(drop=True),
        cold_peak.reset_index(drop=True),
    ],
    axis=1,
)


                                                              
                              
 
                                                        
             
                      
                        
             
                                                              

truth["responds_heat_by_S2_cluster"] = (
    truth["cluster_heat_S2"] > 0
)

truth["responds_cold_by_S2_cluster"] = (
    truth["cluster_cold_S2"] > 0
)

truth["measured_heat_direction"] = (
    truth["heat_peak_direction"]
)

truth["measured_cold_direction"] = (
    truth["cold_peak_direction"]
)


                                                              
                           
                                                              

truth["dynamic_group_id"] = (
    truth["gene_S2"].map(norm)
    + "_"
    + truth["psite_group"].map(norm)
)

truth["protein_id"] = truth["uniprot"].map(norm)


                                                              
                             
                                                              

canonical_cols = [
    "dynamic_group_id",
    "protein_id",
    "gene_S2",
    "uniprot",
    "psite_group",
    "s2_row",
    "s1_row",
    "id",
    "sgd",
    "orf",
    "Name",
    "Description",

    "cluster_heat_S2",
    "cluster_cold_S2",

    "behaviour",

    "heat_peak_time",
    "heat_peak_log2FC",
    "heat_peak_abs_log2FC",
    "heat_peak_direction",

    "cold_peak_time",
    "cold_peak_log2FC",
    "cold_peak_abs_log2FC",
    "cold_peak_direction",

    "responds_heat_by_S2_cluster",
    "responds_cold_by_S2_cluster",

    "r^2_Heat",
    "r^2_Cold",
] + heat_cols + cold_cols

canonical = truth[
    [c for c in canonical_cols if c in truth.columns]
].copy()

canonical.to_csv(
    OUT / "kanshin_dynamic_groups_CANONICAL.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


                                                              
                       
 
                              
                                                        
                                                              

protein_summary = (
    canonical.groupby(
        ["protein_id", "gene_S2", "uniprot"],
        dropna=False,
    )
    .agg(
        n_dynamic_groups=("dynamic_group_id", "nunique"),
        n_heat_cluster_groups=(
            "responds_heat_by_S2_cluster",
            "sum",
        ),
        n_cold_cluster_groups=(
            "responds_cold_by_S2_cluster",
            "sum",
        ),
    )
    .reset_index()
)

protein_summary.to_csv(
    OUT / "kanshin_dynamic_proteins.tsv",
    sep="\t",
    index=False,
)


                                                              
                
                                                              

summary = [
    ("S1 quantified rows", len(s1)),
    ("S2 dynamic groups", len(s2)),
    (
        "S2 groups successfully matched to S1",
        int((joined["_merge"] == "both").sum()),
    ),
    (
        "S2 groups unmatched to S1",
        len(unmatched),
    ),
    (
        "S1/S2 cluster-discordant groups",
        len(cluster_discordant),
    ),
    (
        "S2 unique UniProt accessions",
        s2["uniprot"].nunique(),
    ),
    (
        "S2 unique parsed gene labels",
        s2["gene"].nunique(),
    ),
    (
        "S2 multi-site groups",
        len(multisite),
    ),
    (
        "Atomic residue occurrences after mechanical expansion",
        len(atomic),
    ),
    (
        "Unique gene+residue IDs after mechanical expansion",
        atomic["atomic_site_id"].nunique(),
    ),
    (
        "Duplicated gene+residue IDs",
        len(duplicate_atomic),
    ),
    (
        "S2 dynamic groups lacking Behaviour in S1",
        len(s2_without_behaviour),
    ),
]

summary_df = pd.DataFrame(
    summary,
    columns=["metric", "value"],
)

summary_df.to_csv(
    QC / "QC_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
                     
 
                                                          
                                              
                                                              

assert len(s1) == 2777
assert len(s2) == 347
assert len(unmatched) == 0


                                                              
                           
                                                              

with open(QC / "README_QC.txt", "w") as fh:

    fh.write(
        "Kanshin 2015 phosphoproteomics rebuild\n"
        "=======================================\n\n"
    )

    fh.write(
        "AUTHORITATIVE INCLUSION RULE\n"
        "----------------------------\n"
        "A dynamic phosphosite GROUP is included if and only if it "
        "appears in Supplementary Table S2 (Dynamic_pSites).\n\n"
    )

    fh.write(
        "Supplementary Table S1 is used to obtain time-course values, "
        "identifiers and metadata. The S1 Behaviour column is an "
        "annotation and is NOT used as an inclusion criterion.\n\n"
    )

    fh.write(
        "DIRECTION RULE\n"
        "--------------\n"
        "For each condition, the measured peak is the observed time "
        "point with maximum absolute log2 fold-change. The original "
        "signed log2 fold-change at that time point determines "
        "increase/decrease. Absolute value is used only to select the "
        "peak; its sign is never discarded.\n\n"
    )

    fh.write(
        "COUNTING NOTE\n"
        "-------------\n"
        "The official S2 workbook contains 347 dynamic site-group "
        "rows. Multi-site labels are mechanically expanded only for "
        "QC and later structural mapping. We do not yet force the "
        "paper-level 388-site / 271-protein headline because its "
        "counting convention must be explicitly reconstructed.\n\n"
    )

    fh.write(summary_df.to_string(index=False))
    fh.write("\n")


print()
print("=" * 80)
print("KANSHIN REBUILD COMPLETE")
print("=" * 80)
print(summary_df.to_string(index=False))

print()
print("Canonical dynamic groups:")
print(
    OUT / "kanshin_dynamic_groups_CANONICAL.tsv.gz"
)

print()
print("Protein map:")
print(
    OUT / "kanshin_dynamic_proteins.tsv"
)

print()
print("QC:")
print(
    QC / "QC_SUMMARY.tsv"
)
