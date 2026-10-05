#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path("/")

S1FILE = ROOT / "00_raw" / "Kanshin_2015_Table_S1.xlsx"
S2FILE = ROOT / "00_raw" / "Kanshin_2015_Table_S2.xlsx"

OUT = ROOT / "13_response_classification_audit"
OUT.mkdir(parents=True, exist_ok=True)


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


                                                              
                                    
                                                              

s1 = pd.read_excel(
    S1FILE,
    sheet_name="Phosphopeptides",
)

s2 = pd.read_excel(
    S2FILE,
    sheet_name="Dynamic_pSites",
)

print("S1 rows:", len(s1))
print("S2 rows:", len(s2))


                                                              
                        
                                                              

s1["join_key"] = (
    s1["Uniprot"].map(norm)
    + "|"
    + s1["pSites"].map(norm)
)

                                         
                  
                      
             
 
                                                             
s2["psite_only"] = (
    s2["pSite"]
    .astype(str)
    .str.replace(
        r"^[^_]+_",
        "",
        regex=True,
    )
)

s2["join_key"] = (
    s2["Uniprot"].map(norm)
    + "|"
    + s2["psite_only"].map(norm)
)

if s1["join_key"].duplicated().any():
    dup = s1[
        s1["join_key"].duplicated(False)
    ][
        [
            "join_key",
            "Gene",
            "pSites",
        ]
    ]

    print("\nDUPLICATE S1 KEYS")
    print(dup.to_string(index=False))

    raise RuntimeError(
        "S1 join keys not unique"
    )


ann_cols = [
    "join_key",
    "Gene",
    "Behaviour",
]

                                                              
for c in s1.columns:

    lc = str(c).lower()

    if any(
        token in lc
        for token in [
            "fit",
            "clust",
            "r²",
            "r2",
            "member",
            "behav",
        ]
    ):
        if c not in ann_cols:
            ann_cols.append(c)


x = s2.merge(
    s1[ann_cols],
    on="join_key",
    how="left",
    validate="one_to_one",
    indicator=True,
)

unmatched = x[
    x["_merge"] != "both"
].copy()

print(
    "S2 rows mapped to S1:",
    int((x["_merge"] == "both").sum()),
    "/",
    len(x),
)

if len(unmatched):
    print()
    print("UNMATCHED S2 ROWS")
    print(
        unmatched[
            [
                "Name",
                "Uniprot",
                "pSite",
                "join_key",
                "_merge",
            ]
        ].to_string(index=False)
    )
    raise RuntimeError(
        f"{len(unmatched)} S2 rows genuinely failed S1 mapping"
    )

x = x.drop(columns="_merge")


                                                              
                            
                                                              

x["Behaviour_clean"] = (
    x["Behaviour"]
    .fillna("")
    .astype(str)
    .str.strip()
)

print()
print("=" * 80)
print("SOURCE BEHAVIOUR VALUES")
print("=" * 80)

print(
    x["Behaviour_clean"]
    .value_counts(dropna=False)
    .rename_axis("Behaviour")
    .reset_index(name="n")
    .to_string(index=False)
)


                                                              
                           
                                                              

for cond, col in [
    ("heat", "ClusterID_Heat"),
    ("cold", "ClusterID_Cold"),
]:

    vals = pd.to_numeric(
        x[col],
        errors="coerce",
    ).fillna(0)

    x[f"{cond}_cluster_present"] = (
        vals > 0
    )


                                                              
                                           
 
                                                           
                                             
                                                              

def behaviour_flags(b):

    b = clean(b).lower()

    if b == "heat":
        return True, False

    if b == "cold":
        return False, True

    if b in {
        "bidirectional",
        "temp_independent",
        "temperature_independent",
        "temp independent",
        "temperature independent",
    }:
        return True, True

    return False, False


flags = [
    behaviour_flags(b)
    for b in x["Behaviour_clean"]
]

x["heat_response_from_behaviour"] = [
    a for a, b in flags
]

x["cold_response_from_behaviour"] = [
    b for a, b in flags
]

x["behaviour_classified"] = (
    x["Behaviour_clean"] != ""
)


                                                              
                          
                                                              

cluster_counts = pd.DataFrame([
    {
        "classification": "heat_cluster_present",
        "n": int(
            x["heat_cluster_present"].sum()
        ),
    },
    {
        "classification": "cold_cluster_present",
        "n": int(
            x["cold_cluster_present"].sum()
        ),
    },
    {
        "classification": "both_clusters_present",
        "n": int(
            (
                x["heat_cluster_present"]
                & x["cold_cluster_present"]
            ).sum()
        ),
    },
    {
        "classification": "heat_only_cluster_present",
        "n": int(
            (
                x["heat_cluster_present"]
                & ~x["cold_cluster_present"]
            ).sum()
        ),
    },
    {
        "classification": "cold_only_cluster_present",
        "n": int(
            (
                ~x["heat_cluster_present"]
                & x["cold_cluster_present"]
            ).sum()
        ),
    },
])

cluster_counts.to_csv(
    OUT / "cluster_presence_counts.tsv",
    sep="\t",
    index=False,
)


                                                              
                            
                                                              

behaviour_counts = pd.DataFrame([
    {
        "classification": "heat_response",
        "n": int(
            x[
                "heat_response_from_behaviour"
            ].sum()
        ),
    },
    {
        "classification": "cold_response",
        "n": int(
            x[
                "cold_response_from_behaviour"
            ].sum()
        ),
    },
    {
        "classification": "both_responses",
        "n": int(
            (
                x[
                    "heat_response_from_behaviour"
                ]
                & x[
                    "cold_response_from_behaviour"
                ]
            ).sum()
        ),
    },
    {
        "classification": "heat_only_response",
        "n": int(
            (
                x[
                    "heat_response_from_behaviour"
                ]
                & ~x[
                    "cold_response_from_behaviour"
                ]
            ).sum()
        ),
    },
    {
        "classification": "cold_only_response",
        "n": int(
            (
                ~x[
                    "heat_response_from_behaviour"
                ]
                & x[
                    "cold_response_from_behaviour"
                ]
            ).sum()
        ),
    },
    {
        "classification": "unclassified_behaviour",
        "n": int(
            (~x["behaviour_classified"]).sum()
        ),
    },
])

behaviour_counts.to_csv(
    OUT / "behaviour_response_counts.tsv",
    sep="\t",
    index=False,
)


                                                              
                 
                                      
                                                              

x["cluster_presence_class"] = np.select(
    [
        (
            x["heat_cluster_present"]
            & x["cold_cluster_present"]
        ),
        (
            x["heat_cluster_present"]
            & ~x["cold_cluster_present"]
        ),
        (
            ~x["heat_cluster_present"]
            & x["cold_cluster_present"]
        ),
    ],
    [
        "both",
        "heat_only",
        "cold_only",
    ],
    default="neither",
)


cross = pd.crosstab(
    x["Behaviour_clean"],
    x["cluster_presence_class"],
    margins=True,
)

cross.to_csv(
    OUT / "behaviour_vs_cluster_presence.tsv",
    sep="\t",
)


                                                              
                                              
                                                              

x["cluster_heat_equals_response"] = (
    x["heat_cluster_present"]
    == x["heat_response_from_behaviour"]
)

x["cluster_cold_equals_response"] = (
    x["cold_cluster_present"]
    == x["cold_response_from_behaviour"]
)

discordant = x[
    x["behaviour_classified"]
    & (
        ~x["cluster_heat_equals_response"]
        | ~x["cluster_cold_equals_response"]
    )
].copy()

discordant.to_csv(
    OUT / "cluster_vs_behaviour_discordant.tsv",
    sep="\t",
    index=False,
)


                                                              
                           
                                                              

unclassified = x[
    ~x["behaviour_classified"]
].copy()

unclassified.to_csv(
    OUT / "dynamic_groups_without_behaviour.tsv",
    sep="\t",
    index=False,
)


                                                              
                                               
                                                              

special_cols = [
    c for c in x.columns
    if any(
        token in str(c).lower()
        for token in [
            "fit",
            "clust",
            "r²",
            "r2",
            "member",
            "behav",
        ]
    )
]

with open(
    OUT / "source_selection_columns.txt",
    "w",
) as fh:

    for c in special_cols:
        fh.write(str(c) + "\n")


                                                              
                     
                                                              

x.to_csv(
    OUT / "S2_with_S1_response_annotations.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

print()
print("=" * 80)
print("CLUSTER-ID PRESENCE")
print("=" * 80)
print(
    cluster_counts.to_string(
        index=False
    )
)

print()
print("=" * 80)
print("BEHAVIOUR-BASED RESPONSE")
print("=" * 80)
print(
    behaviour_counts.to_string(
        index=False
    )
)

print()
print("=" * 80)
print("BEHAVIOUR x CLUSTER PRESENCE")
print("=" * 80)
print(cross.to_string())

print()
print("=" * 80)
print("UNCLASSIFIED BEHAVIOUR")
print("=" * 80)

if len(unclassified):

    show = [
        c for c in [
            "Name",
            "Gene",
            "pSite",
            "Behaviour_clean",
            "ClusterID_Heat",
            "ClusterID_Cold",
        ]
        if c in unclassified.columns
    ]

    print(
        unclassified[
            show
        ].to_string(index=False)
    )

else:
    print("NONE")

print()
print(
    "Discordant classified rows:",
    len(discordant),
)

print(
    "Output:",
    OUT,
)
