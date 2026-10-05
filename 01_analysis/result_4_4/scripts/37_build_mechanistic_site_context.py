#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path("/")

RES = (
    ROOT
    / "38_pka_mechanistic_resources"
    / "processed_downloads"
)

OUT = ROOT / "39_mechanistic_site_context"
OUT.mkdir(parents=True, exist_ok=True)

NETWORKIN = (
    ROOT
    / "17_networkin_temporal_discovery"
    / "networkin_temporal_candidate_membership.tsv"
)

DIRECT_ATLAS = (
    ROOT
    / "35_phosphoatlas_pka_replication"
    / "DIRECT_A_PKA_matched_groups.tsv"
)

ORTH_ATLAS = (
    ROOT
    / "35_phosphoatlas_pka_replication"
    / "ORTHOGONAL_NETWORKIN_PKA_matched_groups.tsv"
)

ULIANA = (
    RES
    / "PXD052971"
    / "Tukey_table_2400607.csv"
)

STARV = (
    RES
    / "PXD059338"
    / "Phospho_STY_Sites.txt"
)

CDC19_IV = (
    RES
    / "PXD052995"
    / "final_table_241117.csv"
)

ORTH_LABEL = (
    "TPK2_TPK1_TPK3_group"
    "__ORTHOGONAL_NO_CURATED_AB"
)


                                                              
                            
                                                              

DIRECT = [
    ("ATG1",  "S515"),
    ("CDC19", "S22"),
    ("CHO1",  "S46"),
    ("CKI1",  "S85"),
    ("MAF1",  "S90"),
    ("NTH1",  "S83"),
    ("NTH1",  "S20"),
    ("RGT1",  "S284"),
    ("RGT1",  "S202"),
    ("URA2",  "S1857"),
]


def atom_sites(x):
    return re.findall(
        r"[STY]\d+",
        str(x)
    )


def atom_parts(site):
    m = re.fullmatch(
        r"([STY])(\d+)",
        str(site)
    )
    if not m:
        raise ValueError(site)

    return (
        m.group(1),
        int(m.group(2)),
    )


                                                              
                       
                                                              

rows = []

for gene, site in DIRECT:
    aa, pos = atom_parts(site)

    rows.append({
        "gene": gene,
        "measurement_group": f"{gene}_{site}",
        "atomic_site": site,
        "aa": aa,
        "position": pos,
        "candidate_class": "direct_A",
    })


mem = pd.read_csv(
    NETWORKIN,
    sep="\t",
)

orth = mem[
    mem["analysis_label"] == ORTH_LABEL
].copy()

orth = orth.drop_duplicates(
    [
        "gene",
        "site_id",
    ]
)

if len(orth) != 20:
    raise RuntimeError(
        f"Expected 20 orthogonal groups, got {len(orth)}"
    )


for _, r in orth.iterrows():

    gene = str(
        r["gene"]
    ).strip().upper()

    site_id = str(
        r["site_id"]
    ).strip()

                                                            
    if site_id.upper().startswith(
        gene + "_"
    ):
        group_sites = site_id[
            len(gene) + 1:
        ]
    else:
        group_sites = site_id

    atoms = atom_sites(
        group_sites
    )

    for site in atoms:

        aa, pos = atom_parts(
            site
        )

        rows.append({
            "gene": gene,
            "measurement_group": f"{gene}_{group_sites}",
            "atomic_site": site,
            "aa": aa,
            "position": pos,
            "candidate_class": "orthogonal_NetworKIN",
        })


cand = pd.DataFrame(
    rows
).drop_duplicates()


                                                              
                                                    
                                                              

atlas_frames = []

for path, klass in [
    (
        DIRECT_ATLAS,
        "direct_A",
    ),
    (
        ORTH_ATLAS,
        "orthogonal_NetworKIN",
    ),
]:

    d = pd.read_csv(
        path,
        sep="\t",
    )

    for _, r in d.iterrows():

        gene = str(
            r["gene"]
        ).strip().upper()

        psites = str(
            r["pSites"]
        ).strip()

        atlas_frames.append({
            "gene":
                gene,

            "measurement_group":
                f"{gene}_{psites}",

            "Kanshin_delta5":
                pd.to_numeric(
                    r.get(
                        "Kanshin_heat_minus_cold_5min"
                    ),
                    errors="coerce",
                ),

            "PhosphoAtlas_HS42":
                pd.to_numeric(
                    r.get(
                        "PhosphoAtlas_HS42"
                    ),
                    errors="coerce",
                ),

            "PhosphoAtlas_CS18":
                pd.to_numeric(
                    r.get(
                        "PhosphoAtlas_CS18"
                    ),
                    errors="coerce",
                ),

            "PhosphoAtlas_HS42_minus_CS18":
                pd.to_numeric(
                    r.get(
                        "PhosphoAtlas_heat_minus_cold"
                    ),
                    errors="coerce",
                ),

            "PhosphoAtlas_HS42_q":
                pd.to_numeric(
                    r.get(
                        "HS42_q_min"
                    ),
                    errors="coerce",
                ),

            "PhosphoAtlas_CS18_q":
                pd.to_numeric(
                    r.get(
                        "CS18_q_min"
                    ),
                    errors="coerce",
                ),
        })


atlas = pd.DataFrame(
    atlas_frames
)

cand = cand.merge(
    atlas,
    on=[
        "gene",
        "measurement_group",
    ],
    how="left",
    validate="m:1",
)


                                                              
               
                                                 
                                                              

u = pd.read_csv(
    ULIANA,
    sep=";",
)

u["gene"] = (
    u[
        "Gene.names...primary.."
    ]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)

u["position"] = pd.to_numeric(
    u["position"],
    errors="coerce",
)

                                               
u["aa"] = (
    u["ID"]
    .astype(str)
    .str.extract(
        r"_([STY])_\d+$",
        expand=False,
    )
)

u = u[
    u["position"].notna()
    &
    u["aa"].notna()
].copy()

u["position"] = (
    u["position"]
    .astype(int)
)


ucols = [
    "gene",
    "aa",
    "position",
    "diff.day8_plus.day8_minus",
    "diff.Expo_plus.Expo_minus",
    "diff.HS_plus.HS_minus",
    "p_adj_BH#p_adj.day8_plus.day8_minus",
    "p_adj_BH#p_adj.Expo_plus.Expo_minus",
    "p_adj_BH#p_adj.HS_plus.HS_minus",
]


u = u[
    ucols
].rename(
    columns={
        "diff.day8_plus.day8_minus":
            "Uliana_day8_PKAi_minus_ctrl",

        "diff.Expo_plus.Expo_minus":
            "Uliana_Expo_PKAi_minus_ctrl",

        "diff.HS_plus.HS_minus":
            "Uliana_HS_PKAi_minus_ctrl",

        "p_adj_BH#p_adj.day8_plus.day8_minus":
            "Uliana_day8_q",

        "p_adj_BH#p_adj.Expo_plus.Expo_minus":
            "Uliana_Expo_q",

        "p_adj_BH#p_adj.HS_plus.HS_minus":
            "Uliana_HS_q",
    }
)


                                                       
                                                       
u = (
    u.groupby(
        [
            "gene",
            "aa",
            "position",
        ],
        as_index=False,
    )
    .agg(
        Uliana_day8_PKAi_minus_ctrl=(
            "Uliana_day8_PKAi_minus_ctrl",
            "median",
        ),

        Uliana_Expo_PKAi_minus_ctrl=(
            "Uliana_Expo_PKAi_minus_ctrl",
            "median",
        ),

        Uliana_HS_PKAi_minus_ctrl=(
            "Uliana_HS_PKAi_minus_ctrl",
            "median",
        ),

        Uliana_day8_q=(
            "Uliana_day8_q",
            "min",
        ),

        Uliana_Expo_q=(
            "Uliana_Expo_q",
            "min",
        ),

        Uliana_HS_q=(
            "Uliana_HS_q",
            "min",
        ),
    )
)


cand = cand.merge(
    u,
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="left",
    validate="m:1",
)


                                                              
                                                    
 
                           
                                  
                                        
                                                              

s = pd.read_csv(
    STARV,
    sep="\t",
    low_memory=False,
)


def genes_from_fasta(x):
    if pd.isna(x):
        return []

    return [
        z.upper()
        for z in re.findall(
            r"\bGN=([A-Za-z0-9_.-]+)",
            str(x)
        )
    ]


starv_rows = []

intensity_cols = [
    "Intensity D1_D5_1",
    "Intensity D1_D5_2",
    "Intensity Expo_1",
    "Intensity Expo_2",
]


for c in intensity_cols:
    if c not in s.columns:
        raise RuntimeError(
            f"Missing starvation column: {c}"
        )


for _, r in s.iterrows():

    genes = genes_from_fasta(
        r.get(
            "Fasta headers",
            ""
        )
    )

    if not genes:
        continue

    try:
        pos = int(
            float(
                r["Position"]
            )
        )
    except Exception:
        continue

    aa = str(
        r.get(
            "Amino acid",
            ""
        )
    ).strip().upper()

    if aa not in {
        "S",
        "T",
        "Y",
    }:
        continue

    vals = {}

    for c in intensity_cols:

        v = pd.to_numeric(
            r.get(c),
            errors="coerce",
        )

                                                 
                                       
        if pd.notna(v) and v > 0:
            vals[c] = float(
                np.log2(v)
            )
        else:
            vals[c] = np.nan


    starvation_values = np.array(
        [
            vals[
                "Intensity D1_D5_1"
            ],
            vals[
                "Intensity D1_D5_2"
            ],
        ],
        dtype=float,
    )

    expo_values = np.array(
        [
            vals[
                "Intensity Expo_1"
            ],
            vals[
                "Intensity Expo_2"
            ],
        ],
        dtype=float,
    )


    n_starv = int(
        np.isfinite(
            starvation_values
        ).sum()
    )

    n_expo = int(
        np.isfinite(
            expo_values
        ).sum()
    )


    starv_mean = (
        float(
            np.nanmean(
                starvation_values
            )
        )
        if n_starv
        else np.nan
    )

    expo_mean = (
        float(
            np.nanmean(
                expo_values
            )
        )
        if n_expo
        else np.nan
    )


    localization = pd.to_numeric(
        r.get(
            "Localization prob"
        ),
        errors="coerce",
    )


    for gene in genes:

        starv_rows.append({
            "gene":
                gene,

            "aa":
                aa,

            "position":
                pos,

            "starvation_localization_prob":
                localization,

            "starvation_n":
                n_starv,

            "expo_n":
                n_expo,

            "starvation_mean_log2_intensity":
                starv_mean,

            "expo_mean_log2_intensity":
                expo_mean,

            "starvation_minus_expo_log2":
                (
                    starv_mean
                    - expo_mean
                    if (
                        np.isfinite(
                            starv_mean
                        )
                        and
                        np.isfinite(
                            expo_mean
                        )
                    )
                    else np.nan
                ),
        })


st = pd.DataFrame(
    starv_rows
)


                                                      
                                  
st = (
    st.sort_values(
        "starvation_localization_prob",
        ascending=False,
    )
    .drop_duplicates(
        [
            "gene",
            "aa",
            "position",
        ],
        keep="first",
    )
)


cand = cand.merge(
    st,
    on=[
        "gene",
        "aa",
        "position",
    ],
    how="left",
    validate="m:1",
)


                                                              
                                                           
                                                              

iv = pd.read_csv(
    CDC19_IV,
    sep=";",
)


iv22 = iv[
    iv["ID"] == "CDC19_22"
].copy()


if len(iv22) != 1:
    raise RuntimeError(
        f"Expected one CDC19_22 row; found {len(iv22)}"
    )


r = iv22.iloc[0]


def mean_numeric(cols):

    v = pd.to_numeric(
        r[
            cols
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    good = np.isfinite(
        v
    )

    return (
        float(
            np.mean(
                v[
                    good
                ]
            )
        )
        if good.any()
        else np.nan,
        int(
            good.sum()
        ),
    )


ctrl_1h, n_ctrl_1h = mean_numeric(
    [
        "CTRL_1hr_1",
        "CTRL_1hr_2",
        "CTRL_1hr_3",
    ]
)

pka_1h, n_pka_1h = mean_numeric(
    [
        "PKA_1hr_1",
        "PKA_1hr_2",
        "PKA_1hr_3",
    ]
)

ctrl_on, n_ctrl_on = mean_numeric(
    [
        "CTRL_on_1",
        "CTRL_on_2",
        "CTRL_on_3",
    ]
)

pka_on, n_pka_on = mean_numeric(
    [
        "PKA_on_1",
        "PKA_on_2",
        "PKA_on_3",
    ]
)


cand[
    "CDC19_invitro_PKA_direct"
] = False

cand[
    "CDC19_invitro_ctrl_1h_mean_log2"
] = np.nan

cand[
    "CDC19_invitro_PKA_1h_mean_log2"
] = np.nan

cand[
    "CDC19_invitro_PKA_minus_ctrl_1h"
] = np.nan

cand[
    "CDC19_invitro_ctrl_overnight_mean_log2"
] = np.nan

cand[
    "CDC19_invitro_PKA_overnight_mean_log2"
] = np.nan

cand[
    "CDC19_invitro_PKA_minus_ctrl_overnight"
] = np.nan


mask = (
    (cand["gene"] == "CDC19")
    &
    (cand["atomic_site"] == "S22")
)


cand.loc[
    mask,
    "CDC19_invitro_PKA_direct",
] = True

cand.loc[
    mask,
    "CDC19_invitro_ctrl_1h_mean_log2",
] = ctrl_1h

cand.loc[
    mask,
    "CDC19_invitro_PKA_1h_mean_log2",
] = pka_1h

cand.loc[
    mask,
    "CDC19_invitro_PKA_minus_ctrl_1h",
] = (
    pka_1h
    - ctrl_1h
)

cand.loc[
    mask,
    "CDC19_invitro_ctrl_overnight_mean_log2",
] = ctrl_on

cand.loc[
    mask,
    "CDC19_invitro_PKA_overnight_mean_log2",
] = pka_on

cand.loc[
    mask,
    "CDC19_invitro_PKA_minus_ctrl_overnight",
] = (
    pka_on
    - ctrl_on
)


                                                              
                                                  
                                                              

cand[
    "Atlas_exact_match"
] = cand[
    "PhosphoAtlas_HS42_minus_CS18"
].notna()

cand[
    "Uliana_exact_match"
] = (
    cand[
        "Uliana_Expo_PKAi_minus_ctrl"
    ].notna()
    |
    cand[
        "Uliana_HS_PKAi_minus_ctrl"
    ].notna()
    |
    cand[
        "Uliana_day8_PKAi_minus_ctrl"
    ].notna()
)

cand[
    "starvation_exact_match"
] = cand[
    "starvation_minus_expo_log2"
].notna()


                                                                     
cand[
    "n_external_evidence_layers"
] = (
    cand[
        [
            "Atlas_exact_match",
            "Uliana_exact_match",
            "starvation_exact_match",
            "CDC19_invitro_PKA_direct",
        ]
    ]
    .fillna(False)
    .astype(bool)
    .sum(axis=1)
)


cand = cand.sort_values(
    [
        "n_external_evidence_layers",
        "candidate_class",
        "gene",
        "position",
    ],
    ascending=[
        False,
        True,
        True,
        True,
    ],
)


cand.to_csv(
    OUT
    / "PKA_SITE_CONTEXT_MATRIX.tsv",
    sep="\t",
    index=False,
)


                                                              
                                     
                                                              

keep = [
    "gene",
    "atomic_site",
    "measurement_group",
    "candidate_class",

    "Kanshin_delta5",
    "PhosphoAtlas_HS42_minus_CS18",
    "PhosphoAtlas_HS42_q",

    "Uliana_Expo_PKAi_minus_ctrl",
    "Uliana_Expo_q",

    "Uliana_HS_PKAi_minus_ctrl",
    "Uliana_HS_q",

    "Uliana_day8_PKAi_minus_ctrl",
    "Uliana_day8_q",

    "starvation_localization_prob",
    "starvation_minus_expo_log2",

    "CDC19_invitro_PKA_direct",
    "CDC19_invitro_PKA_minus_ctrl_1h",
    "CDC19_invitro_PKA_minus_ctrl_overnight",

    "n_external_evidence_layers",
]


cand[
    keep
].to_csv(
    OUT
    / "STRUCTURAL_FOLLOWUP_CONTEXT.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 140)
print("PKA SITE CONTEXT MATRIX")
print("=" * 140)

print(
    cand[
        keep
    ].to_string(
        index=False
    )
)

print()
print("Candidate atomic sites:", len(cand))
print(
    "Atlas matches:",
    int(
        cand[
            "Atlas_exact_match"
        ].sum()
    )
)
print(
    "Uliana PKAi exact matches:",
    int(
        cand[
            "Uliana_exact_match"
        ].sum()
    )
)
print(
    "Starvation exact matches:",
    int(
        cand[
            "starvation_exact_match"
        ].sum()
    )
)

print()
print("CDC19 S22 in-vitro PKA:")
print(
    "  ctrl 1h mean:",
    ctrl_1h,
    "n=",
    n_ctrl_1h,
)
print(
    "  PKA 1h mean:",
    pka_1h,
    "n=",
    n_pka_1h,
)
print(
    "  delta:",
    pka_1h - ctrl_1h,
)

print(
    "  ctrl overnight mean:",
    ctrl_on,
    "n=",
    n_ctrl_on,
)
print(
    "  PKA overnight mean:",
    pka_on,
    "n=",
    n_pka_on,
)
print(
    "  delta:",
    pka_on - ctrl_on,
)

print()
print("Output:", OUT)
