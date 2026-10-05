#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(
    "/"
)

CTX = (
    ROOT
    / "39_mechanistic_site_context"
    / "PKA_SITE_CONTEXT_MATRIX.tsv"
)

DIRECT_CTX = (
    ROOT
    / "31_uliana_context_interaction"
    / "DIRECT_A_PKA_context_validation.tsv"
)

OUT = (
    ROOT
    / "39_mechanistic_site_context"
)

x = pd.read_csv(
    CTX,
    sep="\t",
)

d = pd.read_csv(
    DIRECT_CTX,
    sep="\t",
)

print("DIRECT CONTEXT COLUMNS:")
for c in d.columns:
    print(" ", c)


                                                              
                                 
                                                              

d["gene"] = (
    d["gene"]
    .astype(str)
    .str.strip()
    .str.upper()
)

d["atomic_site"] = (
    d["site"]
    .astype(str)
    .str.strip()
)


wanted = [
    "gene",
    "atomic_site",
]

rename = {}


candidate_fields = {
    "interaction_beta":
        "Uliana_context_interaction_beta",

    "interaction_q_directA_overlap":
        "Uliana_context_interaction_q",

    "targeted_Expo_PKAi_minus_ctrl":
        "Uliana_targeted_Expo_PKAi_minus_ctrl",

    "targeted_HS_PKAi_minus_ctrl":
        "Uliana_targeted_HS_PKAi_minus_ctrl",

    "targeted_interaction_beta":
        "Uliana_targeted_context_interaction_beta",

    "targeted_interaction_p":
        "Uliana_targeted_context_interaction_p",

    "observed_targeted_PRM":
        "Uliana_targeted_PRM_observed",
}


for old, new in candidate_fields.items():
    if old in d.columns:
        wanted.append(old)
        rename[old] = new


z = (
    d[wanted]
    .rename(columns=rename)
    .drop_duplicates(
        [
            "gene",
            "atomic_site",
        ]
    )
)


                                                              
                                 
                                                              

y = x.merge(
    z,
    on=[
        "gene",
        "atomic_site",
    ],
    how="left",
    validate="m:1",
)


                                                              
                            
                                                              

y[
    "significant_untargeted_context_interaction"
] = (
    pd.to_numeric(
        y.get(
            "Uliana_context_interaction_q",
            np.nan,
        ),
        errors="coerce",
    )
    < 0.05
)


y[
    "significant_targeted_context_interaction"
] = (
    pd.to_numeric(
        y.get(
            "Uliana_targeted_context_interaction_p",
            np.nan,
        ),
        errors="coerce",
    )
    < 0.05
)


y[
    "site_specific_context_switch"
] = (
    y[
        "significant_untargeted_context_interaction"
    ].fillna(False)
    |
    y[
        "significant_targeted_context_interaction"
    ].fillna(False)
)


                                       
                                               
y[
    "n_mechanistic_evidence_layers"
] = (
    pd.to_numeric(
        y[
            "n_external_evidence_layers"
        ],
        errors="coerce",
    )
    .fillna(0)
    .astype(int)
    +
    y[
        "site_specific_context_switch"
    ]
    .fillna(False)
    .astype(int)
)


y = y.sort_values(
    [
        "n_mechanistic_evidence_layers",
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


y.to_csv(
    OUT
    / "PKA_SITE_CONTEXT_MATRIX_COMPLETE.tsv",
    sep="\t",
    index=False,
)


show = [
    "gene",
    "atomic_site",
    "measurement_group",
    "candidate_class",

    "Kanshin_delta5",
    "PhosphoAtlas_HS42_minus_CS18",

    "Uliana_Expo_PKAi_minus_ctrl",
    "Uliana_Expo_q",

    "Uliana_HS_PKAi_minus_ctrl",
    "Uliana_HS_q",

    "Uliana_day8_PKAi_minus_ctrl",
    "Uliana_day8_q",

    "starvation_minus_expo_log2",

    "CDC19_invitro_PKA_direct",
    "CDC19_invitro_PKA_minus_ctrl_1h",

    "Uliana_context_interaction_beta",
    "Uliana_context_interaction_q",

    "Uliana_targeted_Expo_PKAi_minus_ctrl",
    "Uliana_targeted_HS_PKAi_minus_ctrl",
    "Uliana_targeted_context_interaction_beta",
    "Uliana_targeted_context_interaction_p",

    "site_specific_context_switch",
    "n_mechanistic_evidence_layers",
]

show = [
    c for c in show
    if c in y.columns
]


y[
    show
].to_csv(
    OUT
    / "MECHANISTIC_CONTEXT_COMPLETE_COMPACT.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 150)
print("COMPLETE MECHANISTIC CONTEXT")
print("=" * 150)

print(
    y[
        show
    ].to_string(
        index=False
    )
)

print()
print(
    "context-switching sites:"
)

print(
    y.loc[
        y[
            "site_specific_context_switch"
        ],
        [
            "gene",
            "atomic_site",
        ],
    ].to_string(
        index=False
    )
)

