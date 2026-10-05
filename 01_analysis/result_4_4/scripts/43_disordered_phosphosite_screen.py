#!/usr/bin/env python3

from pathlib import Path
import json
import re
import numpy as np
import pandas as pd


ROOT = Path(
    "/"
)

OUT = (
    ROOT
    / "42_disordered_phosphosite_mechanisms"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


PHASE1 = (
    ROOT
    / "40_structural_tractability"
    / "STRUCTURAL_TRACTABILITY_PHASE1.tsv"
)

CONTEXT = (
    ROOT
    / "39_mechanistic_site_context"
    / "PKA_SITE_CONTEXT_MATRIX_COMPLETE.tsv"
)

COVERAGE = (
    ROOT
    / "40_structural_tractability"
    / "EXACT_EXPERIMENTAL_RESIDUE_COVERAGE.tsv"
)

UNIPROT_DIR = (
    ROOT
    / "40_structural_tractability"
    / "raw"
    / "uniprot"
)


TARGETS = [
    ("ATG1",  "S515"),
    ("MAF1",  "S90"),
    ("MSN4",  "S316"),
    ("TSL1",  "S77"),
    ("NUP60", "S10"),
]


                                                              
         
                                                              

def site_parts(site):

    m = re.fullmatch(
        r"([STY])(\d+)",
        site,
    )

    if not m:
        raise ValueError(site)

    return (
        m.group(1),
        int(m.group(2)),
    )


def get_bounds(feature):

    loc = feature.get(
        "location",
        {}
    )

    try:
        return (
            int(
                loc["start"]["value"]
            ),
            int(
                loc["end"]["value"]
            ),
        )
    except Exception:
        return (
            None,
            None,
        )


def feature_text(feature):

    start, end = get_bounds(
        feature
    )

    return (
        f"{feature.get('type','')}:"
        f"{feature.get('description','')}"
        f"[{start}-{end}]"
    )


def clean_bool(x):

    if isinstance(x, bool):
        return x

    return str(x).lower() in {
        "true",
        "1",
        "yes",
    }


                                                              
        
                                                              

p1 = pd.read_csv(
    PHASE1,
    sep="\t",
)

ctx = pd.read_csv(
    CONTEXT,
    sep="\t",
)

if COVERAGE.exists():

    cov = pd.read_csv(
        COVERAGE,
        sep="\t",
    )

else:

    cov = pd.DataFrame()


                                                              
                                       
                                                              

coord_summary = {}


if len(cov):

    for (
        gene,
        site
    ), d in cov.groupby(
        [
            "gene",
            "site",
        ]
    ):

        mapped = (
            d[
                "site_in_entity_alignment"
            ]
            .fillna(False)
            .astype(bool)
        )

        coords = (
            d[
                "coordinate_present"
            ]
            .fillna(False)
            .astype(bool)
        )

        coord_summary[
            (
                str(gene).upper(),
                str(site),
            )
        ] = {
            "n_experimental_entities":
                len(d),

            "n_sequence_mapped_entities":
                int(
                    mapped.sum()
                ),

            "n_coordinate_resolved_entities":
                int(
                    coords.sum()
                ),
        }


                                                              
                                              
                                                              

rows = []


for gene, site in TARGETS:

    aa, pos = site_parts(
        site
    )


    q = p1[
        (p1["gene"] == gene)
        &
        (p1["site"] == site)
    ]


    if len(q) != 1:
        raise RuntimeError(
            f"{gene} {site}: "
            f"expected one phase1 row, got {len(q)}"
        )


    q = q.iloc[0]

    accession = str(
        q[
            "uniprot_accession"
        ]
    )


                                                              
                              
                                                              

    hits = list(
        UNIPROT_DIR.glob(
            f"{gene}_{site}_{accession}.json"
        )
    )


    if len(hits) != 1:
        raise RuntimeError(
            f"{gene} {site}: "
            f"expected one UniProt JSON, got {hits}"
        )


    uj = json.loads(
        hits[0].read_text()
    )


    seq = (
        uj.get(
            "sequence",
            {}
        )
        .get(
            "value",
            ""
        )
    )


    if not (
        1 <= pos <= len(seq)
    ):
        raise RuntimeError(
            f"{gene} {site}: position out of range"
        )


    observed_aa = seq[
        pos - 1
    ]


    if observed_aa != aa:
        raise RuntimeError(
            f"{gene} {site}: "
            f"expected {aa}, sequence has {observed_aa}"
        )


                                                              
                            
                                                              

    lo20 = max(
        1,
        pos - 20,
    )

    hi20 = min(
        len(seq),
        pos + 20,
    )

    seq41 = seq[
        lo20 - 1:
        hi20
    ]


    lo7 = max(
        1,
        pos - 7,
    )

    hi7 = min(
        len(seq),
        pos + 7,
    )

    seq15 = seq[
        lo7 - 1:
        hi7
    ]


                                        
    def res(offset):

        p = pos + offset

        if (
            p < 1
            or p > len(seq)
        ):
            return ""

        return seq[
            p - 1
        ]


    minus4 = res(-4)
    minus3 = res(-3)
    minus2 = res(-2)
    minus1 = res(-1)

    plus1 = res(1)
    plus2 = res(2)


                                                              
                                            
     
                                                          
                              
                                                              

    basic = (
        seq41.count("K")
        +
        seq41.count("R")
    )

    acidic = (
        seq41.count("D")
        +
        seq41.count("E")
    )

    histidine = seq41.count(
        "H"
    )


    formal_charge_proxy = (
        basic
        - acidic
    )


                                                         
                                             
    phospho_charge_proxy = (
        formal_charge_proxy
        - 2
    )


                                                              
                                      
     
                                                 
                                                              

    upstream_basic_count = sum(
        z in {
            "K",
            "R",
        }
        for z in [
            minus4,
            minus3,
            minus2,
            minus1,
        ]
    )


    basic_at_minus3 = (
        minus3 in {
            "K",
            "R",
        }
    )

    basic_at_minus2 = (
        minus2 in {
            "K",
            "R",
        }
    )

    proline_plus1 = (
        plus1 == "P"
    )


                                                              
                             
                                                              

    at_site = []
    near20 = []

    disorder_at_site = []
    mutagenesis_near = []
    binding_near = []
    modified_at_site = []


    for f in uj.get(
        "features",
        []
    ):

        start, end = get_bounds(
            f
        )

        if (
            start is None
            or end is None
        ):
            continue


        typ = str(
            f.get(
                "type",
                ""
            )
        )

        txt = feature_text(
            f
        )


        overlaps_site = (
            start
            <= pos
            <= end
        )


        overlaps_window = not (
            end < lo20
            or
            start > hi20
        )


        if overlaps_site:

            at_site.append(
                txt
            )


            if (
                typ.lower()
                == "region"
                and
                "disorder"
                in str(
                    f.get(
                        "description",
                        ""
                    )
                ).lower()
            ):
                disorder_at_site.append(
                    txt
                )


            if typ.lower() == "modified residue":

                modified_at_site.append(
                    txt
                )


        if overlaps_window:

            near20.append(
                txt
            )


            if typ.lower() == "mutagenesis":

                mutagenesis_near.append(
                    txt
                )


            if typ.lower() in {
                "binding site",
                "region",
                "motif",
                "short sequence motif",
                "cross-link",
            }:

                binding_near.append(
                    txt
                )


                                                              
                                   
                                                              

    cs = coord_summary.get(
        (
            gene,
            site,
        ),
        {
            "n_experimental_entities":
                0,

            "n_sequence_mapped_entities":
                0,

            "n_coordinate_resolved_entities":
                0,
        }
    )


                                                              
                         
                                                              

    cx = ctx[
        (ctx["gene"] == gene)
        &
        (ctx["atomic_site"] == site)
    ]


    if len(cx) != 1:

        raise RuntimeError(
            f"{gene} {site}: "
            f"expected one mechanistic-context row, "
            f"got {len(cx)}"
        )


    cx = cx.iloc[0]


                                                              
                                     
                                                              

    plddt = pd.to_numeric(
        q.get(
            "alphafold_site_plddt"
        ),
        errors="coerce",
    )


    explicit_disorder = bool(
        disorder_at_site
    )


    unresolved_exp = (
        cs[
            "n_sequence_mapped_entities"
        ]
        > 0
        and
        cs[
            "n_coordinate_resolved_entities"
        ]
        == 0
    )


    if (
        explicit_disorder
        or
        (
            pd.notna(plddt)
            and plddt < 50
        )
    ):

        interpretation_class = (
            "IDR_or_highly_flexible_regulatory_site"
        )

    elif unresolved_exp:

        interpretation_class = (
            "experimentally_unresolved_flexible_site"
        )

    elif (
        pd.notna(plddt)
        and plddt < 70
    ):

        interpretation_class = (
            "low_confidence_flexible_site"
        )

    else:

        interpretation_class = (
            "structured_candidate"
        )


    row = {
        "gene":
            gene,

        "site":
            site,

        "uniprot":
            accession,

        "position":
            pos,

        "sequence_length":
            len(seq),

        "local_sequence_15aa":
            seq15,

        "local_sequence_41aa":
            seq41,

        "residue_minus4":
            minus4,

        "residue_minus3":
            minus3,

        "residue_minus2":
            minus2,

        "residue_minus1":
            minus1,

        "residue_plus1":
            plus1,

        "residue_plus2":
            plus2,

        "upstream_basic_count_minus4_to_minus1":
            upstream_basic_count,

        "basic_at_minus3":
            basic_at_minus3,

        "basic_at_minus2":
            basic_at_minus2,

        "proline_at_plus1":
            proline_plus1,

        "local41_basic_KR":
            basic,

        "local41_acidic_DE":
            acidic,

        "local41_histidine":
            histidine,

        "local41_charge_proxy_unmodified":
            formal_charge_proxy,

        "local41_charge_proxy_phosphorylated":
            phospho_charge_proxy,

        "local41_charge_shift_from_phosphate":
            -2,

        "alphafold_site_plddt":
            plddt,

        "alphafold_local_21aa_median_plddt":
            pd.to_numeric(
                q.get(
                    "alphafold_local_21aa_median_plddt"
                ),
                errors="coerce",
            ),

        "uniprot_explicit_disorder_at_site":
            explicit_disorder,

        "uniprot_features_at_site":
            " | ".join(
                at_site
            ),

        "uniprot_modified_residue_at_site":
            " | ".join(
                modified_at_site
            ),

        "uniprot_features_within_20aa":
            " | ".join(
                near20
            ),

        "mutagenesis_within_20aa":
            " | ".join(
                mutagenesis_near
            ),

        "binding_or_region_features_within_20aa":
            " | ".join(
                binding_near
            ),

        **cs,

        "experimentally_sequence_mapped_but_unresolved":
            unresolved_exp,

        "structural_interpretation_class":
            interpretation_class,

        "Kanshin_delta5":
            cx.get(
                "Kanshin_delta5"
            ),

        "PhosphoAtlas_HS42_minus_CS18":
            cx.get(
                "PhosphoAtlas_HS42_minus_CS18"
            ),

        "Uliana_Expo_PKAi_minus_ctrl":
            cx.get(
                "Uliana_Expo_PKAi_minus_ctrl"
            ),

        "Uliana_Expo_q":
            cx.get(
                "Uliana_Expo_q"
            ),

        "Uliana_HS_PKAi_minus_ctrl":
            cx.get(
                "Uliana_HS_PKAi_minus_ctrl"
            ),

        "Uliana_HS_q":
            cx.get(
                "Uliana_HS_q"
            ),

        "starvation_minus_expo_log2":
            cx.get(
                "starvation_minus_expo_log2"
            ),

        "site_specific_context_switch":
            cx.get(
                "site_specific_context_switch",
                False,
            ),
    }


    rows.append(
        row
    )


out = pd.DataFrame(
    rows
)


out.to_csv(
    OUT
    / "DISORDERED_PHOSPHOSITE_MECHANISTIC_SCREEN.tsv",
    sep="\t",
    index=False,
)


                                                              
                             
                                                              

show = [
    "gene",
    "site",

    "local_sequence_15aa",

    "residue_minus3",
    "residue_minus2",
    "residue_minus1",
    "residue_plus1",

    "upstream_basic_count_minus4_to_minus1",

    "local41_basic_KR",
    "local41_acidic_DE",
    "local41_charge_proxy_unmodified",
    "local41_charge_proxy_phosphorylated",

    "alphafold_site_plddt",
    "alphafold_local_21aa_median_plddt",

    "uniprot_explicit_disorder_at_site",

    "n_experimental_entities",
    "n_sequence_mapped_entities",
    "n_coordinate_resolved_entities",

    "experimentally_sequence_mapped_but_unresolved",

    "uniprot_modified_residue_at_site",
    "mutagenesis_within_20aa",

    "Kanshin_delta5",
    "PhosphoAtlas_HS42_minus_CS18",

    "Uliana_Expo_PKAi_minus_ctrl",
    "Uliana_Expo_q",

    "Uliana_HS_PKAi_minus_ctrl",
    "Uliana_HS_q",

    "starvation_minus_expo_log2",

    "site_specific_context_switch",

    "structural_interpretation_class",
]


out[
    show
].to_csv(
    OUT
    / "DISORDERED_PHOSPHOSITE_COMPACT.tsv",
    sep="\t",
    index=False,
)


print()
print(
    "=" * 150
)
print(
    "DISORDERED / FLEXIBLE PHOSPHOSITE SCREEN"
)
print(
    "=" * 150
)


print(
    out[
        show
    ].to_string(
        index=False
    )
)


print()
print(
    "Output:",
    OUT
)
