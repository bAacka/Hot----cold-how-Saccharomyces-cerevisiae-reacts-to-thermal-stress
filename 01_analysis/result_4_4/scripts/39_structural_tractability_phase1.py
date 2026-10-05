#!/usr/bin/env python3

from pathlib import Path
import json
import re
import time

import numpy as np
import pandas as pd
import requests


ROOT = Path(
    "/"
)

OUT = ROOT / "40_structural_tractability"

RAW_UNIPROT = OUT / "raw" / "uniprot"
RAW_AF = OUT / "raw" / "alphafold"
RAW_RCSB = OUT / "raw" / "rcsb"

for p in [
    OUT,
    RAW_UNIPROT,
    RAW_AF,
    RAW_RCSB,
]:
    p.mkdir(
        parents=True,
        exist_ok=True,
    )


S1FILE = (
    ROOT
    / "00_raw"
    / "Kanshin_2015_Table_S1.xlsx"
)


                                                              
                
 
                                                      
                                              
                                                              

CANDIDATES = [
    {
        "gene": "CDC19",
        "site": "S22",
        "archetype":
            "direct_invitro_PKA_context_specific",
    },
    {
        "gene": "ATG1",
        "site": "S515",
        "archetype":
            "targeted_context_switch",
    },
    {
        "gene": "MAF1",
        "site": "S90",
        "archetype":
            "untargeted_context_switch",
    },
    {
        "gene": "MSN4",
        "site": "S316",
        "archetype":
            "strong_replicated_acute",
    },
    {
        "gene": "TSL1",
        "site": "S77",
        "archetype":
            "opposite_direction_replicated",
    },
    {
        "gene": "NUP60",
        "site": "S10",
        "archetype":
            "strong_thermal_replication",
    },
]


SESSION = requests.Session()

SESSION.headers.update({
    "User-Agent":
        "stefan-yeast-structural-screen/1.0"
})


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def site_parts(site):

    m = re.fullmatch(
        r"([STY])(\d+)",
        site,
    )

    if not m:
        raise ValueError(
            site
        )

    return (
        m.group(1),
        int(
            m.group(2)
        ),
    )


def safe_json_get(
    url,
    timeout=60,
):

    r = SESSION.get(
        url,
        timeout=timeout,
    )

    if r.status_code == 404:
        return None

    r.raise_for_status()

    return r.json()


                                                              
                                                       
                                                              

s1 = pd.read_excel(
    S1FILE,
    sheet_name="Phosphopeptides",
)


def norm_header(x):
    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(x).lower(),
    )


def resolve_column(columns, aliases, label):
    by_norm = {
        norm_header(c): c
        for c in columns
    }

                                      
    for alias in aliases:
        key = norm_header(alias)

        if key in by_norm:
            return by_norm[key]

                            
    candidates = [
        c
        for c in columns
        if any(
            norm_header(alias)
            in norm_header(c)
            for alias in aliases
        )
    ]

    candidates = list(
        dict.fromkeys(candidates)
    )

    if len(candidates) == 1:
        return candidates[0]

    raise RuntimeError(
        f"Could not uniquely resolve {label}. "
        f"Candidates={candidates}; "
        f"columns={list(columns)}"
    )


GENE_COL = resolve_column(
    s1.columns,
    [
        "Gene",
    ],
    "gene column",
)

UNIPROT_COL = resolve_column(
    s1.columns,
    [
        "UniProt",
        "Uniprot",
        "UniProt ID",
        "UniProt accession",
        "Uniprot accession",
    ],
    "UniProt column",
)


print(
    "Resolved S1 columns:"
)
print(
    "  gene    =",
    repr(GENE_COL),
)
print(
    "  UniProt =",
    repr(UNIPROT_COL),
)


s1["Gene_norm"] = (
    s1[GENE_COL]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)


def extract_accessions(series):

    out = []

    for raw in series.dropna():

        text = str(
            raw
        ).strip()

        if not text:
            continue

                                                
        pipe = re.findall(
            r"(?:sp|tr)\|([^|]+)\|",
            text,
        )

        if pipe:
            out.extend(
                pipe
            )
            continue

                                                     
        for x in re.split(
            r"[;,\s]+",
            text,
        ):

            x = x.strip()

            if x:
                out.append(
                    x
                )

    return sorted(
        set(
            out
        )
    )


gene_to_accs = {}

for c in CANDIDATES:

    gene = c[
        "gene"
    ]

    vals = s1.loc[
        s1["Gene_norm"] == gene,
        UNIPROT_COL,
    ]

    gene_to_accs[
        gene
    ] = extract_accessions(
        vals
    )


print(
    "===== KANSHIN GENE -> UNIPROT ====="
)

for gene, accs in gene_to_accs.items():

    print(
        gene,
        "=>",
        accs,
    )


                                                              
                              
                                                              

def feature_bounds(feature):

    loc = feature.get(
        "location",
        {},
    )

    try:
        start = int(
            loc[
                "start"
            ][
                "value"
            ]
        )

        end = int(
            loc[
                "end"
            ][
                "value"
            ]
        )

        return (
            start,
            end,
        )

    except Exception:

        return (
            None,
            None,
        )


def feature_string(feature):

    start, end = feature_bounds(
        feature
    )

    typ = feature.get(
        "type",
        ""
    )

    desc = feature.get(
        "description",
        ""
    )

    return (
        f"{typ}:{desc}"
        f"[{start}-{end}]"
    )


                                                              
                              
                                                              

def parse_af_ca(
    path,
):

    rows = []

    with open(
        path,
        errors="replace",
    ) as fh:

        for line in fh:

            if not line.startswith(
                "ATOM"
            ):
                continue

            atom = line[
                12:16
            ].strip()

            if atom != "CA":
                continue

            resname = line[
                17:20
            ].strip()

            chain = line[
                21:22
            ].strip()

            try:
                resnum = int(
                    line[
                        22:26
                    ]
                )

                bfactor = float(
                    line[
                        60:66
                    ]
                )

                x = float(
                    line[
                        30:38
                    ]
                )

                y = float(
                    line[
                        38:46
                    ]
                )

                z = float(
                    line[
                        46:54
                    ]
                )

            except Exception:
                continue

            rows.append({
                "chain":
                    chain,

                "resnum":
                    resnum,

                "resname3":
                    resname,

                "residue":
                    AA3_TO_1.get(
                        resname,
                        "?"
                    ),

                "plddt":
                    bfactor,

                "x":
                    x,

                "y":
                    y,

                "z":
                    z,
            })

    return pd.DataFrame(
        rows
    )


                                                              
                                       
                                                              

def rcsb_uniprot_search(
    accession,
):

    url = (
        "https://search.rcsb.org/"
        "rcsbsearch/v2/query"
    )

    payload = {
        "query": {
            "type":
                "terminal",

            "service":
                "text",

            "parameters": {
                "attribute":
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers."
                    "database_accession",

                "operator":
                    "exact_match",

                "value":
                    accession,
            },
        },

        "return_type":
            "polymer_entity",

        "request_options": {
            "return_all_hits":
                True,

            "results_content_type": [
                "experimental"
            ],
        },
    }

    r = SESSION.post(
        url,
        json=payload,
        timeout=60,
    )

                              
    if r.status_code == 204:
        return {
            "total_count": 0,
            "result_set": [],
        }

    r.raise_for_status()

    return r.json()


                                                              
                              
                                                              

rows = []


for c in CANDIDATES:

    gene = c[
        "gene"
    ]

    site = c[
        "site"
    ]

    expected_aa, pos = site_parts(
        site
    )

    accs = gene_to_accs.get(
        gene,
        [],
    )


    print()
    print(
        "=" * 90
    )
    print(
        gene,
        site,
        c[
            "archetype"
        ],
    )
    print(
        "=" * 90
    )


    rec = {
        "gene":
            gene,

        "site":
            site,

        "position":
            pos,

        "expected_residue":
            expected_aa,

        "archetype":
            c[
                "archetype"
            ],

        "kanshin_uniprot_candidates":
            ";".join(
                accs
            ),

        "uniprot_mapping_status":
            (
                "unique"
                if len(
                    accs
                ) == 1
                else
                (
                    "missing"
                    if len(
                        accs
                    ) == 0
                    else
                    "ambiguous"
                )
            ),
    }


    if len(
        accs
    ) != 1:

        rows.append(
            rec
        )

        print(
            "Skipping web structural lookup:",
            len(
                accs
            ),
            "UniProt candidates",
        )

        continue


    acc = accs[
        0
    ]

    rec[
        "uniprot_accession"
    ] = acc


                                                              
                                              
                                                              

    uj = safe_json_get(
        f"https://rest.uniprot.org/"
        f"uniprotkb/{acc}.json"
    )


    if uj is None:

        rec[
            "uniprot_fetch"
        ] = "404"

        rows.append(
            rec
        )

        continue


    with open(
        RAW_UNIPROT
        / f"{gene}_{site}_{acc}.json",
        "w",
    ) as fh:

        json.dump(
            uj,
            fh,
            indent=2,
        )


    rec[
        "uniprot_fetch"
    ] = "ok"


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


    rec[
        "sequence_length"
    ] = len(
        seq
    )


    if 1 <= pos <= len(
        seq
    ):

        seq_res = seq[
            pos - 1
        ]

    else:

        seq_res = ""


    rec[
        "uniprot_site_residue"
    ] = seq_res


    rec[
        "sequence_site_matches"
    ] = (
        seq_res
        == expected_aa
    )


    features = uj.get(
        "features",
        []
    )


    overlap = []
    nearby = []
    modified_here = []


    for f in features:

        start, end = feature_bounds(
            f
        )

        if (
            start is None
            or end is None
        ):
            continue


        if (
            start
            <= pos
            <= end
        ):

            overlap.append(
                feature_string(
                    f
                )
            )

            if f.get(
                "type"
            ) == "Modified residue":

                modified_here.append(
                    feature_string(
                        f
                    )
                )


        elif (
            abs(
                pos
                - start
            )
            <= 20
            or
            abs(
                pos
                - end
            )
            <= 20
        ):

            nearby.append(
                feature_string(
                    f
                )
            )


    rec[
        "uniprot_features_at_site"
    ] = " | ".join(
        overlap
    )


    rec[
        "uniprot_features_within_20aa"
    ] = " | ".join(
        nearby
    )


    rec[
        "uniprot_modified_residue_annotation"
    ] = " | ".join(
        modified_here
    )


                                                              
               
                                                              

    af_url = (
        "https://alphafold.ebi.ac.uk/"
        f"api/prediction/{acc}"
    )

    af = safe_json_get(
        af_url
    )


    if af is None:

        rec[
            "alphafold_available"
        ] = False

    else:

        if isinstance(
            af,
            dict,
        ):
            af = [
                af
            ]


        with open(
            RAW_AF
            / f"{gene}_{site}_{acc}_metadata.json",
            "w",
        ) as fh:

            json.dump(
                af,
                fh,
                indent=2,
            )


        covering = []

        for model in af:

            start = (
                model.get(
                    "uniprotStart"
                )
                or
                model.get(
                    "sequenceStart"
                )
                or
                1
            )

            end = (
                model.get(
                    "uniprotEnd"
                )
                or
                model.get(
                    "sequenceEnd"
                )
                or
                len(
                    seq
                )
            )

            try:
                start = int(
                    start
                )
                end = int(
                    end
                )

            except Exception:
                continue

            if (
                start
                <= pos
                <= end
            ):
                covering.append(
                    model
                )


        rec[
            "alphafold_available"
        ] = bool(
            covering
        )


        if covering:

            model = covering[
                0
            ]

            rec[
                "alphafold_entry_id"
            ] = model.get(
                "entryId",
                ""
            )

            rec[
                "alphafold_global_metric"
            ] = model.get(
                "globalMetricValue",
                np.nan,
            )


            pdb_url = model.get(
                "pdbUrl",
                ""
            )

            rec[
                "alphafold_pdb_url"
            ] = pdb_url


            if pdb_url:

                p = (
                    RAW_AF
                    / (
                        f"{gene}_"
                        f"{site}_"
                        f"{acc}.pdb"
                    )
                )

                r = SESSION.get(
                    pdb_url,
                    timeout=120,
                )

                r.raise_for_status()

                p.write_bytes(
                    r.content
                )


                ca = parse_af_ca(
                    p
                )


                site_ca = ca[
                    ca[
                        "resnum"
                    ] == pos
                ]


                if len(
                    site_ca
                ) == 1:

                    sr = site_ca.iloc[
                        0
                    ]

                    rec[
                        "alphafold_site_residue"
                    ] = sr[
                        "residue"
                    ]

                    rec[
                        "alphafold_site_residue_matches"
                    ] = (
                        sr[
                            "residue"
                        ]
                        == expected_aa
                    )

                    rec[
                        "alphafold_site_plddt"
                    ] = float(
                        sr[
                            "plddt"
                        ]
                    )


                    local = ca[
                        ca[
                            "resnum"
                        ].between(
                            pos - 10,
                            pos + 10,
                        )
                    ]


                    rec[
                        "alphafold_local_21aa_median_plddt"
                    ] = float(
                        local[
                            "plddt"
                        ].median()
                    )


                    rec[
                        "alphafold_local_21aa_min_plddt"
                    ] = float(
                        local[
                            "plddt"
                        ].min()
                    )


                                                              
                    xyz = ca[
                        [
                            "x",
                            "y",
                            "z",
                        ]
                    ].to_numpy(
                        dtype=float
                    )

                    target = sr[
                        [
                            "x",
                            "y",
                            "z",
                        ]
                    ].to_numpy(
                        dtype=float
                    )

                    dist = np.sqrt(
                        (
                            (
                                xyz
                                - target
                            )
                            ** 2
                        ).sum(
                            axis=1
                        )
                    )


                    rec[
                        "alphafold_CA_neighbors_within_10A"
                    ] = int(
                        (
                            (
                                dist
                                <= 10
                            )
                            &
                            (
                                dist
                                > 0
                            )
                        ).sum()
                    )

                else:

                    rec[
                        "alphafold_site_lookup_count"
                    ] = len(
                        site_ca
                    )


                                                              
                                          
                                                              

    try:

        rcsb = rcsb_uniprot_search(
            acc
        )

        with open(
            RAW_RCSB
            / f"{gene}_{site}_{acc}.json",
            "w",
        ) as fh:

            json.dump(
                rcsb,
                fh,
                indent=2,
            )


        identifiers = [
            x.get(
                "identifier",
                ""
            )
            for x in rcsb.get(
                "result_set",
                []
            )
        ]


        pdb_ids = sorted(
            {
                x.split(
                    "_",
                    1,
                )[0].upper()

                for x in identifiers

                if x
            }
        )


        rec[
            "n_experimental_PDB_entries"
        ] = len(
            pdb_ids
        )

        rec[
            "experimental_PDB_entries"
        ] = ";".join(
            pdb_ids
        )


    except Exception as e:

        rec[
            "rcsb_error"
        ] = repr(
            e
        )


    rows.append(
        rec
    )


    time.sleep(
        0.2
    )


                                                              
           
                                                              

out = pd.DataFrame(
    rows
)


                                 
                               
def local_confidence(v):

    if pd.isna(
        v
    ):
        return "NA"

    v = float(
        v
    )

    if v >= 90:
        return "very_high"

    if v >= 70:
        return "confident"

    if v >= 50:
        return "low"

    return "very_low"


out[
    "alphafold_site_confidence_class"
] = out.get(
    "alphafold_site_plddt",
    pd.Series(
        np.nan,
        index=out.index,
    ),
).apply(
    local_confidence
)


out.to_csv(
    OUT
    / "STRUCTURAL_TRACTABILITY_PHASE1.tsv",
    sep="\t",
    index=False,
)


print()
print(
    "=" * 150
)
print(
    "STRUCTURAL TRACTABILITY — PHASE 1"
)
print(
    "=" * 150
)


show = [
    "gene",
    "site",
    "archetype",
    "uniprot_accession",
    "sequence_length",
    "uniprot_site_residue",
    "sequence_site_matches",
    "uniprot_features_at_site",
    "uniprot_modified_residue_annotation",
    "alphafold_available",
    "alphafold_site_plddt",
    "alphafold_local_21aa_median_plddt",
    "alphafold_site_confidence_class",
    "alphafold_CA_neighbors_within_10A",
    "n_experimental_PDB_entries",
    "experimental_PDB_entries",
]


show = [
    c
    for c in show
    if c in out.columns
]


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
    / "STRUCTURAL_TRACTABILITY_PHASE1.tsv"
)
