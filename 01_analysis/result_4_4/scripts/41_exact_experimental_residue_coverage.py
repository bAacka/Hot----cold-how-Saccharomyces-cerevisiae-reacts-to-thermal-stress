#!/usr/bin/env python3

from pathlib import Path
import json
import math

import pandas as pd
import requests


try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except ImportError:
    raise SystemExit(
        "Biopython is required for mmCIF parsing.\n"
        "Install with:\n"
        "  python -m pip install biopython"
    )


ROOT = Path(
    "/"
)

OUT = ROOT / "40_structural_tractability"

RAW = (
    OUT
    / "raw"
    / "exact_residue_coverage"
)

RAW.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                        
                                                              

TARGETS = [
    {
        "gene": "CDC19",
        "site": "S22",
        "uniprot": "P00549",
        "position": 22,
        "expected_residue": "S",
        "expected_entities": [
            "1A3W_1",
            "1A3X_1",
        ],
    },

    {
        "gene": "MAF1",
        "site": "S90",
        "uniprot": "P41910",
        "position": 90,
        "expected_residue": "S",
        "expected_entities": [
            "6TUT_18",
        ],
    },

    {
        "gene": "NUP60",
        "site": "S10",
        "uniprot": "P39705",
        "position": 10,
        "expected_residue": "S",
        "expected_entities": [
            "8ZZA_30",
            "8ZZB_30",
            "8ZZC_30",
            "9A8M_4",
            "9A8N_4",
        ],
    },
]


SESSION = requests.Session()

SESSION.headers.update({
    "User-Agent":
        "stefan-yeast-exact-structural-coverage/1.0"
})


GRAPHQL = """
query ExactUniProtToPDB(
  $id: String!,
  $range: [Int!]
) {
  alignments(
    from: UNIPROT,
    to: PDB_ENTITY,
    queryId: $id,
    range: $range
  ) {
    query_sequence
    target_alignments {
      target_id
      target_sequence
      coverage {
        query_coverage
        query_length
        target_coverage
        target_length
      }
      aligned_regions {
        query_begin
        query_end
        target_begin
        target_end
      }
    }
  }
}
"""


def as_list(x):
    if x is None:
        return []

    if isinstance(
        x,
        list,
    ):
        return x

    return [x]


def graphql_alignment(
    uniprot,
    position,
):

    payload = {
        "query":
            GRAPHQL,

        "variables": {
            "id":
                uniprot,

            "range": [
                position,
                position,
            ],
        },
    }

    r = SESSION.post(
        "https://sequence-coordinates.rcsb.org/graphql",
        json=payload,
        timeout=90,
    )

    r.raise_for_status()

    data = r.json()

    if data.get(
        "errors"
    ):
        raise RuntimeError(
            json.dumps(
                data[
                    "errors"
                ],
                indent=2,
            )
        )

    return data


def map_uniprot_position(
    target_alignment,
    position,
):
    """
    Protein-protein aligned regions map one residue
    at a time within an aligned block.

    Return the PDB entity sequence position corresponding
    to the requested UniProt residue.
    """

    for region in target_alignment.get(
        "aligned_regions",
        []
    ):

        qb = int(
            region[
                "query_begin"
            ]
        )

        qe = int(
            region[
                "query_end"
            ]
        )

        tb = int(
            region[
                "target_begin"
            ]
        )

        te = int(
            region[
                "target_end"
            ]
        )


        if (
            qb <= position <= qe
        ):

            qspan = (
                qe - qb
            )

            tspan = (
                te - tb
            )

            if qspan != tspan:
                raise RuntimeError(
                    "Non-1:1 protein alignment block "
                    f"for {target_alignment['target_id']}: "
                    f"{region}"
                )

            return (
                tb
                + (
                    position
                    - qb
                )
            )

    return None


def download_cif(
    pdb_id,
):

    dest = (
        RAW
        / f"{pdb_id}.cif"
    )

    if dest.exists():
        return dest

    url = (
        "https://files.rcsb.org/download/"
        f"{pdb_id}.cif"
    )

    r = SESSION.get(
        url,
        timeout=180,
    )

    r.raise_for_status()

    dest.write_bytes(
        r.content
    )

    return dest


def get_col(
    cif,
    key,
):

    if key not in cif:
        return []

    return as_list(
        cif[
            key
        ]
    )


def coordinate_rows_for_site(
    cif_path,
    entity_id,
    entity_seq_pos,
):

    cif = MMCIF2Dict(
        str(
            cif_path
        )
    )


    groups = get_col(
        cif,
        "_atom_site.group_PDB",
    )

    entity_ids = get_col(
        cif,
        "_atom_site.label_entity_id",
    )

    seq_ids = get_col(
        cif,
        "_atom_site.label_seq_id",
    )

    label_asym = get_col(
        cif,
        "_atom_site.label_asym_id",
    )

    auth_asym = get_col(
        cif,
        "_atom_site.auth_asym_id",
    )

    auth_seq = get_col(
        cif,
        "_atom_site.auth_seq_id",
    )

    comp = get_col(
        cif,
        "_atom_site.label_comp_id",
    )

    atom = get_col(
        cif,
        "_atom_site.label_atom_id",
    )

    alt = get_col(
        cif,
        "_atom_site.label_alt_id",
    )


    n = len(
        groups
    )


    required = {
        "entity":
            len(
                entity_ids
            ),

        "seq":
            len(
                seq_ids
            ),

        "label_asym":
            len(
                label_asym
            ),

        "auth_asym":
            len(
                auth_asym
            ),

        "auth_seq":
            len(
                auth_seq
            ),

        "comp":
            len(
                comp
            ),

        "atom":
            len(
                atom
            ),
    }


    bad = {
        k: v
        for k, v in required.items()
        if v != n
    }


    if bad:
        raise RuntimeError(
            f"Unexpected _atom_site column lengths "
            f"in {cif_path}: "
            f"n={n}, bad={bad}"
        )


    rows = []


    for i in range(
        n
    ):

        if groups[i] not in {
            "ATOM",
            "HETATM",
        }:
            continue

        if str(
            entity_ids[i]
        ) != str(
            entity_id
        ):
            continue

        try:
            seq = int(
                seq_ids[i]
            )
        except Exception:
            continue

        if seq != int(
            entity_seq_pos
        ):
            continue


        rows.append({
            "label_asym_id":
                label_asym[i],

            "auth_asym_id":
                auth_asym[i],

            "auth_seq_id":
                auth_seq[i],

            "residue_name":
                comp[i],

            "atom_name":
                atom[i],

            "alt_id":
                (
                    alt[i]
                    if i < len(
                        alt
                    )
                    else ""
                ),
        })


    return rows


                                                              
     
                                                              

records = []


for target in TARGETS:

    gene = target[
        "gene"
    ]

    site = target[
        "site"
    ]

    accession = target[
        "uniprot"
    ]

    position = target[
        "position"
    ]

    expected_entities = set(
        target[
            "expected_entities"
        ]
    )


    print()
    print(
        "=" * 100
    )
    print(
        gene,
        site,
        accession,
        "UniProt position",
        position,
    )
    print(
        "=" * 100
    )


    data = graphql_alignment(
        accession,
        position,
    )


    raw_json = (
        RAW
        / f"{gene}_{site}_{accession}_alignment.json"
    )

    raw_json.write_text(
        json.dumps(
            data,
            indent=2,
        )
    )


    align = (
        data.get(
            "data",
            {}
        )
        .get(
            "alignments"
        )
    )


    if align is None:

        print(
            "No alignment object returned."
        )

        for entity in sorted(
            expected_entities
        ):

            pdb_id, entity_id = (
                entity.split(
                    "_",
                    1,
                )
            )

            records.append({
                "gene":
                    gene,

                "site":
                    site,

                "uniprot_accession":
                    accession,

                "uniprot_position":
                    position,

                "pdb_entity":
                    entity,

                "pdb_id":
                    pdb_id,

                "entity_id":
                    entity_id,

                "site_in_entity_alignment":
                    False,

                "coordinate_present":
                    False,
            })

        continue


    targets = {
        x[
            "target_id"
        ].upper():
            x

        for x in align.get(
            "target_alignments",
            []
        )
    }


    print(
        "Entities returned for exact residue:",
        sorted(
            targets
        ),
    )


    for expected in sorted(
        expected_entities
    ):

        pdb_id, entity_id = (
            expected.split(
                "_",
                1,
            )
        )


        rec = {
            "gene":
                gene,

            "site":
                site,

            "uniprot_accession":
                accession,

            "uniprot_position":
                position,

            "expected_residue":
                target[
                    "expected_residue"
                ],

            "pdb_entity":
                expected,

            "pdb_id":
                pdb_id,

            "entity_id":
                entity_id,
        }


        ta = targets.get(
            expected.upper()
        )


        if ta is None:

            rec[
                "site_in_entity_alignment"
            ] = False

            rec[
                "coordinate_present"
            ] = False

            records.append(
                rec
            )

            print(
                expected,
                "=> exact UniProt site NOT in entity alignment",
            )

            continue


        entity_pos = map_uniprot_position(
            ta,
            position,
        )


        rec[
            "site_in_entity_alignment"
        ] = (
            entity_pos
            is not None
        )

        rec[
            "entity_sequence_position"
        ] = entity_pos


        rec[
            "aligned_regions"
        ] = json.dumps(
            ta.get(
                "aligned_regions",
                []
            ),
            separators=(
                ",",
                ":",
            ),
        )


        cov = ta.get(
            "coverage",
            {}
        ) or {}


        for k in [
            "query_coverage",
            "query_length",
            "target_coverage",
            "target_length",
        ]:

            rec[
                f"coverage_{k}"
            ] = cov.get(
                k
            )


        if entity_pos is None:

            rec[
                "coordinate_present"
            ] = False

            records.append(
                rec
            )

            print(
                expected,
                "=> alignment object exists but site is not mapped",
            )

            continue


        cif = download_cif(
            pdb_id
        )


        atoms = coordinate_rows_for_site(
            cif,
            entity_id,
            entity_pos,
        )


        rec[
            "coordinate_present"
        ] = bool(
            atoms
        )

        rec[
            "n_atoms_at_site"
        ] = len(
            atoms
        )


        if atoms:

            rec[
                "label_asym_ids"
            ] = ";".join(
                sorted(
                    {
                        str(
                            z[
                                "label_asym_id"
                            ]
                        )
                        for z in atoms
                    }
                )
            )

            rec[
                "auth_asym_ids"
            ] = ";".join(
                sorted(
                    {
                        str(
                            z[
                                "auth_asym_id"
                            ]
                        )
                        for z in atoms
                    }
                )
            )

            rec[
                "auth_seq_ids"
            ] = ";".join(
                sorted(
                    {
                        str(
                            z[
                                "auth_seq_id"
                            ]
                        )
                        for z in atoms
                    }
                )
            )

            rec[
                "residue_names"
            ] = ";".join(
                sorted(
                    {
                        str(
                            z[
                                "residue_name"
                            ]
                        )
                        for z in atoms
                    }
                )
            )

            rec[
                "n_structural_instances"
            ] = len(
                {
                    (
                        z[
                            "label_asym_id"
                        ],
                        z[
                            "auth_asym_id"
                        ],
                    )
                    for z in atoms
                }
            )


        records.append(
            rec
        )


        print(
            expected,
            "=> entity pos",
            entity_pos,
            "| coordinates:",
            bool(
                atoms
            ),
            "| atoms:",
            len(
                atoms
            ),
            "| auth chains:",
            rec.get(
                "auth_asym_ids",
                "",
            ),
            "| auth residue:",
            rec.get(
                "auth_seq_ids",
                "",
            ),
        )


out = pd.DataFrame(
    records
)


out.to_csv(
    OUT
    / "EXACT_EXPERIMENTAL_RESIDUE_COVERAGE.tsv",
    sep="\t",
    index=False,
)


print()
print(
    "=" * 120
)
print(
    "EXACT EXPERIMENTAL RESIDUE COVERAGE"
)
print(
    "=" * 120
)


show = [
    "gene",
    "site",
    "pdb_entity",
    "site_in_entity_alignment",
    "entity_sequence_position",
    "coordinate_present",
    "n_atoms_at_site",
    "n_structural_instances",
    "label_asym_ids",
    "auth_asym_ids",
    "auth_seq_ids",
    "residue_names",
]


show = [
    c for c in show
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
    / "EXACT_EXPERIMENTAL_RESIDUE_COVERAGE.tsv"
)
