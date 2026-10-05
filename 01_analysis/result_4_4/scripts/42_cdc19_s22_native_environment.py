#!/usr/bin/env python3

from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
import requests

from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley


ROOT = Path(
    "/"
)

CIFDIR = (
    ROOT
    / "40_structural_tractability"
    / "raw"
    / "exact_residue_coverage"
)

OUT = (
    ROOT
    / "41_cdc19_s22_structure"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


STRUCTURES = [
    "1A3W",
    "1A3X",
]

TARGET_RESI = 22
TARGET_RESNAME = "SER"

                                          
                                                
SER_MAX_ASA = 155.0


session = requests.Session()
session.headers.update({
    "User-Agent":
        "stefan-yeast-cdc19-structure/1.0"
})


def heavy(atom):
    element = str(
        atom.element
    ).strip().upper()

    return element != "H"


def residue_label(residue):

    chain = residue.get_parent()

    hetflag, resseq, icode = (
        residue.id
    )

    return (
        chain.id,
        str(resseq),
        str(icode).strip(),
        residue.resname,
    )


def find_target_residues(structure):

    hits = []

    for model in structure:

                                          
        for chain in model:

            for residue in chain:

                if (
                    residue.id[0] == " "
                    and
                    residue.id[1] == TARGET_RESI
                    and
                    residue.resname == TARGET_RESNAME
                ):
                    hits.append(
                        residue
                    )

        break

    return hits


def min_atom_distance(
    target_atom,
    residue,
):

    best = None

    for atom in residue.get_atoms():

        if not heavy(
            atom
        ):
            continue

        d = float(
            np.linalg.norm(
                target_atom.coord
                - atom.coord
            )
        )

        if (
            best is None
            or d < best[0]
        ):
            best = (
                d,
                atom.name,
            )

    return best


def entry_metadata(pdb_id):

    url = (
        "https://data.rcsb.org/rest/v1/core/entry/"
        f"{pdb_id}"
    )

    r = session.get(
        url,
        timeout=60,
    )

    r.raise_for_status()

    j = r.json()

    info = (
        j.get(
            "rcsb_entry_info",
            {}
        )
        or {}
    )

    exptl = (
        j.get(
            "exptl",
            []
        )
        or []
    )

    struct = (
        j.get(
            "struct",
            {}
        )
        or {}
    )

    return {
        "pdb_id":
            pdb_id,

        "title":
            struct.get(
                "title",
                ""
            ),

        "experimental_method":
            ";".join(
                str(x.get("method", ""))
                for x in exptl
            ),

        "resolution_A":
            (
                info.get(
                    "resolution_combined",
                    [np.nan]
                )[0]
                if info.get(
                    "resolution_combined"
                )
                else np.nan
            ),
    }


instance_rows = []
contact_rows = []
polar_rows = []
hetero_rows = []
metadata_rows = []


for pdb_id in STRUCTURES:

    cif = (
        CIFDIR
        / f"{pdb_id}.cif"
    )

    if not cif.exists():
        raise RuntimeError(
            f"Missing mmCIF: {cif}"
        )


    print()
    print(
        "=" * 100
    )
    print(
        pdb_id
    )
    print(
        "=" * 100
    )


                                                              
              
                                                              

    try:
        meta = entry_metadata(
            pdb_id
        )

    except Exception as e:

        meta = {
            "pdb_id":
                pdb_id,

            "metadata_error":
                repr(e),
        }

    metadata_rows.append(
        meta
    )


                                                              
                       
                                                              

    parser = MMCIFParser(
        QUIET=True,
    )

    structure = parser.get_structure(
        pdb_id,
        str(cif),
    )


                                   
    sr = ShrakeRupley(
        probe_radius=1.4,
        n_points=960,
    )

    sr.compute(
        structure,
        level="R",
    )


    targets = find_target_residues(
        structure
    )


    print(
        "S22 structural instances:",
        len(targets),
    )


    if not targets:
        continue


    model = next(
        structure.get_models()
    )


    all_residues = [
        residue
        for chain in model
        for residue in chain
    ]


    for target in targets:

        chain = target.get_parent()

        if "OG" not in target:
            raise RuntimeError(
                f"{pdb_id} chain {chain.id} "
                "S22 lacks OG atom"
            )


        og = target[
            "OG"
        ]

        cb = (
            target["CB"]
            if "CB" in target
            else None
        )


        residue_sasa = float(
            getattr(
                target,
                "sasa",
                np.nan,
            )
        )


        rsa = (
            residue_sasa
            / SER_MAX_ASA
            if np.isfinite(
                residue_sasa
            )
            else np.nan
        )


        instance = {
            "pdb_id":
                pdb_id,

            "chain":
                chain.id,

            "auth_residue":
                TARGET_RESI,

            "resname":
                target.resname,

            "residue_SASA_A2":
                residue_sasa,

            "approx_relative_SASA":
                rsa,

            "OG_x":
                float(
                    og.coord[0]
                ),

            "OG_y":
                float(
                    og.coord[1]
                ),

            "OG_z":
                float(
                    og.coord[2]
                ),
        }


                                                              
                                             
                                                              

        local_contacts = []

        for residue in all_residues:

            if residue is target:
                continue

            best = min_atom_distance(
                og,
                residue,
            )

            if best is None:
                continue

            dist, atom_name = best

            if dist <= 8.0:

                other_chain = (
                    residue
                    .get_parent()
                    .id
                )

                hetflag, resseq, icode = (
                    residue.id
                )

                row = {
                    "pdb_id":
                        pdb_id,

                    "target_chain":
                        chain.id,

                    "target_site":
                        "S22",

                    "neighbor_chain":
                        other_chain,

                    "neighbor_resseq":
                        resseq,

                    "neighbor_icode":
                        str(
                            icode
                        ).strip(),

                    "neighbor_resname":
                        residue.resname,

                    "neighbor_hetflag":
                        str(
                            hetflag
                        ),

                    "nearest_atom":
                        atom_name,

                    "OG_min_distance_A":
                        dist,

                    "same_chain":
                        (
                            other_chain
                            == chain.id
                        ),
                }

                local_contacts.append(
                    row
                )


                if dist <= 5.0:
                    contact_rows.append(
                        row
                    )


                                             
                if dist <= 6.0:

                    charged_or_polar = {
                        "ASP",
                        "GLU",
                        "LYS",
                        "ARG",
                        "HIS",
                        "ASN",
                        "GLN",
                        "SER",
                        "THR",
                        "TYR",
                        "CYS",
                    }

                    if (
                        residue.resname
                        in charged_or_polar
                    ):

                        polar_rows.append(
                            row
                        )


                                                      
                if hetflag != " ":

                    hetero_rows.append(
                        row
                    )


        if local_contacts:

            lc = pd.DataFrame(
                local_contacts
            )


            instance[
                "nearest_nonself_distance_A"
            ] = float(
                lc[
                    "OG_min_distance_A"
                ].min()
            )


            instance[
                "n_residue_contacts_5A"
            ] = int(
                (
                    lc[
                        "OG_min_distance_A"
                    ]
                    <= 5.0
                ).sum()
            )


            instance[
                "n_residue_neighbors_8A"
            ] = len(
                lc
            )


            instance[
                "n_interchain_neighbors_8A"
            ] = int(
                (
                    ~lc[
                        "same_chain"
                    ]
                ).sum()
            )


            instance[
                "n_interchain_contacts_5A"
            ] = int(
                (
                    (
                        ~lc[
                            "same_chain"
                        ]
                    )
                    &
                    (
                        lc[
                            "OG_min_distance_A"
                        ]
                        <= 5.0
                    )
                ).sum()
            )


            nearest = (
                lc.sort_values(
                    "OG_min_distance_A"
                )
                .head(
                    8
                )
            )


            instance[
                "nearest_residues"
            ] = " | ".join(
                (
                    nearest[
                        "neighbor_chain"
                    ].astype(str)
                    + ":"
                    + nearest[
                        "neighbor_resname"
                    ].astype(str)
                    + nearest[
                        "neighbor_resseq"
                    ].astype(str)
                    + "@"
                    + nearest[
                        "OG_min_distance_A"
                    ].map(
                        lambda x:
                            f"{x:.2f}"
                    )
                )
            )


                                                              
                                                     
                                                              

        instance[
            "possible_interchain_interface_by_5A"
        ] = (
            instance.get(
                "n_interchain_contacts_5A",
                0,
            )
            > 0
        )


        instance_rows.append(
            instance
        )


        print(
            "chain",
            chain.id,
            "| SASA:",
            round(
                residue_sasa,
                2
            ),
            "| RSA:",
            round(
                rsa,
                3
            ),
            "| contacts <=5A:",
            instance.get(
                "n_residue_contacts_5A",
                0
            ),
            "| interchain <=5A:",
            instance.get(
                "n_interchain_contacts_5A",
                0
            ),
        )

        print(
            " nearest:",
            instance.get(
                "nearest_residues",
                ""
            ),
        )


                                                              
               
                                                              

instances = pd.DataFrame(
    instance_rows
)

contacts = pd.DataFrame(
    contact_rows
)

polar = pd.DataFrame(
    polar_rows
)

hetero = pd.DataFrame(
    hetero_rows
)

metadata = pd.DataFrame(
    metadata_rows
)


instances.to_csv(
    OUT
    / "CDC19_S22_STRUCTURE_INSTANCES.tsv",
    sep="\t",
    index=False,
)

contacts.to_csv(
    OUT
    / "CDC19_S22_OG_CONTACTS_5A.tsv",
    sep="\t",
    index=False,
)

polar.to_csv(
    OUT
    / "CDC19_S22_POLAR_NEIGHBORS_6A.tsv",
    sep="\t",
    index=False,
)

hetero.to_csv(
    OUT
    / "CDC19_S22_HETERO_NEIGHBORS_8A.tsv",
    sep="\t",
    index=False,
)

metadata.to_csv(
    OUT
    / "CDC19_STRUCTURE_METADATA.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                   
                                                              

if len(
    contacts
):

    c = contacts.copy()

    c[
        "neighbor_identity"
    ] = (
        c[
            "neighbor_resname"
        ].astype(str)
        +
        c[
            "neighbor_resseq"
        ].astype(str)
    )


    persistence = (
        c.groupby(
            "neighbor_identity",
            as_index=False,
        )
        .agg(
            structures=(
                "pdb_id",
                lambda z:
                    ";".join(
                        sorted(
                            set(
                                z
                            )
                        )
                    ),
            ),

            target_chains=(
                "target_chain",
                lambda z:
                    ";".join(
                        sorted(
                            set(
                                z
                            )
                        )
                    ),
            ),

            min_distance_A=(
                "OG_min_distance_A",
                "min",
            ),

            max_distance_A=(
                "OG_min_distance_A",
                "max",
            ),

            median_distance_A=(
                "OG_min_distance_A",
                "median",
            ),

            n_observations=(
                "OG_min_distance_A",
                "size",
            ),
        )
    )


    persistence[
        "n_structures"
    ] = (
        persistence[
            "structures"
        ]
        .str.split(
            ";"
        )
        .str.len()
    )


    persistence = persistence.sort_values(
        [
            "n_structures",
            "n_observations",
            "median_distance_A",
        ],
        ascending=[
            False,
            False,
            True,
        ],
    )


    persistence.to_csv(
        OUT
        / "CDC19_S22_CONTACT_PERSISTENCE.tsv",
        sep="\t",
        index=False,
    )


                                                              
                 
                                                              

print()
print(
    "=" * 120
)
print(
    "CDC19 S22 NATIVE STRUCTURAL ENVIRONMENT"
)
print(
    "=" * 120
)


print()
print(
    "STRUCTURE METADATA"
)

print(
    metadata.to_string(
        index=False
    )
)


print()
print(
    "S22 INSTANCES"
)

print(
    instances.to_string(
        index=False
    )
)


print()
print(
    "CLOSE CONTACTS <= 5 A"
)

if len(
    contacts
):
    print(
        contacts.sort_values(
            [
                "pdb_id",
                "target_chain",
                "OG_min_distance_A",
            ]
        ).to_string(
            index=False
        )
    )
else:
    print(
        "none"
    )


print()
print(
    "POLAR / CHARGED NEIGHBORS <= 6 A"
)

if len(
    polar
):
    print(
        polar.sort_values(
            [
                "pdb_id",
                "target_chain",
                "OG_min_distance_A",
            ]
        ).to_string(
            index=False
        )
    )
else:
    print(
        "none"
    )


print()
print(
    "HETERO / LIGAND NEIGHBORS <= 8 A"
)

if len(
    hetero
):
    print(
        hetero.sort_values(
            [
                "pdb_id",
                "target_chain",
                "OG_min_distance_A",
            ]
        ).to_string(
            index=False
        )
    )
else:
    print(
        "none"
    )


print()
print(
    "Output:",
    OUT,
)
