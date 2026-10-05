#!/usr/bin/env python3

from pathlib import Path
import json
import re
import requests
import pandas as pd


ROOT = Path(
    "/"
)

OUT = ROOT / "40_structural_tractability"

PHASE1 = (
    OUT
    / "STRUCTURAL_TRACTABILITY_PHASE1.tsv"
)

RCSB_RAW = (
    OUT
    / "raw"
    / "rcsb"
)

PDBe_RAW = (
    OUT
    / "raw"
    / "pdbe_mapping"
)

PDBe_RAW.mkdir(
    parents=True,
    exist_ok=True,
)


session = requests.Session()

session.headers.update({
    "User-Agent":
        "stefan-yeast-structure-audit/1.0"
})


x = pd.read_csv(
    PHASE1,
    sep="\t",
)


def rcsb_entities(
    gene,
    site,
    accession,
):

    p = (
        RCSB_RAW
        / f"{gene}_{site}_{accession}.json"
    )

    if not p.exists():
        return []

    j = json.loads(
        p.read_text()
    )

    out = []

    for item in j.get(
        "result_set",
        []
    ):

        identifier = str(
            item.get(
                "identifier",
                ""
            )
        )

        m = re.fullmatch(
            r"([0-9A-Za-z]{4})_(\d+)",
            identifier,
        )

        if m:

            out.append(
                (
                    m.group(1).upper(),
                    int(m.group(2)),
                )
            )

    return sorted(
        set(
            out
        )
    )


def compact_keys(
    obj,
    prefix="",
    depth=0,
    limit=4,
):

    """
    Print the JSON shape without assuming a PDBe
    response schema prematurely.
    """

    rows = []

    if depth > limit:
        return rows

    if isinstance(
        obj,
        dict,
    ):

        for k, v in obj.items():

            key = (
                f"{prefix}.{k}"
                if prefix
                else str(k)
            )

            rows.append(
                (
                    key,
                    type(v).__name__,
                )
            )

            rows.extend(
                compact_keys(
                    v,
                    key,
                    depth + 1,
                    limit,
                )
            )

    elif isinstance(
        obj,
        list,
    ) and obj:

        rows.extend(
            compact_keys(
                obj[0],
                prefix + "[0]",
                depth + 1,
                limit,
            )
        )

    return rows


records = []


for _, r in x.iterrows():

    gene = str(
        r["gene"]
    )

    site = str(
        r["site"]
    )

    accession = str(
        r.get(
            "uniprot_accession",
            ""
        )
    )

    pdb_text = str(
        r.get(
            "experimental_PDB_entries",
            ""
        )
    )

    if (
        not pdb_text
        or pdb_text.lower() == "nan"
    ):
        continue


    site_pos = int(
        r["position"]
    )


    entities = rcsb_entities(
        gene,
        site,
        accession,
    )


    print()
    print(
        "=" * 100
    )
    print(
        gene,
        site,
        accession,
    )
    print(
        "RCSB polymer entities:",
        entities,
    )
    print(
        "=" * 100
    )


    expected_pdbs = {
        z.strip().upper()
        for z in pdb_text.split(";")
        if z.strip()
    }


    for pdb_id, entity_id in entities:

        if pdb_id not in expected_pdbs:
            continue


        url = (
            "https://www.ebi.ac.uk/pdbe/api/"
            "pdb/entry/uniprot_mapping/"
            f"{pdb_id.lower()}/{entity_id}"
        )


        print()
        print(
            f"{pdb_id} entity {entity_id}"
        )
        print(
            url
        )


        try:

            resp = session.get(
                url,
                timeout=60,
            )


            print(
                "HTTP:",
                resp.status_code,
            )


            rec = {
                "gene":
                    gene,

                "site":
                    site,

                "uniprot_accession":
                    accession,

                "uniprot_position":
                    site_pos,

                "pdb_id":
                    pdb_id,

                "entity_id":
                    entity_id,

                "http_status":
                    resp.status_code,

                "url":
                    url,
            }


            if resp.status_code == 200:

                data = resp.json()

                raw_path = (
                    PDBe_RAW
                    / (
                        f"{gene}_{site}_"
                        f"{pdb_id}_entity{entity_id}.json"
                    )
                )

                raw_path.write_text(
                    json.dumps(
                        data,
                        indent=2,
                    )
                )


                print(
                    "JSON SHAPE:"
                )

                shape = compact_keys(
                    data
                )

                for key, typ in shape[:100]:

                    print(
                        " ",
                        key,
                        "=>",
                        typ,
                    )


                print()
                print(
                    "FIRST 3000 JSON CHARACTERS:"
                )

                txt = json.dumps(
                    data,
                    indent=2,
                )

                print(
                    txt[:3000]
                )


                rec[
                    "raw_json"
                ] = str(
                    raw_path
                )

            else:

                print(
                    resp.text[:1000]
                )


            records.append(
                rec
            )


        except Exception as e:

            print(
                "ERROR:",
                repr(e),
            )

            records.append({
                "gene":
                    gene,

                "site":
                    site,

                "uniprot_accession":
                    accession,

                "uniprot_position":
                    site_pos,

                "pdb_id":
                    pdb_id,

                "entity_id":
                    entity_id,

                "error":
                    repr(e),
            })


out = pd.DataFrame(
    records
)

out.to_csv(
    OUT
    / "EXPERIMENTAL_STRUCTURE_MAPPING_AUDIT.tsv",
    sep="\t",
    index=False,
)


print()
print(
    "=" * 100
)
print(
    "AUDIT SUMMARY"
)
print(
    "=" * 100
)

if len(out):

    print(
        out.to_string(
            index=False
        )
    )

print()
print(
    "Output:",
    OUT
    / "EXPERIMENTAL_STRUCTURE_MAPPING_AUDIT.tsv"
)
