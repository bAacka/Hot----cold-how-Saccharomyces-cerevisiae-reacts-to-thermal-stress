#!/usr/bin/env python3

from pathlib import Path
import hashlib
import sys

import pandas as pd
import requests


ROOT = Path("/")

INV = (
    ROOT
    / "38_pka_mechanistic_resources"
    / "inventory"
    / "ALL_PRIDE_FILE_INVENTORY.tsv"
)

OUT = (
    ROOT
    / "38_pka_mechanistic_resources"
    / "processed_downloads"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                            
WANTED = {

    "PXD052971": [
        "R_analysis_final_240607.R",
        "R_analysis_targeted_240608.R",
        "phospho_large_240607.csv",
        "peptide_annotation.csv",
        "PRM_peptide_results_long_240608.csv",
        "Tukey_table_2400607.csv",
    ],

    "PXD052995": [
        "Annotation_file.xlsx",
        "Phospho_Sites.txt",
        "Supp_Table_3_FU_241117.xlsx",
        "Transition_results_230915.csv",
        "data_analysis_240610.R",
        "final_table_241117.csv",
        "parameters.txt",
        "summary.txt",
    ],

    "PXD059338": [
        "PXD059338_community_annotated.sdrf.tsv",
        "Annotation_file.xlsx",
        "data_analysis.R",
        "Phospho_STY_Sites.txt",
        "parameters.txt",
        "summary.txt",
    ],
}


df = pd.read_csv(
    INV,
    sep="\t",
)


session = requests.Session()

session.headers.update({
    "User-Agent":
        "stefan-yeast-rebuild/1.0"
})


manifest = []


for project, names in WANTED.items():

    project_out = OUT / project

    project_out.mkdir(
        parents=True,
        exist_ok=True,
    )

    sub = df[
        df["project"] == project
    ].copy()

    for name in names:

        hit = sub[
            sub["file_name"] == name
        ]

        if len(hit) != 1:
            raise RuntimeError(
                f"{project} / {name}: "
                f"expected exactly one inventory row, found {len(hit)}"
            )

        r = hit.iloc[0]

        ftp = str(
            r["ftp"]
        )

        if not ftp.startswith(
            "ftp://ftp.pride.ebi.ac.uk/"
        ):
            raise RuntimeError(
                f"Unexpected PRIDE URL: {ftp}"
            )

        url = ftp.replace(
            "ftp://ftp.pride.ebi.ac.uk/",
            "https://ftp.pride.ebi.ac.uk/",
            1,
        )

        dest = project_out / name

        print(
            f"Downloading {project}/{name}"
        )

        with session.get(
            url,
            stream=True,
            timeout=120,
        ) as resp:

            resp.raise_for_status()

            h = hashlib.sha256()

            nbytes = 0

            with open(
                dest,
                "wb",
            ) as fh:

                for chunk in resp.iter_content(
                    chunk_size=1024 * 1024
                ):

                    if not chunk:
                        continue

                    fh.write(
                        chunk
                    )

                    h.update(
                        chunk
                    )

                    nbytes += len(
                        chunk
                    )

        expected = r.get(
            "file_size_bytes",
            None,
        )

        if pd.notna(
            expected
        ):

            expected = int(
                expected
            )

            if nbytes != expected:
                raise RuntimeError(
                    f"Size mismatch for {dest}: "
                    f"{nbytes} != {expected}"
                )

        manifest.append({
            "project":
                project,

            "file_name":
                name,

            "bytes":
                nbytes,

            "sha256":
                h.hexdigest(),

            "source":
                url,
        })


manifest = pd.DataFrame(
    manifest
)


manifest.to_csv(
    OUT / "DOWNLOAD_MANIFEST.tsv",
    sep="\t",
    index=False,
)


print()
print(
    manifest[
        [
            "project",
            "file_name",
            "bytes",
            "sha256",
        ]
    ].to_string(
        index=False
    )
)

print()
print(
    "Output:",
    OUT,
)
