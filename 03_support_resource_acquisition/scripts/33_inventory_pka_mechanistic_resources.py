#!/usr/bin/env python3

from pathlib import Path
import json
import re
import sys

import pandas as pd
import requests


ROOT = Path("/")
OUT = ROOT / "38_pka_mechanistic_resources"
PRIDE_OUT = OUT / "pride"
INV_OUT = OUT / "inventory"

PRIDE_OUT.mkdir(parents=True, exist_ok=True)
INV_OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"

PROJECTS = {
    "PXD052971": "PKA-dependent phospho enrichment: exponential/day8/heat shock",
    "PXD052995": "in-vitro PKA phosphosite mapping",
    "PXD059338": "starvation phosphoproteomics",
}

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "stefan-yeast-rebuild/1.0"
})


def get_json(url, params=None):
    r = SESSION.get(
        url,
        params=params,
        timeout=60,
    )
    r.raise_for_status()
    return r.json(), dict(r.headers)


def unwrap_items(obj):
    """
    Support current bare-list PRIDE responses and older HAL-style shapes.
    """
    if isinstance(obj, list):
        return obj

    if not isinstance(obj, dict):
        return []

    if "_embedded" in obj and isinstance(obj["_embedded"], dict):
        for value in obj["_embedded"].values():
            if isinstance(value, list):
                return value

    for key in ("files", "items", "content"):
        if key in obj and isinstance(obj[key], list):
            return obj[key]

    return []


def pick(d, *keys):
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None


def nested_value(x):
    if isinstance(x, dict):
        return (
            x.get("value")
            or x.get("name")
            or x.get("label")
            or x.get("accession")
        )
    return x


def file_locations(rec):
    locations = pick(
        rec,
        "publicFileLocations",
        "fileLocations",
        "locations",
    )

    ftp = []
    aspera = []
    http = []

    if isinstance(locations, list):
        for x in locations:
            if not isinstance(x, dict):
                continue

            value = (
                x.get("value")
                or x.get("location")
                or x.get("url")
                or ""
            )

            typ = str(
                x.get("name")
                or x.get("type")
                or x.get("protocol")
                or ""
            ).lower()

            if not value:
                continue

            if "ftp" in typ or str(value).startswith("ftp"):
                ftp.append(str(value))
            elif "aspera" in typ:
                aspera.append(str(value))
            elif str(value).startswith(("http://", "https://")):
                http.append(str(value))

    return (
        ";".join(ftp),
        ";".join(aspera),
        ";".join(http),
    )


def get_all_files(accession):
    rows = []
    page = 0
    page_size = 100

    while True:
        obj, headers = get_json(
            f"{BASE}/projects/{accession}/files",
            params={
                "page": page,
                "pageSize": page_size,
            },
        )

        items = unwrap_items(obj)

        if not items:
            break

        rows.extend(items)

        if len(items) < page_size:
            break

        page += 1

        if page > 100:
            raise RuntimeError("Pagination runaway")

    return rows


all_flat = []
project_summary = []


for accession, role in PROJECTS.items():

    print()
    print("=" * 100)
    print(accession, "-", role)
    print("=" * 100)

    proj, _ = get_json(
        f"{BASE}/projects/{accession}"
    )

    with open(
        PRIDE_OUT / f"{accession}_project.json",
        "w",
    ) as fh:
        json.dump(
            proj,
            fh,
            indent=2,
        )

    files = get_all_files(
        accession
    )

    with open(
        PRIDE_OUT / f"{accession}_files.json",
        "w",
    ) as fh:
        json.dump(
            files,
            fh,
            indent=2,
        )

    flat = []

    for rec in files:

        ftp, aspera, http = file_locations(rec)

        category = nested_value(
            pick(
                rec,
                "fileCategory",
                "category",
            )
        )

        size = pick(
            rec,
            "fileSizeBytes",
            "fileSize",
            "size",
        )

        try:
            size = int(size)
        except Exception:
            size = None

        row = {
            "project":
                accession,

            "project_role":
                role,

            "file_name":
                pick(
                    rec,
                    "fileName",
                    "name",
                ),

            "file_category":
                category,

            "file_size_bytes":
                size,

            "file_size_MB":
                (
                    size / 1024**2
                    if size is not None
                    else None
                ),

            "file_accession":
                pick(
                    rec,
                    "accession",
                    "fileAccession",
                ),

            "ftp":
                ftp,

            "aspera":
                aspera,

            "http":
                http,
        }

        flat.append(row)
        all_flat.append(row)

    df = pd.DataFrame(flat)

    df.to_csv(
        INV_OUT / f"{accession}_FILE_INVENTORY.tsv",
        sep="\t",
        index=False,
    )

    n = len(df)

    total_gb = (
        df["file_size_bytes"]
        .fillna(0)
        .sum()
        / 1024**3
        if n
        else 0
    )

    print("files:", n)
    print("total size GB:", round(total_gb, 3))

    if n:
        print()
        print("categories:")
        print(
            df["file_category"]
            .fillna("UNKNOWN")
            .value_counts()
            .to_string()
        )

        print()
        print("candidate processed/result files:")

        mask = (
            ~df["file_name"]
            .fillna("")
            .str.lower()
            .str.endswith(".raw")
        )

        shown = df.loc[
            mask,
            [
                "file_name",
                "file_category",
                "file_size_MB",
            ],
        ].sort_values(
            [
                "file_category",
                "file_name",
            ]
        )

        print(
            shown.to_string(
                index=False
            )
        )

    project_summary.append({
        "project":
            accession,
        "role":
            role,
        "n_files":
            n,
        "total_GB":
            total_gb,
    })


all_df = pd.DataFrame(
    all_flat
)

all_df.to_csv(
    INV_OUT / "ALL_PRIDE_FILE_INVENTORY.tsv",
    sep="\t",
    index=False,
)


summary = pd.DataFrame(
    project_summary
)

summary.to_csv(
    INV_OUT / "PRIDE_PROJECT_SUMMARY.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                         
                                                              

name = (
    all_df["file_name"]
    .fillna("")
    .str.lower()
)

category = (
    all_df["file_category"]
    .fillna("")
    .str.upper()
)

interesting = (
    category.isin(
        [
            "RESULT",
            "OTHER",
            "PEAK",
        ]
    )
    |
    name.str.endswith(
        (
            ".xlsx",
            ".xls",
            ".csv",
            ".tsv",
            ".txt",
            ".zip",
            ".gz",
            ".r",
            ".rmd",
            ".py",
            ".pdf",
        )
    )
)

                                                            
interesting &= ~name.str.endswith(
    ".raw"
)

cand = all_df[
    interesting
].copy()

cand = cand.sort_values(
    [
        "project",
        "file_category",
        "file_name",
    ]
)

cand.to_csv(
    INV_OUT / "CANDIDATE_FILES_FOR_DOWNLOAD.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 100)
print("CANDIDATE NON-RAW FILES ACROSS ALL THREE DATASETS")
print("=" * 100)

print(
    cand[
        [
            "project",
            "file_name",
            "file_category",
            "file_size_MB",
        ]
    ].to_string(
        index=False
    )
)

print()
print("Output:", OUT)
