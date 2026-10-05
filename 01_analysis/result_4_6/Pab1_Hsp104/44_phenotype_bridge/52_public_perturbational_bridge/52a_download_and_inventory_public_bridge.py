#!/usr/bin/env python3

from pathlib import Path
from urllib.parse import urljoin
import hashlib
import json
import re
import zipfile

import pandas as pd
import requests
from bs4 import BeautifulSoup


ROOT = Path(
    "/"
)

OUT = (
    ROOT /
    "44_phenotype_bridge" /
    "52_public_perturbational_bridge"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                
                                                              

SKY1_ARTICLE = (
    "https://www.nature.com/articles/"
    "s41467-019-11550-w"
)

GUK1_ARTICLE = (
    "https://journals.plos.org/plosgenetics/"
    "article?id=10.1371/journal.pgen.1008951"
)


session = requests.Session()

session.headers.update({
    "User-Agent": (
        "Mozilla/5.0 "
        "(X11; Linux x86_64) "
        "AppleWebKit/537.36 "
        "Chrome/140 Safari/537.36"
    )
})


                                                              
         
                                                              

def sha256(path):

    h = hashlib.sha256()

    with open(
        path,
        "rb",
    ) as fh:

        while True:

            b = fh.read(
                1024 * 1024
            )

            if not b:
                break

            h.update(
                b
            )

    return h.hexdigest()


def valid_xlsx(path):

    if (
        not path.exists()
        or path.stat().st_size < 5000
    ):
        return False

    if not zipfile.is_zipfile(
        path
    ):
        return False

    try:

        with zipfile.ZipFile(
            path
        ) as z:

            return (
                "[Content_Types].xml"
                in z.namelist()
            )

    except Exception:

        return False


def download(
    url,
    dest,
):

    print()
    print(
        "DOWNLOAD:",
        url,
    )

    r = session.get(
        url,
        timeout=180,
        allow_redirects=True,
    )

    print(
        "HTTP:",
        r.status_code,
    )

    print(
        "final URL:",
        r.url,
    )

    print(
        "bytes:",
        len(r.content),
    )

    r.raise_for_status()

    dest.write_bytes(
        r.content
    )

    if not valid_xlsx(
        dest
    ):

        head = (
            dest.read_bytes()[:300]
            if dest.exists()
            else b""
        )

        raise RuntimeError(
            f"Not a valid XLSX: {dest}\n"
            f"head={head!r}"
        )

    print(
        "VALID XLSX:",
        dest.name,
        f"{dest.stat().st_size / 1024 / 1024:.3f} MiB",
    )

    return r.url


def text_of_row(
    row,
):

    return " | ".join(
        ""
        if pd.isna(v)
        else str(v)
        for v in row
    )


def workbook_inventory(
    path,
    regex,
    max_hits_per_sheet=25,
):

    print()
    print(
        "=" * 120
    )

    print(
        "WORKBOOK:",
        path.name
    )

    print(
        "=" * 120
    )

    xl = pd.ExcelFile(
        path,
        engine="openpyxl",
    )

    print(
        "sheets:",
        xl.sheet_names,
    )

    rows = []

    for sheet in xl.sheet_names:

                                    
                                                                      
        x = pd.read_excel(
            path,
            sheet_name=sheet,
            header=None,
            engine="openpyxl",
        )

        print()
        print(
            f"[{sheet}] shape={x.shape}"
        )

        mask = (
            x.astype(str)
             .apply(
                 lambda c:
                     c.str.contains(
                         regex,
                         regex=True,
                         case=False,
                         na=False,
                     )
             )
             .any(
                 axis=1
             )
        )

        idx = list(
            x.index[
                mask
            ]
        )

        print(
            "relevant rows:",
            len(idx),
        )

        for i in idx[
            :max_hits_per_sheet
        ]:

            lo = max(
                0,
                i - 1,
            )

            hi = min(
                len(x),
                i + 2,
            )

            for j in range(
                lo,
                hi,
            ):

                rows.append({
                    "file":
                        path.name,

                    "sheet":
                        sheet,

                    "row_0based":
                        int(j),

                    "matched_row_0based":
                        int(i),

                    "text":
                        text_of_row(
                            x.iloc[j]
                        ),
                })

                print(
                    f"row {j}:",
                    text_of_row(
                        x.iloc[j]
                    )[:1000],
                )

            print(
                "-" * 80
            )

    return pd.DataFrame(
        rows
    )


                                                              
                                           
                                                              

print(
    "=" * 120
)

print(
    "RESOLVE SHATTUCK 2019 SOURCE DATA"
)

print(
    "=" * 120
)


r = session.get(
    SKY1_ARTICLE,
    timeout=60,
)

print(
    "article HTTP:",
    r.status_code,
)

r.raise_for_status()

soup = BeautifulSoup(
    r.text,
    "html.parser",
)


xlsx_candidates = []

for a in soup.find_all(
    "a",
    href=True,
):

    href = a.get(
        "href",
        ""
    )

    text = " ".join(
        a.stripped_strings
    )

    parent_text = " ".join(
        a.parent.stripped_strings
    ) if a.parent else ""

    if (
        ".xlsx" in href.lower()
        or
        "source data" in text.lower()
        or
        "source data" in parent_text.lower()
    ):

        xlsx_candidates.append(
            (
                text,
                parent_text[:500],
                urljoin(
                    SKY1_ARTICLE,
                    href,
                ),
            )
        )


print(
    "candidate links:"
)

for x in xlsx_candidates:

    print(
        x
    )


                                                             
source_candidates = [
    u
    for txt, parent, u
    in xlsx_candidates
    if (
        ".xlsx" in u.lower()
        and
        (
            "source data"
            in (
                txt + " " + parent
            ).lower()
            or
            "source"
            in u.lower()
        )
    )
]


                                                               
                                                  
if not source_candidates:

    source_candidates = [
        u
        for txt, parent, u
        in xlsx_candidates
        if ".xlsx" in u.lower()
    ]


source_candidates = list(
    dict.fromkeys(
        source_candidates
    )
)


if not source_candidates:

    raise RuntimeError(
        "Could not resolve Shattuck Source Data XLSX "
        "from Nature article HTML."
    )


SHATTUCK = (
    OUT /
    "Shattuck2019_Source_Data.xlsx"
)


downloaded = False

errors = []

for url in source_candidates:

    try:

        final_url = download(
            url,
            SHATTUCK,
        )

        downloaded = True
        break

    except Exception as e:

        errors.append(
            (
                url,
                repr(e),
            )
        )


if not downloaded:

    raise RuntimeError(
        "All Shattuck XLSX candidates failed:\n"
        +
        "\n".join(
            f"{u}: {e}"
            for u, e in errors
        )
    )


                                                              
                                               
                                                              

print()
print(
    "=" * 120
)

print(
    "RESOLVE ANDERSSON 2021 GUK1 SUPPLEMENTARY XLSX"
)

print(
    "=" * 120
)


r = session.get(
    GUK1_ARTICLE,
    timeout=60,
)

print(
    "article HTTP:",
    r.status_code,
)

r.raise_for_status()

soup = BeautifulSoup(
    r.text,
    "html.parser",
)


                                                            
links = []

for a in soup.find_all(
    "a",
    href=True,
):

    href = a[
        "href"
    ]

    txt = " ".join(
        a.stripped_strings
    )

                                                           
                                                    
    context = []

    node = a

    for _ in range(
        4
    ):

        if node is None:
            break

        context.append(
            " ".join(
                node.stripped_strings
            )
        )

        node = node.parent

    context = " ".join(
        context
    )

    links.append({
        "text":
            txt,

        "context":
            context[:3000],

        "url":
            urljoin(
                GUK1_ARTICLE,
                href,
            ),
    })


targets = {
    "S4":
        OUT /
        "Andersson2021_S4_Guk1.xlsx",

    "S6":
        OUT /
        "Andersson2021_S6_Guk1_HSP104.xlsx",

    "S10":
        OUT /
        "Andersson2021_S10_Guk1_timelapse.xlsx",
}


resolved = {}


for label in targets:

    rx = re.compile(
        rf"\b{label}\s+Table\b",
        flags=re.I,
    )

    cand = [
        q
        for q in links
        if rx.search(
            q[
                "context"
            ]
        )
        and (
            "xlsx"
            in (
                q["text"]
                +
                " "
                +
                q["context"]
                +
                " "
                +
                q["url"]
            ).lower()
            or
            "supplementary"
            in q["url"].lower()
            or
            "file"
            in q["url"].lower()
        )
    ]


    print()
    print(
        label,
        "candidates:",
        len(
            cand
        ),
    )

    for q in cand[
        :10
    ]:

        print(
            q[
                "text"
            ],
            q[
                "url"
            ],
        )


    success = False

    seen = set()

    for q in cand:

        url = q[
            "url"
        ]

        if url in seen:
            continue

        seen.add(
            url
        )

        try:

            final_url = download(
                url,
                targets[
                    label
                ],
            )

            resolved[
                label
            ] = final_url

            success = True
            break

        except Exception:

                                                    
                                                   
            continue


    if not success:

        raise RuntimeError(
            f"Could not resolve/download {label} Table "
            "from PLOS article."
        )


                                                              
            
                                                              

manifest = []


for label, path in [
    (
        "Shattuck2019_Source_Data",
        SHATTUCK,
    ),
    (
        "Andersson2021_S4",
        targets["S4"],
    ),
    (
        "Andersson2021_S6",
        targets["S6"],
    ),
    (
        "Andersson2021_S10",
        targets["S10"],
    ),
]:

    manifest.append({
        "label":
            label,

        "path":
            str(
                path
            ),

        "bytes":
            path.stat().st_size,

        "sha256":
            sha256(
                path
            ),
    })


manifest = pd.DataFrame(
    manifest
)

manifest.to_csv(
    OUT /
    "52a_PUBLIC_DATA_MANIFEST.tsv",
    sep="\t",
    index=False,
)


                                                              
                                     
                                                              

RX_SG = (
    r"HSP104|hsp104|Pab1|"
    r"Fig\.?\s*6|Figure\s*6|"
    r"recovery|heat"
)

RX_GUK = (
    r"HSP104|hsp104|Guk1|guk1|"
    r"SSA1|SSA2|SSA4|"
    r"recovery|heat|focus|foci|"
    r"time"
)


inv = []


inv.append(
    workbook_inventory(
        SHATTUCK,
        RX_SG,
        max_hits_per_sheet=40,
    )
)


for label in [
    "S4",
    "S6",
    "S10",
]:

    inv.append(
        workbook_inventory(
            targets[
                label
            ],
            RX_GUK,
            max_hits_per_sheet=40,
        )
    )


inv = pd.concat(
    inv,
    ignore_index=True,
)


inv.to_csv(
    OUT /
    "52a_RELEVANT_WORKBOOK_ROWS.tsv",
    sep="\t",
    index=False,
)


                                                              
                                         
                                                              

print()
print(
    "=" * 120
)

print(
    "LOCAL FROZEN GUK1 ARTIFACTS"
)

print(
    "=" * 120
)


local_hits = []


PHENO = (
    ROOT /
    "44_phenotype_bridge"
)


for p in PHENO.rglob(
    "*"
):

    if not p.is_file():
        continue

    name = p.name.lower()

    if (
        "46n" in name
        or
        (
            "guk1"
            in name
            and
            (
                "phenotype"
                in name
                or
                "auc"
                in name
                or
                "recovery"
                in name
            )
        )
    ):

        local_hits.append(
            p
        )


for p in sorted(
    local_hits
):

    print(
        p
    )


pd.DataFrame({
    "path":
        [
            str(p)
            for p
            in sorted(
                local_hits
            )
        ]
}).to_csv(
    OUT /
    "52a_LOCAL_GUK1_ARTIFACT_INVENTORY.tsv",
    sep="\t",
    index=False,
)


                                                              
         
                                                              

print()
print(
    "=" * 120
)

print(
    "PUBLIC DATA MANIFEST"
)

print(
    "=" * 120
)

print(
    manifest.to_string(
        index=False
    )
)


print()
print(
    "=" * 120
)

print(
    "WRITTEN"
)

print(
    "=" * 120
)

for f in [
    "52a_PUBLIC_DATA_MANIFEST.tsv",
    "52a_RELEVANT_WORKBOOK_ROWS.tsv",
    "52a_LOCAL_GUK1_ARTIFACT_INVENTORY.tsv",
]:

    print(
        OUT / f
    )
