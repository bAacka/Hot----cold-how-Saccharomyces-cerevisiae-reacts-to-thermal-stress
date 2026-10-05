#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime
import inspect
import re

import numpy as np
import pandas as pd

from aicspylibczi import CziFile


ROOT = Path(
    "/"
)

DATA = (
    ROOT
    / "44_phenotype_bridge"
    / "downloads"
    / "190321_WT"
)

OUT = (
    ROOT
    / "44_phenotype_bridge"
    / "dataset_audit"
    / "190321_WT"
    / "timestamps"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


files = sorted(
    DATA.glob("*.czi")
)

if len(files) != 1:
    raise RuntimeError(
        f"Expected exactly one CZI; found {len(files)}"
    )

path = files[0]

czi = CziFile(
    path
)


print("CZI:")
print(path)

print()
print("read_subblock_metadata signature:")
print(
    inspect.signature(
        czi.read_subblock_metadata
    )
)


                                                              
                         
 
                                                              
                                 
                                                              

records = []

for s in range(4):

    for t in range(16):

        try:

            result = czi.read_subblock_metadata(
                S=s,
                T=t,
                C=0,
                Z=5,
            )

        except TypeError:

                                                                         
            result = czi.read_subblock_metadata()


                                                              
                                                            
                                                              

        items = []

        if isinstance(
            result,
            list,
        ):
            items = result

        else:
            items = [result]


        found_strings = []

        def collect_strings(x):

            if x is None:
                return

            if isinstance(
                x,
                str,
            ):
                found_strings.append(
                    x
                )
                return

            if isinstance(
                x,
                bytes,
            ):
                found_strings.append(
                    x.decode(
                        "utf-8",
                        errors="ignore",
                    )
                )
                return

            if isinstance(
                x,
                dict,
            ):

                for k, v in x.items():

                    found_strings.append(
                        str(k)
                    )

                    collect_strings(v)

                return

            if isinstance(
                x,
                (tuple, list),
            ):

                for v in x:
                    collect_strings(v)

                return

            found_strings.append(
                str(x)
            )


        for item in items:
            collect_strings(
                item
            )


        text = "\n".join(
            found_strings
        )


                                                              
                                       
                                                              

        iso_hits = re.findall(
            r"20\d{2}-\d{2}-\d{2}T"
            r"\d{2}:\d{2}:\d{2}"
            r"(?:\.\d+)?Z?",
            text,
        )


                                                              
                                                   
                                                             
                                                              

        numeric_hits = []

        patterns = [
            r'AcquisitionTime[^0-9+\-]*([+\-]?\d+(?:\.\d+)?)',
            r'TimeStamp[^0-9+\-]*([+\-]?\d+(?:\.\d+)?)',
            r'Timestamp[^0-9+\-]*([+\-]?\d+(?:\.\d+)?)',
        ]

        for pattern in patterns:

            for hit in re.findall(
                pattern,
                text,
                flags=re.I,
            ):

                try:
                    numeric_hits.append(
                        float(hit)
                    )
                except ValueError:
                    pass


        records.append({
            "scene":
                s,

            "time_index":
                t,

            "n_iso_hits":
                len(iso_hits),

            "iso_hits":
                "|".join(
                    sorted(
                        set(
                            iso_hits
                        )
                    )
                ),

            "numeric_time_hits":
                "|".join(
                    str(x)
                    for x in sorted(
                        set(
                            numeric_hits
                        )
                    )
                ),

            "metadata_excerpt":
                text[:1000].replace(
                    "\n",
                    " ",
                ),
        })


df = pd.DataFrame(
    records
)

df.to_csv(
    OUT
    / "SUBBLOCK_TIMESTAMP_RAW.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                           
                                                              

parsed = []

for rec in records:

    iso = (
        rec["iso_hits"].split("|")
        if rec["iso_hits"]
        else []
    )

    dt = None

    if iso:

        s = iso[0]

        if s.endswith("Z"):
            s = (
                s[:-1]
                + "+00:00"
            )

        try:
            dt = datetime.fromisoformat(
                s
            )
        except ValueError:
            dt = None


    parsed.append({
        "scene":
            rec["scene"],

        "time_index":
            rec["time_index"],

        "timestamp":
            (
                dt.isoformat()
                if dt is not None
                else ""
            ),
    })


tdf = pd.DataFrame(
    parsed
)


                                                              
                                                              
                                                              

elapsed_rows = []

for scene, g in tdf.groupby(
    "scene"
):

    g = g.copy()

    valid = g[
        g[
            "timestamp"
        ] != ""
    ].copy()

    if len(valid):

        valid[
            "_dt"
        ] = pd.to_datetime(
            valid[
                "timestamp"
            ],
            utc=True,
        )

        t0 = valid[
            "_dt"
        ].min()

        valid[
            "elapsed_min"
        ] = (
            valid[
                "_dt"
            ]
            - t0
        ).dt.total_seconds() / 60.0

        lookup = dict(
            zip(
                valid[
                    "time_index"
                ],
                valid[
                    "elapsed_min"
                ],
            )
        )

    else:

        lookup = {}


    for t in range(16):

        elapsed_rows.append({
            "scene":
                scene,

            "time_index":
                t,

            "timestamp":
                g.loc[
                    g[
                        "time_index"
                    ] == t,
                    "timestamp",
                ].iloc[0],

            "elapsed_min_from_scene_first":
                lookup.get(
                    t,
                    np.nan,
                ),
        })


elapsed = pd.DataFrame(
    elapsed_rows
)


elapsed.to_csv(
    OUT
    / "FRAME_TIMESTAMPS.tsv",
    sep="\t",
    index=False,
)


                                                              
                                                   
 
               
                            
                                                              

INCREMENT_SEC = 599.300278

nominal = pd.DataFrame({
    "time_index":
        np.arange(
            16,
            dtype=int,
        ),
})

nominal[
    "elapsed_sec_nominal"
] = (
    nominal[
        "time_index"
    ]
    * INCREMENT_SEC
)

nominal[
    "elapsed_min_nominal"
] = (
    nominal[
        "elapsed_sec_nominal"
    ]
    / 60.0
)


                                                      
nominal[
    "biological_min_if_T0_is_0"
] = nominal[
    "elapsed_min_nominal"
]


                                                              
nominal[
    "authors_nominal_min_if_T1_is_10"
] = (
    nominal[
        "time_index"
    ]
    * 10
)


nominal.to_csv(
    OUT
    / "NOMINAL_TIME_MAPPING.tsv",
    sep="\t",
    index=False,
)


print()
print("=" * 100)
print("TIMESTAMP EXTRACTION")
print("=" * 100)

print(
    elapsed.to_string(
        index=False
    )
)


print()
print("=" * 100)
print("NOMINAL MAPPING")
print("=" * 100)

print(
    nominal.to_string(
        index=False
    )
)


print()
print(
    "Output:",
    OUT
)
