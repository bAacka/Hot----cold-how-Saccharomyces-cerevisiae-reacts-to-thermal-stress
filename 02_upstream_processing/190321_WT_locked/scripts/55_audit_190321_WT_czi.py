#!/usr/bin/env python3

from pathlib import Path
import xml.etree.ElementTree as ET

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
        f"Expected exactly one WT CZI; found {len(files)}:\n"
        + "\n".join(str(x) for x in files)
    )

path = files[0]

print("CZI:")
print(path)
print("bytes:", path.stat().st_size)


czi = CziFile(
    path
)

print()
print("=" * 100)
print("CORE STRUCTURE")
print("=" * 100)

print("dims:", czi.dims)
print("size:", czi.size)
print("pixel_type:", czi.pixel_type)
print("shape_is_consistent:", czi.shape_is_consistent)
print("is_mosaic:", czi.is_mosaic())

dims_shape = czi.get_dims_shape()

print()
print("get_dims_shape():")

for i, x in enumerate(dims_shape):
    print(f"[{i}]", x)


                                                              
                   
                                                              

meta = czi.meta

xml_text = ET.tostring(
    meta,
    encoding="unicode",
)

(
    OUT
    / "190321_WT_METADATA.xml"
).write_text(
    xml_text
)


def lname(tag):
    return str(tag).split("}")[-1]


                                                              
          
                                                              

channels = []

for elem in meta.iter():

    if lname(elem.tag) != "Channel":
        continue

    rec = {
        "Id": elem.attrib.get("Id", ""),
        "Name": elem.attrib.get("Name", ""),
    }

    for child in elem:

        text = (
            child.text.strip()
            if child.text
            else ""
        )

        if text:
            rec[lname(child.tag)] = text

    channels.append(rec)


print()
print("=" * 100)
print("CHANNELS")
print("=" * 100)

if channels:

    ch = (
        pd.DataFrame(channels)
        .drop_duplicates()
    )

    print(
        ch.to_string(index=False)
    )

    ch.to_csv(
        OUT / "CHANNEL_METADATA.tsv",
        sep="\t",
        index=False,
    )


                                                              
                  
                                                              

scales = []

for elem in meta.iter():

    if lname(elem.tag) != "Distance":
        continue

    dim = elem.attrib.get(
        "Id",
        ""
    )

    if dim not in {"X", "Y", "Z"}:
        continue

    value = None
    unit = None

    for child in elem:

        tag = lname(child.tag)

        if tag == "Value":
            value = child.text

        elif tag in {
            "DefaultUnitFormat",
            "Unit",
        }:
            unit = child.text

    scales.append({
        "dimension": dim,
        "raw_value": value,
        "metadata_unit": unit,
    })


sc = pd.DataFrame(scales).drop_duplicates()

print()
print("=" * 100)
print("PHYSICAL SCALING")
print("=" * 100)

print(
    sc.to_string(index=False)
)

sc.to_csv(
    OUT / "PHYSICAL_SCALING.tsv",
    sep="\t",
    index=False,
)


                                                              
                            
 
                                                            
                                                              

timing_names = {
    "StartTime",
    "TimeSpan",
    "Interval",
    "Increment",
    "Duration",
    "TimeIncrement",
    "TimeInterval",
}


timing_rows = []

for elem in meta.iter():

    name = lname(
        elem.tag
    )

    if name not in timing_names:
        continue

    text = (
        elem.text.strip()
        if elem.text
        else ""
    )

    timing_rows.append({
        "tag":
            name,

        "text":
            text,

        "attributes":
            repr(
                dict(
                    elem.attrib
                )
            ),
    })


timing = pd.DataFrame(
    timing_rows
)


print()
print("=" * 100)
print("TIMING-LIKE XML ELEMENTS")
print("=" * 100)

if len(timing):

    print(
        timing.to_string(
            index=False
        )
    )

    timing.to_csv(
        OUT / "TIMING_METADATA_CANDIDATES.tsv",
        sep="\t",
        index=False,
    )

else:

    print(
        "No timing-like elements found in main XML metadata."
    )


print()
print("=" * 100)
print("DONE")
print("=" * 100)

print("Output:", OUT)
