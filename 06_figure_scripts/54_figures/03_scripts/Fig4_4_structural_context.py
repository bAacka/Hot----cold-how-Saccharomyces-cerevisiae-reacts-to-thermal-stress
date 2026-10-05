#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Patch, Circle

from Bio.PDB import PDBParser, MMCIFParser


                                                              
       
                                                              

ROOT = Path("/")

SEARCH_DIRS = [
    ROOT / "38_pka_mechanistic_resources",
    ROOT / "39_mechanistic_site_context",
    ROOT / "40_structural_tractability",
    ROOT / "41_cdc19_s22_structure",
    ROOT / "42_disordered_phosphosite_mechanisms",
    ROOT / "43_structural_mechanism_freeze",
]

FIG = ROOT / "54_figures"

WORK = (
    FIG
    / "02_curated_inputs"
    / "Fig4_4"
)

REND = (
    FIG
    / "04_panels"
    / "Fig4_4_structures"
)

FINAL = FIG / "05_final"
LEGENDS = FIG / "06_legends"

for d in [WORK, REND, FINAL, LEGENDS]:
    d.mkdir(parents=True, exist_ok=True)


                                                              
                                
                                                              

SITES = {
    "ATG1": {
        "site": 515,
        "label": "Atg1-S515",
        "expected_site_plddt": 29.36,
        "expected_local_median": 28.77,
    },
    "MAF1": {
        "site": 90,
        "label": "Maf1-S90",
        "expected_site_plddt": 50.22,
        "expected_local_median": 50.22,
    },
    "MSN4": {
        "site": 316,
        "label": "Msn4-S316",
        "expected_site_plddt": 35.72,
        "expected_local_median": 34.84,
    },
}

CDC19_SITE = 22

CDC19_SASA_RANGE = (6.9, 7.8)
CDC19_RSA_RANGE = (0.044, 0.050)
CDC19_NEIGHBORS_5A = 9


                                                              
         
                                                              

def fail(msg):
    raise RuntimeError(msg)


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def norm(x):
    return re.sub(
        r"[^a-z0-9]+",
        "_",
        str(x).lower(),
    ).strip("_")


def download(url, dest):
    dest = Path(dest)

    print("DOWNLOAD:", url)

    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0"
        },
    )

    with urllib.request.urlopen(
        req,
        timeout=60,
    ) as r:
        dest.write_bytes(
            r.read()
        )

    return dest


def table_files():

    out = []

    for d in SEARCH_DIRS:
        if not d.exists():
            continue

        for p in d.rglob("*"):

            if not p.is_file():
                continue

            name = p.name.lower()

            if name.endswith(
                (
                    ".tsv",
                    ".tsv.gz",
                    ".csv",
                    ".csv.gz",
                )
            ):
                out.append(p)

    return sorted(out)


def read_table(path):

    if ".tsv" in path.name.lower():
        return pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )

    return pd.read_csv(
        path,
        low_memory=False,
    )


                                                              
                            
                                                              

def local_uniprot_for_gene(gene):

    candidates = []

    for p in table_files():

        try:
            df = read_table(p)
        except Exception:
            continue

        if len(df) == 0:
            continue

                                                  
        rowmask = pd.Series(
            False,
            index=df.index,
        )

        for c in df.columns:

            if df[c].dtype != object:
                continue

            s = (
                df[c]
                .astype(str)
                .str.upper()
            )

            rowmask |= (
                s == gene.upper()
            )

        sub = df[rowmask]

        if len(sub) == 0:
            continue

        acc_cols = [
            c for c in df.columns
            if (
                "uniprot" in norm(c)
                or "accession" in norm(c)
            )
        ]

        for c in acc_cols:

            for value in (
                sub[c]
                .dropna()
                .astype(str)
            ):

                value = value.strip()

                if re.fullmatch(
                    r"[A-Z0-9]{6,10}",
                    value,
                ):
                    candidates.append(
                        (
                            p,
                            c,
                            value,
                        )
                    )

    accessions = sorted(
        set(
            x[2]
            for x in candidates
        )
    )

    if len(accessions) == 1:
        return accessions[0]

    if len(accessions) > 1:

        print(
            f"Multiple local UniProt candidates for {gene}:",
            accessions,
        )

                                                         
        six = [
            x for x in accessions
            if len(x) == 6
        ]

        if len(six) == 1:
            return six[0]

    return None


def uniprot_rest_for_gene(gene):

    query = (
        f"(gene_exact:{gene}) "
        f"AND (organism_id:559292)"
    )

    url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        + urllib.parse.urlencode(
            {
                "query": query,
                "format": "tsv",
                "fields": "accession,gene_names",
                "size": 10,
            }
        )
    )

    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Mozilla/5.0"
        },
    )

    with urllib.request.urlopen(
        req,
        timeout=60,
    ) as r:
        txt = (
            r.read()
            .decode("utf-8")
        )

    lines = [
        x for x in txt.splitlines()
        if x.strip()
    ]

    if len(lines) < 2:
        fail(
            f"UniProt could not resolve {gene}"
        )

    df = pd.read_csv(
        pd.io.common.StringIO(txt),
        sep="\t",
    )

    if "Entry" not in df.columns:
        fail(
            f"Unexpected UniProt response for {gene}"
        )

    return str(
        df.iloc[0]["Entry"]
    )


def resolve_accession(gene):

    acc = local_uniprot_for_gene(
        gene
    )

    if acc:
        print(
            f"{gene}: local UniProt = {acc}"
        )
        return acc

    print(
        f"{gene}: no unique local accession; querying UniProt"
    )

    acc = uniprot_rest_for_gene(
        gene
    )

    print(
        f"{gene}: UniProt = {acc}"
    )

    return acc


accessions = {
    gene: resolve_accession(gene)
    for gene in SITES
}

                                                               
if (
    "ATG1" in accessions
    and accessions["ATG1"] != "P53104"
):
    print(
        "WARNING: resolved ATG1 accession differs "
        f"from expected P53104: {accessions[ATG1]}"
    )


                                                              
                                
                                                              

def find_local_structure(tokens):

    tokens = [
        t.lower()
        for t in tokens
    ]

    hits = []

    for d in SEARCH_DIRS:

        if not d.exists():
            continue

        for p in d.rglob("*"):

            if not p.is_file():
                continue

            if not p.name.lower().endswith(
                (
                    ".pdb",
                    ".cif",
                    ".mmcif",
                )
            ):
                continue

            name = p.name.lower()

            if any(
                token in name
                for token in tokens
            ):
                hits.append(p)

    if hits:
        return sorted(hits)[0]

    return None


def ensure_rcsb(pdb_id):

    pdb_id = pdb_id.upper()

    local = find_local_structure(
        [pdb_id.lower()]
    )

    if local:
        print(
            f"{pdb_id}: local structure = {local}"
        )
        return local

    dest = (
        WORK
        / f"{pdb_id}.pdb"
    )

    if not dest.exists():
        download(
            f"https://files.rcsb.org/download/{pdb_id}.pdb",
            dest,
        )

    return dest


pdb_1a3w = ensure_rcsb(
    "1A3W"
)

pdb_1a3x = ensure_rcsb(
    "1A3X"
)


                                                              
                                  
                                                              

def ensure_alphafold(
    gene,
    accession,
):

    local = find_local_structure(
        [
            gene.lower(),
            accession.lower(),
            f"af-{accession.lower()}",
        ]
    )

    if local:
        print(
            f"{gene}: local AlphaFold model = {local}"
        )
        return local

                                                              
    versions = [
        "v6",
        "v4",
        "v3",
        "v2",
    ]

    for version in versions:

        dest = (
            WORK
            / f"AF-{accession}-F1-model_{version}.pdb"
        )

        if dest.exists():
            return dest

        url = (
            "https://alphafold.ebi.ac.uk/files/"
            f"AF-{accession}-F1-model_{version}.pdb"
        )

        try:
            download(
                url,
                dest,
            )

            if (
                dest.exists()
                and dest.stat().st_size > 1000
            ):
                return dest

        except Exception as e:

            if dest.exists():
                dest.unlink()

            print(
                f"  {version} unavailable: {e}"
            )

    fail(
        f"Could not locate/download AlphaFold model "
        f"for {gene} ({accession})"
    )


af_models = {
    gene: ensure_alphafold(
        gene,
        accession,
    )
    for gene, accession
    in accessions.items()
}


                                                              
                   
                                                              

def parser_for(path):

    if path.suffix.lower() in {
        ".cif",
        ".mmcif",
    }:
        return MMCIFParser(
            QUIET=True
        )

    return PDBParser(
        QUIET=True
    )


def structure(path):

    return parser_for(
        path
    ).get_structure(
        path.stem,
        str(path),
    )


def chain_with_site(
    path,
    site,
):

    st = structure(path)

    hits = []

    for model in st:

        for chain in model:

            for residue in chain:

                het, resseq, icode = (
                    residue.id
                )

                if (
                    het == " "
                    and resseq == site
                ):
                    hits.append(
                        chain.id
                    )

        break

    if not hits:
        fail(
            f"Residue {site} not found in {path}"
        )

    return hits[0]


cdc_chain_1a3w = chain_with_site(
    pdb_1a3w,
    CDC19_SITE,
)

cdc_chain_1a3x = chain_with_site(
    pdb_1a3x,
    CDC19_SITE,
)


                                                              
                                  
                                                              

def plddt_profile(
    path,
    site,
    flank=50,
):

    st = structure(path)

    rows = []

                                        
    model = next(iter(st))

    for chain in model:

        chain_rows = []

        for residue in chain:

            het, resseq, icode = (
                residue.id
            )

            if het != " ":
                continue

            if "CA" not in residue:
                continue

            if not (
                site - flank
                <= resseq
                <= site + flank
            ):
                continue

            ca = residue["CA"]

            chain_rows.append(
                {
                    "residue": int(resseq),
                    "relative_position": int(
                        resseq - site
                    ),
                    "plddt": float(
                        ca.get_bfactor()
                    ),
                    "chain": chain.id,
                }
            )

        if any(
            x["residue"] == site
            for x in chain_rows
        ):
            rows = chain_rows
            break

    if not rows:
        fail(
            f"Could not extract ±{flank} pLDDT "
            f"around site {site} from {path}"
        )

    return pd.DataFrame(rows)


profiles = []

validation = []

for gene, spec in SITES.items():

    site = spec['site']

    df = plddt_profile(
        af_models[gene],
        site,
        flank=50,
    )

    df.insert(
        0,
        "gene",
        gene,
    )

    profiles.append(df)

    site_plddt = float(
        df.loc[
            df["residue"] == site,
            "plddt",
        ].iloc[0]
    )

    local21 = df[
        df["relative_position"]
        .between(
            -10,
            10,
        )
    ]

    local_median = float(
        local21["plddt"]
        .median()
    )

    validation.append(
        {
            "gene": gene,
            "site": site,
            "site_plddt": site_plddt,
            "local21_median": local_median,
            "expected_site_plddt":
                spec["expected_site_plddt"],
            "expected_local_median":
                spec["expected_local_median"],
        }
    )

                                
    if abs(
        site_plddt
        - spec["expected_site_plddt"]
    ) > 2.0:
        fail(
            f"{gene}: model pLDDT at site = "
            f"{site_plddt:.2f}, expected ≈ "
            f"{spec[expected_site_plddt]:.2f}. "
            "Likely wrong AlphaFold model/version."
        )


profiles = pd.concat(
    profiles,
    ignore_index=True,
)

validation = pd.DataFrame(
    validation
)

profiles.to_csv(
    WORK
    / "Fig4_4_local_pLDDT_profiles.tsv",
    sep="\t",
    index=False,
)

validation.to_csv(
    WORK
    / "Fig4_4_pLDDT_validation.tsv",
    sep="\t",
    index=False,
)


                                                              
                       
                                                              

manifest_rows = []

for role, path in [
    ("CDC19_1A3W", pdb_1a3w),
    ("CDC19_1A3X", pdb_1a3x),
]:

    manifest_rows.append(
        {
            "role": role,
            "path": str(path),
            "sha256": sha256(path),
        }
    )

for gene, path in af_models.items():

    manifest_rows.append(
        {
            "role":
                f"{gene}_AlphaFold",
            "path": str(path),
            "sha256": sha256(path),
            "uniprot": accessions[gene],
        }
    )

pd.DataFrame(
    manifest_rows
).to_csv(
    WORK
    / "Fig4_4_structure_manifest.tsv",
    sep="\t",
    index=False,
)


                                                              
       
                                                              

pymol = shutil.which(
    "pymol"
)

if pymol is None:
    fail(
        "PyMOL executable not found. Install pymol-open-source "
        "in the active environment, then rerun."
    )


                                                              
             
                                                              

COMMON_PML = """
bg_color white
set ray_opaque_background, off
set antialias, 2
set orthoscopic, on
set depth_cue, 0
set fog, 0
set ambient, 0.55
set direct, 0.45
set specular, 0.20
set shininess, 25
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set ray_trace_mode, 1
set ray_shadows, 0
"""


def run_pymol(
    text,
    name,
):

    pml = (
        REND
        / f"{name}.pml"
    )

    pml.write_text(
        text,
        encoding="utf-8",
    )

    subprocess.run(
        [
            pymol,
            "-cq",
            str(pml),
        ],
        check=True,
    )


                                                              
                           
                                                              

whole_png = (
    REND
    / "A_Cdc19_1A3W_whole.png"
)

pml = f"""
{COMMON_PML}

load {pdb_1a3w}, cdc19

remove solvent
hide everything

show cartoon, cdc19
color grey70, cdc19

show surface, cdc19
color grey85, cdc19
set transparency, 0.62, cdc19

select s22, cdc19 and chain {cdc_chain_1a3w} and resi 22 and resn SER

# Deliberately exaggerate the phosphosite for figure readability.
show sticks, s22
set stick_radius, 0.38, s22
color magenta, s22

# Two large reference spheres on the actual Ser22 residue.
select s22mark, s22 and (name CA or name OG)
show spheres, s22mark
set sphere_scale, 1.20, s22mark
color magenta, s22mark

orient cdc19
zoom cdc19, 4

viewport 1600, 1200
ray 1600, 1200
png {whole_png}, dpi=300

quit
"""

run_pymol(
    pml,
    "A_Cdc19_1A3W_whole",
)


                                                              
                                 
                                                              

def closeup_pml(
    structure_path,
    chain,
    out_png,
    object_name,
):

    return f"""
{COMMON_PML}

load {structure_path}, {object_name}

remove solvent
hide everything

select site, {object_name} and chain {chain} and resi 22 and resn SER
select local5, byres (site around 5)
select local8, byres (site around 8)

show cartoon, local8
color grey80, local8

show sticks, local5
set stick_radius, 0.15, local5
color grey55, local5

# Deliberately exaggerate the phosphosite.
show sticks, site
set stick_radius, 0.36, site
color magenta, site

select sitemark, site and (name CA or name OG)
show spheres, sitemark
set sphere_scale, 0.95, sitemark
color magenta, sitemark

orient site
zoom local8, 3

viewport 1400, 1100
ray 1400, 1100
png {out_png}, dpi=300

quit
"""


b1_png = (
    REND
    / "B1_Cdc19_S22_1A3W.png"
)

b2_png = (
    REND
    / "B2_Cdc19_S22_1A3X.png"
)

run_pymol(
    closeup_pml(
        pdb_1a3w,
        cdc_chain_1a3w,
        b1_png,
        "cdc1",
    ),
    "B1_Cdc19_S22_1A3W",
)

run_pymol(
    closeup_pml(
        pdb_1a3x,
        cdc_chain_1a3x,
        b2_png,
        "cdc2",
    ),
    "B2_Cdc19_S22_1A3X",
)


                                                              
                                  
                                                              

AF_COLORS = """
set_color af_vlow, [1.00, 0.49, 0.27]
set_color af_low,  [1.00, 0.86, 0.35]
set_color af_conf, [0.35, 0.78, 0.93]
set_color af_high, [0.10, 0.32, 0.75]
"""


def alphafold_pml(
    model_path,
    site,
    out_png,
    obj,
):

    lo = max(
        1,
        site - 45,
    )

    hi = site + 45

    return f"""
{COMMON_PML}
{AF_COLORS}

load {model_path}, {obj}

hide everything

select local, {obj} and resi {lo}-{hi}
select site, {obj} and resi {site}

show cartoon, local

select vlow, local and b < 50
select low, local and b >= 50 and b < 70
select conf, local and b >= 70 and b < 90
select high, local and b >= 90

color af_vlow, vlow
color af_low, low
color af_conf, conf
color af_high, high

# Deliberately exaggerate the phosphosite.
show sticks, site
set stick_radius, 0.34, site
color magenta, site

select sitemark, site and (name CA or name OG)
show spheres, sitemark
set sphere_scale, 0.90, sitemark
color magenta, sitemark

orient local
zoom local, 3

viewport 1300, 1000
ray 1300, 1000
png {out_png}, dpi=300

quit
"""


c_pngs = {}

for gene, spec in SITES.items():

    out = (
        REND
        / f"C_{gene}_{spec['site']}_AlphaFold.png"
    )

    run_pymol(
        alphafold_pml(
            af_models[gene],
            spec['site'],
            out,
            gene.lower(),
        ),
        f"C_{gene}_{spec['site']}_AlphaFold",
    )

    c_pngs[gene] = out


                                                              
                       
                                                              

plt.rcParams.update(
    {
        "font.family":
            "DejaVu Sans",
        "font.size": 9.0,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

fig = plt.figure(
    figsize=(18.0, 15.0),
)

outer = fig.add_gridspec(
    3,
    12,
    height_ratios=[
        1.25,
        0.95,
        0.80,
    ],
    hspace=0.25,
    wspace=0.22,
)


                                                              
   
                                                              

axA = fig.add_subplot(
    outer[0, 0:6]
)

img = mpimg.imread(
    whole_png
)

axA.imshow(img)
axA.axis("off")

                                                              
                                                            
 
                                                               
                                                             
                                                        
                                                              

rgb = np.asarray(img)[..., :3].astype(float)

if rgb.max() > 1.5:
    rgb = rgb / 255.0

r = rgb[..., 0]
g = rgb[..., 1]
b = rgb[..., 2]

magenta_mask = (
    (r > 0.55)
    & (b > 0.45)
    & (g < 0.55)
    & ((r - g) > 0.20)
    & ((b - g) > 0.12)
)

yy, xx = np.where(magenta_mask)

if len(xx) < 5:
    raise RuntimeError(
        "Could not locate rendered S22 magenta pixels "
        "in Panel A."
    )

                                                                 
                                                       
site_x_px = float(np.median(xx))
site_y_px = float(np.median(yy))

img_h, img_w = rgb.shape[:2]

site_x = site_x_px / (img_w - 1)
site_y = 1.0 - site_y_px / (img_h - 1)

                                                       
halo = Circle(
    (site_x, site_y),
    radius=0.052,
    transform=axA.transAxes,
    fill=False,
    edgecolor="#B43A77",
    linewidth=2.2,
    zorder=10,
)

axA.add_patch(halo)

                                                                
                                                                 
text_x = max(site_x - 0.20, 0.10)
text_y = min(site_y + 0.13, 0.86)

axA.annotate(
    "S22",
    xy=(site_x, site_y),
    xycoords=axA.transAxes,
    xytext=(text_x, text_y),
    textcoords=axA.transAxes,
    ha="right",
    va="center",
    fontsize=11.0,
    fontweight="bold",
    color="#B43A77",
    arrowprops=dict(
        arrowstyle="-",
        color="#B43A77",
        linewidth=1.8,
        shrinkA=5,
        shrinkB=18,
        connectionstyle="arc3,rad=0.0",
    ),
    zorder=11,
)

axA.text(
    0.01,
    0.98,
    "A",
    transform=axA.transAxes,
    fontsize=16,
    fontweight="bold",
    va="top",
)

axA.text(
    0.50,
    0.015,
    "Cdc19 • 1A3W",
    transform=axA.transAxes,
    ha="center",
    va="bottom",
    fontsize=9.5,
)


                                                              
   
                                                              

Bgrid = outer[
    0,
    6:12
].subgridspec(
    1,
    2,
    wspace=0.04,
)

for i, (
    img_path,
    pdb_label,
) in enumerate(
    [
        (
            b1_png,
            "1A3W",
        ),
        (
            b2_png,
            "1A3X",
        ),
    ]
):

    ax = fig.add_subplot(
        Bgrid[0, i]
    )

    ax.imshow(
        mpimg.imread(
            img_path
        )
    )

    ax.axis("off")

    ax.text(
        0.50,
        0.015,
        pdb_label,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.5,
    )

    if i == 0:
        ax.text(
            -0.03,
            0.98,
            "B",
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )


                                                              
   
                                                              

Cgrid = outer[
    1,
    :
].subgridspec(
    1,
    3,
    wspace=0.04,
)

for i, gene in enumerate(
    [
        "ATG1",
        "MAF1",
        "MSN4",
    ]
):

    ax = fig.add_subplot(
        Cgrid[0, i]
    )

    ax.imshow(
        mpimg.imread(
            c_pngs[gene]
        )
    )

    ax.axis("off")

    ax.text(
        0.50,
        0.015,
        SITES[gene]["label"],
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=9.5,
    )

    if i == 0:
        ax.text(
            -0.025,
            0.98,
            "C",
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )


                                                              
                                        
                                                              

plddt_handles = [
    Patch(
        facecolor="#FF7D45",
        edgecolor="none",
        label="pLDDT < 50",
    ),
    Patch(
        facecolor="#FFDB59",
        edgecolor="none",
        label="50–70",
    ),
    Patch(
        facecolor="#59C7ED",
        edgecolor="none",
        label="70–90",
    ),
    Patch(
        facecolor="#1A52BF",
        edgecolor="none",
        label="> 90",
    ),
]

legend_ax = fig.add_axes(
    [0.31, 0.341, 0.38, 0.025]
)
legend_ax.axis("off")

legend_ax.legend(
    handles=plddt_handles,
    loc="center",
    ncol=4,
    frameon=False,
    fontsize=8.2,
    handlelength=1.35,
    handleheight=0.9,
    columnspacing=1.6,
    handletextpad=0.45,
)


                                                              
                    
                                                              

Dgrid = outer[
    2,
    :
].subgridspec(
    1,
    3,
    wspace=0.18,
)

for i, gene in enumerate(
    [
        "ATG1",
        "MAF1",
        "MSN4",
    ]
):

    ax = fig.add_subplot(
        Dgrid[0, i]
    )

    d = profiles[
        profiles["gene"]
        == gene
    ].copy()

    ax.plot(
        d["relative_position"],
        d["plddt"],
        lw=1.6,
        color="#555555",
    )

    ax.axvline(
        0,
        color="#B43A77",
        lw=1.5,
    )

    ax.axhline(
        50,
        color="#AAAAAA",
        ls="--",
        lw=0.7,
    )

    ax.axhline(
        70,
        color="#AAAAAA",
        ls=":",
        lw=0.7,
    )

    ax.scatter(
        [0],
        [
            float(
                d.loc[
                    d["relative_position"] == 0,
                    "plddt",
                ].iloc[0]
            )
        ],
        s=42,
        color="#B43A77",
        zorder=3,
    )

    ax.set_xlim(
        -50,
        50,
    )

    ax.set_ylim(
        0,
        100,
    )

    ax.set_xticks(
        [
            -50,
            -25,
            0,
            25,
            50,
        ]
    )

    ax.set_xlabel(
        "Позиция спрямо фосфосайта"
    )

    if i == 0:

        ax.set_ylabel(
            "pLDDT"
        )

        ax.text(
            -0.16,
            1.04,
            "D",
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="top",
        )

    else:

        ax.tick_params(
            axis="y",
            labelleft=False,
        )

    ax.text(
        0.03,
        0.94,
        SITES[gene]["label"],
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.0,
    )

    ax.spines["top"].set_visible(
        False
    )

    ax.spines["right"].set_visible(
        False
    )


                   

fig.subplots_adjust(
    left=0.04,
    right=0.985,
    top=0.985,
    bottom=0.055,
)


                                                              
        
                                                              

OUTBASE = (
    FINAL
    / "Fig4_4_structural_context"
)

fig.savefig(
    str(OUTBASE) + ".png",
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE) + ".pdf",
    bbox_inches="tight",
)

fig.savefig(
    str(OUTBASE) + ".svg",
    bbox_inches="tight",
)

plt.close(fig)


                                                              
         
                                                              

val = (
    validation
    .set_index("gene")
)

legend = f"""Фигура 4.4. Структурен контекст на приоритетни фосфорилационни позиции.

(A) Експерименталната кристалографска структура 1A3W на Cdc19. Cdc19-S22 е означен в магента върху полупрозрачната молекулна повърхност.

(B) Локалното структурно обкръжение на Cdc19-S22 в независимите кристалографски структури 1A3W и 1A3X. В анализираните структури остатъкът показва ниска solvent accessibility (SASA приблизително {CDC19_SASA_RANGE[0]:.1f}–{CDC19_SASA_RANGE[1]:.1f} Å²; RSA {CDC19_RSA_RANGE[0]:.3f}–{CDC19_RSA_RANGE[1]:.3f}) и {CDC19_NEIGHBORS_5A} съседни остатъка в радиус до 5 Å. Структурният резултат подкрепя хипотеза за състояние-зависима достъпност, но не демонстрира директно конформационен преход.

(C) Локален AlphaFold структурен контекст около Atg1-S515, Maf1-S90 и Msn4-S316. Моделите са оцветени по pLDDT, а анализираният остатък е маркиран в магента. Ниската увереност на локалните модели не трябва да се интерпретира като експериментално определена третична структура.

(D) AlphaFold pLDDT профили в ±50-аминокиселинен прозорец около трите фосфосайта. pLDDT на самите позиции е {val.loc["ATG1","site_plddt"]:.2f} за Atg1-S515, {val.loc["MAF1","site_plddt"]:.2f} за Maf1-S90 и {val.loc["MSN4","site_plddt"]:.2f} за Msn4-S316. Прекъснатите хоризонтални линии при pLDDT 50 и 70 са визуални confidence ориентири и не представляват директна мярка за intrinsic disorder. Maf1-S90 и Msn4-S316 са допълнително анотирани като разположени в неструктурирани региони в използваните функционални анотации.

Фигурата илюстрира два различни физични регулаторни контекста: слабо достъпен сайт в компактна експериментална структура при Cdc19-S22 и гъвкави/нискоуверени регулаторни региони при Atg1-S515, Maf1-S90 и Msn4-S316. Данните не подкрепят един универсален структурен механизъм на температурно-зависима фосфорегулация.
"""

legend_path = (
    LEGENDS
    / "Fig4_4_legend_bg.txt"
)

legend_path.write_text(
    legend,
    encoding="utf-8",
)


                                                              
                
                                                              

print()
print("===== FIGURE 4.4 COMPLETE =====")

for gene in [
    "ATG1",
    "MAF1",
    "MSN4",
]:

    row = val.loc[gene]

    print(
        f"{gene}: "
        f"UniProt={accessions[gene]} | "
        f"site pLDDT={row.site_plddt:.2f} | "
        f"local21 median={row.local21_median:.2f}"
    )

print()
print(str(OUTBASE) + ".pdf")
print(legend_path)
