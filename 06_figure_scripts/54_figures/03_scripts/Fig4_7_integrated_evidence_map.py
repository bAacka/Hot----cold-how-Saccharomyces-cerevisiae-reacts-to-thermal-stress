#!/usr/bin/env python3

from pathlib import Path
import hashlib

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


                                                              
       
                                                              

ROOT = Path("/")

OUTDIR = ROOT / "54_figures/05_final"
LEGENDDIR = ROOT / "54_figures/06_legends"

OUTDIR.mkdir(parents=True, exist_ok=True)
LEGENDDIR.mkdir(parents=True, exist_ok=True)

BASENAME = "Fig4_7_integrated_evidence_map"


                                                              
              
                                                              

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.0,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

TEXT = "#222222"
MUTED = "#666666"
LINE = "#555555"
LIGHT_LINE = "#CFCFCF"

HEAT = "#A85A5A"
COLD = "#557AA8"
PKA = "#7C5A91"

STRUCT = "#66717B"
YEF3 = "#8B7150"
HSP104 = "#4D7463"

BOX_FACE = "#FAFAFA"
LAYER_FACE = "#FCFCFC"

DIRECT = "#333333"
ASSOC = "#888888"


                                                              
         
                                                              

def box(
    ax,
    x,
    y,
    w,
    h,
    text,
    *,
    edge=LINE,
    face=BOX_FACE,
    lw=1.2,
    fontsize=9.0,
    weight="normal",
    linestyle="-",
    radius=0.012,
    text_color=TEXT,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
        linestyle=linestyle,
        zorder=2,
    )
    ax.add_patch(patch)

    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        color=text_color,
        linespacing=1.18,
        zorder=3,
    )

    return patch


def layer(ax, y, h):
    patch = FancyBboxPatch(
        (0.045, y),
        0.91,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.010",
        linewidth=0.8,
        edgecolor="#DDDDDD",
        facecolor=LAYER_FACE,
        zorder=0,
    )
    ax.add_patch(patch)


def arrow(
    ax,
    x1,
    y1,
    x2,
    y2,
    *,
    color=DIRECT,
    lw=1.4,
    linestyle="-",
    mutation_scale=12,
    connectionstyle="arc3,rad=0",
):
    a = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=mutation_scale,
        linewidth=lw,
        color=color,
        linestyle=linestyle,
        connectionstyle=connectionstyle,
        shrinkA=0,
        shrinkB=0,
        zorder=1,
    )
    ax.add_patch(a)
    return a


def connector(
    ax,
    x1,
    y1,
    x2,
    y2,
    *,
    color=ASSOC,
    lw=1.3,
    linestyle="--",
):
    ax.plot(
        [x1, x2],
        [y1, y2],
        color=color,
        linewidth=lw,
        linestyle=linestyle,
        solid_capstyle="round",
        zorder=1,
    )


def badge(
    ax,
    x,
    y,
    text,
    *,
    edge=LIGHT_LINE,
    face="#FFFFFF",
    fontsize=7.4,
    color=MUTED,
):
    ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=color,
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": face,
            "edgecolor": edge,
            "linewidth": 0.8,
        },
        zorder=4,
    )


def panel_letter(ax, letter, y):
    ax.text(
        0.020,
        y,
        letter,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=TEXT,
    )


                                                              
        
                                                              

fig, ax = plt.subplots(figsize=(16.5, 11.2))

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")


                                                              
                                            
                                                              

layer(ax, 0.705, 0.255)
panel_letter(ax, "A", 0.940)


box(
    ax,
    0.085,
    0.840,
    0.135,
    0.060,
    "Топлинен\nстрес",
    edge=HEAT,
    face="#FBF5F5",
    weight="bold",
)

box(
    ax,
    0.085,
    0.750,
    0.135,
    0.060,
    "Студов\nстрес",
    edge=COLD,
    face="#F4F7FB",
    weight="bold",
)


box(
    ax,
    0.335,
    0.790,
    0.235,
    0.095,
    "Различими субстратни\nфосфорилационни траектории",
    edge=STRUCT,
    face="#F7F8F9",
    weight="bold",
)

arrow(
    ax,
    0.220,
    0.870,
    0.335,
    0.845,
    color=HEAT,
)

arrow(
    ax,
    0.220,
    0.780,
    0.335,
    0.830,
    color=COLD,
)


box(
    ax,
    0.690,
    0.790,
    0.220,
    0.095,
    "PKA-свързана\nсубстратна програма",
    edge=PKA,
    face="#F8F4FA",
    lw=1.6,
    weight="bold",
)

arrow(
    ax,
    0.570,
    0.837,
    0.690,
    0.837,
    color=PKA,
    lw=1.7,
)

badge(
    ax,
    0.800,
    0.744,
    "независима репликация",
    edge="#CABBD1",
    face="#FCFAFD",
    color=PKA,
)


                                                              
                                              
 
            
                                                       
               
                                                              

layer(ax, 0.385, 0.285)
panel_letter(ax, "B", 0.650)


                               

box(
    ax,
    0.080,
    0.545,
    0.180,
    0.065,
    "Локален структурен\nконтекст",
    edge=STRUCT,
    face="#F7F8F9",
    weight="bold",
)

box(
    ax,
    0.330,
    0.555,
    0.225,
    0.070,
    "Cdc19-S22\nограничена достъпност",
    edge=STRUCT,
    face="#FFFFFF",
)

box(
    ax,
    0.330,
    0.435,
    0.225,
    0.080,
    "Atg1-S515 · Maf1-S90 · Msn4-S316\nгъвкави регулаторни региони",
    edge=STRUCT,
    face="#FFFFFF",
)

arrow(
    ax,
    0.260,
    0.578,
    0.330,
    0.590,
    color=STRUCT,
    lw=1.2,
)

arrow(
    ax,
    0.260,
    0.565,
    0.330,
    0.475,
    color=STRUCT,
    lw=1.2,
)

badge(
    ax,
    0.443,
    0.530,
    "state-dependent accessibility hypothesis",
    fontsize=7.0,
)

                                 

box(
    ax,
    0.650,
    0.535,
    0.145,
    0.070,
    "Yef3-T972",
    edge=YEF3,
    face="#FAF8F4",
    lw=1.5,
    weight="bold",
)

box(
    ax,
    0.810,
    0.520,
    0.135,
    0.100,
    "SG/RNP-свързано\nфосфопротеомно\nсъстояние",
    edge=YEF3,
    face="#FAF8F4",
)

connector(
    ax,
    0.795,
    0.570,
    0.810,
    0.570,
    color=YEF3,
    lw=1.7,
    linestyle="--",
)

badge(
    ax,
    0.723,
    0.470,
    "PRE-CAUSAL",
    edge="#D5C7B6",
    face="#FCFAF7",
    color=YEF3,
)

ax.text(
    0.875,
    0.455,
    "асоциация след global-state корекция",
    ha="center",
    va="center",
    fontsize=7.0,
    color=MUTED,
)


                                                              
                              
                                                              

layer(ax, 0.075, 0.270)
panel_letter(ax, "C", 0.325)


             

box(
    ax,
    0.080,
    0.225,
    0.170,
    0.065,
    "Guk1-7 структури",
    edge=STRUCT,
    face="#F7F8F9",
    weight="bold",
)

box(
    ax,
    0.355,
    0.225,
    0.195,
    0.065,
    "Късна постстресова\nрезолюция",
    edge=STRUCT,
    face="#FFFFFF",
)

arrow(
    ax,
    0.250,
    0.258,
    0.355,
    0.258,
    color=STRUCT,
    lw=1.3,
)


             

box(
    ax,
    0.080,
    0.115,
    0.170,
    0.065,
    "Pab1 stress granules",
    edge=STRUCT,
    face="#F7F8F9",
    weight="bold",
)

box(
    ax,
    0.355,
    0.115,
    0.195,
    0.065,
    "Постстресова\nрезолюция",
    edge=STRUCT,
    face="#FFFFFF",
)

arrow(
    ax,
    0.250,
    0.148,
    0.355,
    0.148,
    color=STRUCT,
    lw=1.3,
)


        

box(
    ax,
    0.725,
    0.165,
    0.150,
    0.080,
    "Hsp104",
    edge=HSP104,
    face="#F4F8F6",
    lw=1.7,
    fontsize=10,
    weight="bold",
)

                                

arrow(
    ax,
    0.725,
    0.218,
    0.550,
    0.258,
    color=ASSOC,
    lw=1.4,
    linestyle="--",
    connectionstyle="arc3,rad=0.06",
)

ax.text(
    0.630,
    0.272,
    "descriptive support",
    ha="center",
    va="center",
    fontsize=7.0,
    color=MUTED,
)

                              

arrow(
    ax,
    0.725,
    0.190,
    0.550,
    0.148,
    color=HSP104,
    lw=1.8,
    linestyle="-",
    connectionstyle="arc3,rad=-0.06",
)

ax.text(
    0.633,
    0.132,
    "perturbational support",
    ha="center",
    va="center",
    fontsize=7.0,
    color=HSP104,
)


                                                              
                                
                                                              

ax.plot(
    [0.080, 0.115],
    [0.035, 0.035],
    color=DIRECT,
    linewidth=1.6,
)

ax.text(
    0.125,
    0.035,
    "плътна линия = директно / възпроизводимо подкрепена връзка",
    ha="left",
    va="center",
    fontsize=7.1,
    color=MUTED,
)

ax.plot(
    [0.475, 0.510],
    [0.035, 0.035],
    color=ASSOC,
    linewidth=1.5,
    linestyle="--",
)

ax.text(
    0.520,
    0.035,
    "прекъсната линия = асоциация, описателна подкрепа или хипотеза",
    ha="left",
    va="center",
    fontsize=7.1,
    color=MUTED,
)

ax.text(
    0.500,
    0.008,
    "Слоевете обобщават различни нива на доказателства; "
    "между A–C не се предполага линейна причинна каскада.",
    ha="center",
    va="bottom",
    fontsize=7.2,
    color="#555555",
    fontstyle="italic",
)


                                                              
        
                                                              

fig.subplots_adjust(
    left=0.02,
    right=0.98,
    top=0.985,
    bottom=0.035,
)

png = OUTDIR / f"{BASENAME}.png"
pdf = OUTDIR / f"{BASENAME}.pdf"
svg = OUTDIR / f"{BASENAME}.svg"

fig.savefig(
    png,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

fig.savefig(
    pdf,
    bbox_inches="tight",
    facecolor="white",
    metadata={
        "Title": "",
        "Author": "",
        "Subject": "",
        "Keywords": "",
    },
)

fig.savefig(
    svg,
    bbox_inches="tight",
    facecolor="white",
    metadata={
        "Title": "",
        "Description": "",
        "Creator": "",
        "Date": None,
    },
)

plt.close(fig)


                                                              
               
                                                              

legend = """Фигура 4.7. Интегрирана карта на доказателствата от фосфопротеомния, структурния и фенотипния анализ. (A) Топлинният и студовият стрес формират различими субстратни фосфорилационни траектории, сред които PKA-свързаната програма получава независима експериментална и външна подкрепа. (B) Приоритетните фосфосайтове се намират в различен локален структурен контекст. Yef3-T972 независимо проследява stress-granule/RNP-свързано фосфопротеомно състояние, но остава предкаузален кандидат. (C) Фенотипните данни подкрепят участие на Hsp104 в постстресовата резолюция. За Guk1-7 тази връзка е описателна, докато при Pab1 stress granules е подкрепена чрез директна генетична пертурбация. Прекъснатите връзки обозначават асоциативна, описателна или хипотетична подкрепа и не представляват доказана причинна последователност.
"""

legend_path = LEGENDDIR / "Fig4_7_legend_bg.txt"
legend_path.write_text(
    legend,
    encoding="utf-8",
)


                                                              
            
                                                              

script_path = Path(__file__).resolve()

sha256 = hashlib.sha256(
    script_path.read_bytes()
).hexdigest()

manifest = OUTDIR / f"{BASENAME}.manifest.tsv"

manifest.write_text(
    "\n".join([
        "field\tvalue",
        f"figure\t{BASENAME}",
        "figure_type\tevidence-constrained schematic",
        "new_quantitative_analysis\tno",
        "causal_chain_PKA_Yef3_Hsp104_Guk1\tno",
        "Yef3_T972_status\tpre-causal association",
        "Guk1_Hsp104_status\tdescriptive support",
        "Pab1_Hsp104_status\tperturbational support",
        f"script_sha256\t{sha256}",
    ]) + "\n",
    encoding="utf-8",
)


print("WROTE:")
print(png)
print(pdf)
print(svg)
print(legend_path)
print(manifest)
