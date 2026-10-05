#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
import hashlib
import json

import numpy as np
import pandas as pd


ROOT = Path(
    "/"
)

D = (
    ROOT /
    "44_phenotype_bridge" /
    "52_public_perturbational_bridge"
)

GUK1 = (
    ROOT /
    "44_phenotype_bridge" /
    "46_phenotype_state" /
    "190418_WT_hsp104D_tracked_focus_T0_T8" /
    "46t_final_hsp104_validation" /
    "resolved_genotype_early_late_summary.tsv"
)

HSP104_FREEZE = (
    ROOT /
    "44_phenotype_bridge" /
    "46_phenotype_state" /
    "47_hsp104_evidence_freeze" /
    "47_HSP104_EVIDENCE_FREEZE.txt"
)


                                                              
                                                     
                                                            
                                                     
                                                              

hits = sorted(
    D.glob(
        "52b_Shattuck_Figure6_raw*.tsv"
    )
)

if len(hits) != 1:
    raise RuntimeError(
        f"Expected exactly one Figure 6 raw TSV; got: {hits}"
    )

FIG6 = hits[0]

for p in [
    FIG6,
    GUK1,
    HSP104_FREEZE,
]:
    if not p.is_file():
        raise RuntimeError(
            f"Missing required input: {p}"
        )


                                                              
               
 
                       
 
                      
                       
 
                
                          
                     
                     
                                        
 
             
                        
                          
 
           
             
 
                  
               
                                                              

x = pd.read_csv(
    FIG6,
    sep="\t",
    index_col=0,
)

x.index = pd.to_numeric(
    x.index,
    errors="raise",
).astype(int)


def col(k):
    if str(k) in x.columns:
        return str(k)
    if k in x.columns:
        return k
    raise KeyError(k)


def vals(row, cols):
    z = pd.to_numeric(
        x.loc[
            row,
            [col(c) for c in cols],
        ],
        errors="coerce",
    ).dropna()

    return z.to_numpy(
        dtype=float
    )


LEFT = list(range(2, 7))
RIGHT = list(range(12, 17))


                                                              
                          
 
                                                       
 
                                                     
                                                                
                                                              

def exact_partition_test(a, b):

    a = np.asarray(
        a,
        dtype=float,
    )

    b = np.asarray(
        b,
        dtype=float,
    )

    pool = np.concatenate(
        [a, b]
    )

    na = len(a)

    obs = (
        b.mean()
        -
        a.mean()
    )

    null = []

    n = len(pool)

    for idx in combinations(
        range(n),
        na,
    ):

        take = np.zeros(
            n,
            dtype=bool,
        )

        take[list(idx)] = True

        aa = pool[take]
        bb = pool[~take]

        null.append(
            bb.mean()
            -
            aa.mean()
        )

    null = np.asarray(
        null,
        dtype=float,
    )

    p = np.mean(
        np.abs(null)
        >=
        abs(obs) - 1e-12
    )

    return (
        float(obs),
        float(p),
        int(len(null)),
    )


def cliffs_delta(a, b):
    """
    Positive:
        values in b tend to exceed values in a.
    """

    a = np.asarray(
        a,
        dtype=float,
    )

    b = np.asarray(
        b,
        dtype=float,
    )

    scores = []

    for aa in a:
        for bb in b:
            scores.append(
                np.sign(
                    bb - aa
                )
            )

    return float(
        np.mean(scores)
    )


                                                              
                                    
                   
                                                              

time_rows = {
    "30min_heatshock":
        20,

    "1h_recovery":
        21,

    "2h_recovery":
        22,
}

primary_rows = []

for label, row in time_rows.items():

    hsp = vals(
        row,
        LEFT,
    )

    mut = vals(
        row,
        RIGHT,
    )

    if (
        len(hsp) != 5
        or
        len(mut) != 5
    ):
        raise RuntimeError(
            f"Unexpected replicate count for {label}: "
            f"HSP104={len(hsp)}, hsp104D={len(mut)}"
        )

    delta, p, nperm = (
        exact_partition_test(
            hsp,
            mut,
        )
    )

    primary_rows.append({
        "phenotype":
            "Pab1_foci_percent",

        "condition":
            label,

        "HSP104_n":
            len(hsp),

        "hsp104D_n":
            len(mut),

        "HSP104_mean":
            float(
                hsp.mean()
            ),

        "HSP104_sd":
            float(
                hsp.std(
                    ddof=1
                )
            ),

        "hsp104D_mean":
            float(
                mut.mean()
            ),

        "hsp104D_sd":
            float(
                mut.std(
                    ddof=1
                )
            ),

        "delta_hsp104D_minus_HSP104":
            delta,

        "fold_hsp104D_over_HSP104":
            (
                float(
                    mut.mean()
                    /
                    hsp.mean()
                )
                if hsp.mean() > 0
                else np.nan
            ),

        "cliffs_delta":
            cliffs_delta(
                hsp,
                mut,
            ),

        "exact_randomization_p_two_sided":
            p,

        "n_exact_partitions":
            nperm,
    })


primary = pd.DataFrame(
    primary_rows
)

PRIMARY_OUT = (
    D /
    "52c_PAB1_HSP104_PRIMARY_EFFECTS.tsv"
)

primary.to_csv(
    PRIMARY_OUT,
    sep="\t",
    index=False,
)


                                                              
                               
                                                             
 
                            
 
          
           
 
                 
            
 
             
            
                                                              

sky1 = vals(
    6,
    RIGHT,
)

sky1_kd = vals(
    14,
    RIGHT,
)

gfp = vals(
    22,
    RIGHT,
)

if not (
    len(sky1)
    ==
    len(sky1_kd)
    ==
    len(gfp)
    ==
    5
):
    raise RuntimeError(
        "Unexpected Figure 6B replicate count."
    )


def comparison_row(
    name,
    a_name,
    a,
    b_name,
    b,
):
    delta, p, nperm = (
        exact_partition_test(
            a,
            b,
        )
    )

    return {
        "comparison":
            name,

        "group_A":
            a_name,

        "group_B":
            b_name,

        "group_A_mean":
            float(
                np.mean(a)
            ),

        "group_B_mean":
            float(
                np.mean(b)
            ),

        "B_minus_A":
            delta,

        "exact_randomization_p_two_sided":
            p,

        "cliffs_delta_B_vs_A":
            cliffs_delta(
                a,
                b,
            ),

        "n_exact_partitions":
            nperm,
    }


secondary = pd.DataFrame([
    comparison_row(
        "Sky1_overexpression_vs_GFP_in_hsp104D_2h",
        "Sky1-GFP",
        sky1,
        "GFP",
        gfp,
    ),

    comparison_row(
        "Sky1_K187M_vs_Sky1_in_hsp104D_2h",
        "Sky1-GFP",
        sky1,
        "Sky1K187M-GFP",
        sky1_kd,
    ),
])

SECONDARY_OUT = (
    D /
    "52c_PAB1_SKY1_RESCUE_EFFECTS.tsv"
)

secondary.to_csv(
    SECONDARY_OUT,
    sep="\t",
    index=False,
)


                                                              
                    
                                                              

g = pd.read_csv(
    GUK1,
    sep="\t",
)

g = g.set_index(
    "genotype"
)

if not {
    "WT",
    "hsp104D",
}.issubset(
    g.index
):
    raise RuntimeError(
        "Missing WT/hsp104D Guk1 rows."
    )


wt_early = float(
    g.loc[
        "WT",
        "mean_early_AUC",
    ]
)

mut_early = float(
    g.loc[
        "hsp104D",
        "mean_early_AUC",
    ]
)

wt_late = float(
    g.loc[
        "WT",
        "mean_late_AUC",
    ]
)

mut_late = float(
    g.loc[
        "hsp104D",
        "mean_late_AUC",
    ]
)

wt_t15 = float(
    g.loc[
        "WT",
        "mean_late_T15",
    ]
)

mut_t15 = float(
    g.loc[
        "hsp104D",
        "mean_late_T15",
    ]
)


                                                              
                                          
                                                              

guards = {
    "WT early AUC":
        (
            wt_early,
            -0.465496,
        ),

    "hsp104D early AUC":
        (
            mut_early,
            -0.628083,
        ),

    "WT late AUC":
        (
            wt_late,
            -0.657923,
        ),

    "hsp104D late AUC":
        (
            mut_late,
            -0.140885,
        ),

    "WT T15/T8":
        (
            wt_t15,
            -0.242901,
        ),

    "hsp104D T15/T8":
        (
            mut_t15,
            -0.103475,
        ),
}

for label, (
    observed,
    expected,
) in guards.items():

    if not np.isclose(
        observed,
        expected,
        atol=1e-5,
        rtol=0,
    ):
        raise RuntimeError(
            f"Upstream Guk1 freeze drift: "
            f"{label}: "
            f"{observed} != {expected}"
        )


                                                              
                                
                                                              

pub = (
    primary
    .set_index(
        "condition"
    )
)

pab1_hs_hsp = float(
    pub.loc[
        "30min_heatshock",
        "HSP104_mean",
    ]
)

pab1_hs_mut = float(
    pub.loc[
        "30min_heatshock",
        "hsp104D_mean",
    ]
)

pab1_1h_hsp = float(
    pub.loc[
        "1h_recovery",
        "HSP104_mean",
    ]
)

pab1_1h_mut = float(
    pub.loc[
        "1h_recovery",
        "hsp104D_mean",
    ]
)

pab1_2h_hsp = float(
    pub.loc[
        "2h_recovery",
        "HSP104_mean",
    ]
)

pab1_2h_mut = float(
    pub.loc[
        "2h_recovery",
        "hsp104D_mean",
    ]
)

pab1_1h_p = float(
    pub.loc[
        "1h_recovery",
        "exact_randomization_p_two_sided",
    ]
)

pab1_2h_p = float(
    pub.loc[
        "2h_recovery",
        "exact_randomization_p_two_sided",
    ]
)


                                                              
                               
                                                              

ledger = pd.DataFrame([

    {
        "evidence_id":
            "PUBLIC_PAB1_HSP104_FORMATION",

        "system":
            "Shattuck2019_Pab1",

        "perturbation":
            "HSP104_loss",

        "phase":
            "30min_heatshock",

        "result":
            "FORMATION_PRESERVED",

        "effect":
            (
                pab1_hs_mut
                -
                pab1_hs_hsp
            ),

        "p_value":
            float(
                pub.loc[
                    "30min_heatshock",
                    "exact_randomization_p_two_sided",
                ]
            ),

        "interpretation":
            (
                "Hsp104 loss does not materially prevent "
                "Pab1 focus formation during acute heat shock."
            ),
    },

    {
        "evidence_id":
            "PUBLIC_PAB1_HSP104_1H_RECOVERY",

        "system":
            "Shattuck2019_Pab1",

        "perturbation":
            "HSP104_loss",

        "phase":
            "1h_recovery",

        "result":
            "STRONG_RECOVERY_DEFECT",

        "effect":
            (
                pab1_1h_mut
                -
                pab1_1h_hsp
            ),

        "p_value":
            pab1_1h_p,

        "interpretation":
            (
                "Pab1 foci persist strongly in hsp104D "
                "during recovery."
            ),
    },

    {
        "evidence_id":
            "PUBLIC_PAB1_HSP104_2H_RECOVERY",

        "system":
            "Shattuck2019_Pab1",

        "perturbation":
            "HSP104_loss",

        "phase":
            "2h_recovery",

        "result":
            "VERY_STRONG_RECOVERY_DEFECT",

        "effect":
            (
                pab1_2h_mut
                -
                pab1_2h_hsp
            ),

        "p_value":
            pab1_2h_p,

        "interpretation":
            (
                "HSP104 cells have almost fully dissolved "
                "Pab1 foci while hsp104D cells retain them "
                "in the large majority of cells."
            ),
    },

    {
        "evidence_id":
            "PUBLIC_PAB1_SKY1_RESCUE",

        "system":
            "Shattuck2019_Pab1",

        "perturbation":
            "Sky1_overexpression_in_hsp104D",

        "phase":
            "2h_recovery",

        "result":
            "PARTIAL_RESCUE",

        "effect":
            float(
                sky1.mean()
                -
                gfp.mean()
            ),

        "p_value":
            float(
                secondary.loc[
                    secondary[
                        "comparison"
                    ].eq(
                        "Sky1_overexpression_vs_GFP_in_hsp104D_2h"
                    ),
                    "exact_randomization_p_two_sided",
                ].iloc[0]
            ),

        "interpretation":
            (
                "Elevated active Sky1 partially compensates "
                "for Hsp104 deficiency in Pab1 granule "
                "resolution."
            ),
    },

    {
        "evidence_id":
            "INTERNAL_GUK1_HSP104_EARLY",

        "system":
            "Andersson_raw_Guk1_reanalysis",

        "perturbation":
            "HSP104_loss",

        "phase":
            "T1_T8",

        "result":
            "NO_SUPPORT",

        "effect":
            (
                wt_early
                -
                mut_early
            ),

        "p_value":
            np.nan,

        "interpretation":
            (
                "Frozen early Guk1 AUC does not support "
                "an Hsp104-dependent early recovery defect."
            ),
    },

    {
        "evidence_id":
            "INTERNAL_GUK1_HSP104_LATE",

        "system":
            "Andersson_raw_Guk1_reanalysis",

        "perturbation":
            "HSP104_loss",

        "phase":
            "T8_T15",

        "result":
            "SUPPORTIVE_LATE_DEFECT",

        "effect":
            (
                wt_late
                -
                mut_late
            ),

        "p_value":
            np.nan,

        "interpretation":
            (
                "WT shows substantially greater late-phase "
                "loss of Guk1 focus enrichment than hsp104D."
            ),
    },

    {
        "evidence_id":
            "PUBLISHED_GUK1_HSP104",

        "system":
            "Andersson2021_Guk1",

        "perturbation":
            "HSP104_loss",

        "phase":
            "recovery",

        "result":
            "DIRECT_PUBLISHED_SUPPORT",

        "effect":
            np.nan,

        "p_value":
            np.nan,

        "interpretation":
            (
                "Published study reports that HSP104 deletion "
                "does not abolish Guk1 ring-like aggregate "
                "formation but makes those structures more "
                "stable during recovery."
            ),
    },
])

LEDGER_OUT = (
    D /
    "52c_PUBLIC_PERTURBATIONAL_BRIDGE_EVIDENCE.tsv"
)

ledger.to_csv(
    LEDGER_OUT,
    sep="\t",
    index=False,
)


                                                              
                      
                                                              

result = f"""
PUBLIC PERTURBATIONAL BRIDGE — FINAL RESULT
===========================================

STATUS
------
SUPPORTED.

PRE-SPECIFIED QUESTION
----------------------
Do independently published stress-RNP and Guk1 aggregate recovery
phenotypes respond coherently to perturbation of shared
protein-quality-control machinery?

Primary shared perturbation:

    HSP104 loss / impairment

PUBLIC STRESS-RNP RESULT
------------------------
Source:
Shattuck et al. 2019, Nature Communications
DOI: 10.1038/s41467-019-11550-w

Phenotype:
percentage of cells retaining Pab1 stress-granule foci.

GFP-only control background:

After 30 min heat shock:

    HSP104:
        {pab1_hs_hsp:.6f} %

    hsp104D:
        {pab1_hs_mut:.6f} %

Interpretation:
Pab1 focus formation is essentially intact in hsp104D.

After 1 h recovery:

    HSP104:
        {pab1_1h_hsp:.6f} %

    hsp104D:
        {pab1_1h_mut:.6f} %

    hsp104D - HSP104:
        {pab1_1h_mut - pab1_1h_hsp:.6f} percentage points

    exact two-sided randomization P:
        {pab1_1h_p:.8f}

After 2 h recovery:

    HSP104:
        {pab1_2h_hsp:.6f} %

    hsp104D:
        {pab1_2h_mut:.6f} %

    hsp104D - HSP104:
        {pab1_2h_mut - pab1_2h_hsp:.6f} percentage points

    exact two-sided randomization P:
        {pab1_2h_p:.8f}

Interpretation:
Hsp104 loss produces a profound defect in Pab1 stress-granule
dissolution during recovery despite essentially normal focus
formation during heat shock.

SECONDARY PUBLIC MECHANISTIC RESULT
-----------------------------------
At 2 h recovery in hsp104D:

    Sky1-GFP:
        {sky1.mean():.6f} % cells with Pab1 foci

    GFP control:
        {gfp.mean():.6f} % cells with Pab1 foci

    Sky1(K187M)-GFP:
        {sky1_kd.mean():.6f} % cells with Pab1 foci

Active Sky1 partially compensates for Hsp104 deficiency, whereas
kinase-dead Sky1 does not reproduce that rescue.

This supports overlapping / compensatory recovery machinery rather
than a single linear pathway.

GUK1 RESULT
-----------
Source:
Andersson et al. 2021
DOI: 10.1371/journal.pgen.1008951

Published direct observation:
HSP104 deletion does not abolish formation of the nucleolar
ring-like Guk1-7 aggregate structures but makes them more stable
during recovery.

Independent raw-image reanalysis:

Early frozen T1:T8 AUC:

    WT:
        {wt_early:.6f}

    hsp104D:
        {mut_early:.6f}

    WT - hsp104D:
        {wt_early - mut_early:.6f}

The early endpoint does NOT support an Hsp104-dependent recovery
defect.

Late T8:T15 AUC:

    WT:
        {wt_late:.6f}

    hsp104D:
        {mut_late:.6f}

    WT - hsp104D:
        {wt_late - mut_late:.6f}

Late endpoint log2(T15/T8):

    WT:
        {wt_t15:.6f}

    hsp104D:
        {mut_t15:.6f}

The internal raw-image result therefore supports a later
Hsp104-dependent resolution component.

CROSS-DATASET CONCLUSION
------------------------
The pre-specified public perturbational bridge is SUPPORTED.

Two distinct heat-induced assemblies:

    Pab1 stress granules / stress-RNP state

and

    Guk1-7 misfolded-protein aggregates

both show preserved or substantial initial assembly but impaired
subsequent resolution when Hsp104 function is absent.

Therefore the data support:

    a shared Hsp104-dependent heat-recovery / disaggregation layer.

This is stronger than a generic literature analogy because the same
molecular perturbation is connected to recovery phenotypes in two
independent experimental systems.

CLAIM BOUNDARY
--------------
This analysis does NOT establish:

    stress granule -> Guk1 aggregate

It does NOT establish:

    YEF3-T972 -> Guk1 recovery

It does NOT establish:

    YEF3-T972 -> Hsp104

It does NOT establish that Pab1 and Guk1 are processed by an
identical molecular mechanism.

The appropriate interpretation is convergence on shared
protein-quality-control / recovery machinery.

PROJECT ROLE
------------
This result connects the independently observed YEF3-T972 /
stress-RNP phosphoproteomic state to the Guk1 phenotype at the level
of a biologically validated shared recovery process.

YEF3-T972 itself remains pre-causal.

Direct site-specific causality would still require perturbing T972
and measuring the frozen Guk1 phenotype.
""".strip() + "\n"

RESULT_OUT = (
    D /
    "52c_PUBLIC_BRIDGE_RESULT.txt"
)

RESULT_OUT.write_text(
    result
)


                                                              
                         
                                                              

machine = {

    "status":
        "SUPPORTED",

    "primary_shared_perturbation":
        "HSP104_loss",

    "pab1": {
        "heatshock_HSP104_percent":
            pab1_hs_hsp,

        "heatshock_hsp104D_percent":
            pab1_hs_mut,

        "recovery_1h_HSP104_percent":
            pab1_1h_hsp,

        "recovery_1h_hsp104D_percent":
            pab1_1h_mut,

        "recovery_1h_exact_p":
            pab1_1h_p,

        "recovery_2h_HSP104_percent":
            pab1_2h_hsp,

        "recovery_2h_hsp104D_percent":
            pab1_2h_mut,

        "recovery_2h_exact_p":
            pab1_2h_p,
    },

    "guk1": {
        "early_WT_AUC":
            wt_early,

        "early_hsp104D_AUC":
            mut_early,

        "late_WT_AUC":
            wt_late,

        "late_hsp104D_AUC":
            mut_late,

        "late_WT_log2_T15_vs_T8":
            wt_t15,

        "late_hsp104D_log2_T15_vs_T8":
            mut_t15,
    },

    "supported_claim":
        (
            "Pab1 stress-RNP recovery and Guk1-7 aggregate "
            "recovery converge on a shared Hsp104-dependent "
            "heat-recovery/disaggregation layer."
        ),

    "t972_causal":
        False,

    "t972_status":
        "pre-causal",

    "prohibited_claims": [
        "stress granule causes Guk1 aggregation",
        "YEF3-T972 causes Guk1 recovery",
        "YEF3-T972 acts through Hsp104",
        "Pab1 and Guk1 use identical molecular mechanisms",
    ],
}

JSON_OUT = (
    D /
    "52c_PUBLIC_BRIDGE_RESULT.json"
)

JSON_OUT.write_text(
    json.dumps(
        machine,
        indent=2,
    )
    + "\n"
)


                                                              
                     
                                                              

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

            h.update(b)

    return h.hexdigest()


files = [
    FIG6,
    GUK1,
    HSP104_FREEZE,
    PRIMARY_OUT,
    SECONDARY_OUT,
    LEDGER_OUT,
    RESULT_OUT,
    JSON_OUT,
]

manifest = []

for p in files:

    manifest.append({
        "path":
            str(p),

        "sha256":
            sha256(p),

        "bytes":
            p.stat().st_size,
    })

MANIFEST_OUT = (
    D /
    "52c_PUBLIC_BRIDGE_MANIFEST.tsv"
)

pd.DataFrame(
    manifest
).to_csv(
    MANIFEST_OUT,
    sep="\t",
    index=False,
)


                                                              
                     
                                                              

COMPLETE = (
    D /
    "52c_COMPLETE.txt"
)

COMPLETE.write_text(
    "Public HSP104 perturbational bridge frozen COMPLETE.\n"
)


                                                              
                      
                                                              

print()
print("=" * 80)
print("PAB1 PRIMARY PUBLIC HSP104 EFFECT")
print("=" * 80)

print(
    primary.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6f}",
    )
)

print()
print("=" * 80)
print("PAB1 SKY1 RESCUE")
print("=" * 80)

print(
    secondary.to_string(
        index=False,
        float_format=lambda z:
            f"{z:.6f}",
    )
)

print()
print("=" * 80)
print("FINAL CROSS-DATASET RESULT")
print("=" * 80)

print(result)

print(
    "52_PUBLIC_BRIDGE_STATUS=SUPPORTED"
)

