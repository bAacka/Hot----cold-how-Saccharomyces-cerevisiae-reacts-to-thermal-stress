#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd


ROOT = Path(
    "/"
)

BRIDGE = ROOT / "44_phenotype_bridge"
STATE = BRIDGE / "46_phenotype_state"

FINAL = (
    STATE /
    "190418_WT_hsp104D_tracked_focus_T0_T8" /
    "46t_final_hsp104_validation"
)

MAP = (
    STATE /
    "46r_hsp104_genotype_resolution" /
    "46r4_HSP104_SCENE_MAPPING_FREEZE.tsv"
)

OUT = (
    STATE /
    "47_hsp104_evidence_freeze"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


                                                              
                 
                                                              

SUMMARY = (
    FINAL /
    "resolved_genotype_early_late_summary.tsv"
)

CELLS = (
    FINAL /
    "fixed_cell_early_and_late_endpoints.tsv"
)

TRAJ = (
    FINAL /
    "resolved_genotype_T0_T15_trajectory.tsv"
)

FINAL_TXT = (
    FINAL /
    "FINAL_HSP104_VALIDATION.txt"
)

required = [
    SUMMARY,
    CELLS,
    TRAJ,
    FINAL_TXT,
    MAP,
]

for p in required:
    if not p.is_file():
        raise RuntimeError(
            f"Missing required input: {p}"
        )


                                                              
                          
                                                              

s = pd.read_csv(
    SUMMARY,
    sep="\t",
)

if set(s["genotype"]) != {
    "WT",
    "hsp104D",
}:
    raise RuntimeError(
        "Expected WT and hsp104D summary rows."
    )

q = s.set_index(
    "genotype"
)

wt_early = float(
    q.loc[
        "WT",
        "mean_early_AUC",
    ]
)

mut_early = float(
    q.loc[
        "hsp104D",
        "mean_early_AUC",
    ]
)

wt_late = float(
    q.loc[
        "WT",
        "mean_late_AUC",
    ]
)

mut_late = float(
    q.loc[
        "hsp104D",
        "mean_late_AUC",
    ]
)

wt_end = float(
    q.loc[
        "WT",
        "mean_late_T15",
    ]
)

mut_end = float(
    q.loc[
        "hsp104D",
        "mean_late_T15",
    ]
)

early_delta = (
    wt_early -
    mut_early
)

late_delta = (
    wt_late -
    mut_late
)


                                                              
                                                       
                                      
                                                              

expected = {
    "wt_early":
        -0.465496,

    "mut_early":
        -0.628083,

    "wt_late":
        -0.657923,

    "mut_late":
        -0.140885,

    "wt_end":
        -0.242901,

    "mut_end":
        -0.103475,
}

observed = {
    "wt_early":
        wt_early,

    "mut_early":
        mut_early,

    "wt_late":
        wt_late,

    "mut_late":
        mut_late,

    "wt_end":
        wt_end,

    "mut_end":
        mut_end,
}

for k in expected:

    if not np.isclose(
        observed[k],
        expected[k],
        atol=1e-5,
        rtol=0,
    ):
        raise RuntimeError(
            f"Frozen-result mismatch for {k}: "
            f"{observed[k]} vs {expected[k]}"
        )


                                                              
                      
                                                              

m = pd.read_csv(
    MAP,
    sep="\t",
)

wanted = (
    m[
        m["include_primary"]
        .astype(str)
        .str.upper()
        .eq("YES")
    ]
    [
        [
            "scene",
            "genotype",
            "evidence",
            "confidence",
        ]
    ]
    .sort_values("scene")
)

mapping = {
    int(r.scene):
        str(r.genotype)
    for r in wanted.itertuples()
}

if mapping != {
    0: "hsp104D",
    3: "WT",
}:
    raise RuntimeError(
        f"Unexpected frozen scene map: {mapping}"
    )


                                                              
                 
                                                              

rows = [

    {
        "evidence_id":
            "HSP104_INTERNAL_PRIMARY_EARLY",

        "evidence_class":
            "internal_image_analysis",

        "analysis_status":
            "PRIMARY_FROZEN",

        "result":
            "NO_SUPPORT_FOR_PREDICTED_HSP104_EFFECT",

        "detail":
            (
                "T1:T8 integrated focus/cytosol reduction "
                "was slightly greater in hsp104D than WT."
            ),

        "WT_value":
            wt_early,

        "hsp104D_value":
            mut_early,

        "WT_minus_hsp104D":
            early_delta,

        "interpretation":
            (
                "No evidence that Hsp104 is required for "
                "the initial T1:T8 reduction in Guk1-7 "
                "focus enrichment."
            ),

        "inferential_scope":
            "descriptive_one_scene_per_genotype",

        "source":
            str(SUMMARY),
    },

    {
        "evidence_id":
            "HSP104_INTERNAL_SECONDARY_LATE",

        "evidence_class":
            "internal_image_analysis",

        "analysis_status":
            "SECONDARY_POST_PRIMARY_FIXED_CELLS",

        "result":
            "SUPPORTIVE",

        "detail":
            (
                "Same preselected cells showed substantially "
                "greater T8:T15 reduction in WT."
            ),

        "WT_value":
            wt_late,

        "hsp104D_value":
            mut_late,

        "WT_minus_hsp104D":
            late_delta,

        "interpretation":
            (
                "WT undergoes greater late-phase loss of "
                "Guk1-7 focus enrichment than hsp104D."
            ),

        "inferential_scope":
            "descriptive_one_scene_per_genotype",

        "source":
            str(SUMMARY),
    },

    {
        "evidence_id":
            "HSP104_INTERNAL_LATE_ENDPOINT",

        "evidence_class":
            "internal_image_analysis",

        "analysis_status":
            "SECONDARY_ENDPOINT",

        "result":
            "SUPPORTIVE",

        "detail":
            (
                "T15 relative to T8 remained more reduced "
                "in WT than hsp104D."
            ),

        "WT_value":
            wt_end,

        "hsp104D_value":
            mut_end,

        "WT_minus_hsp104D":
            wt_end - mut_end,

        "interpretation":
            (
                "Late endpoint is consistent with more "
                "complete WT aggregate/focus resolution."
            ),

        "inferential_scope":
            "descriptive_one_scene_per_genotype",

        "source":
            str(SUMMARY),
    },

    {
        "evidence_id":
            "HSP104_PUBLISHED_DIRECT_GUK1",

        "evidence_class":
            "published_same_system",

        "analysis_status":
            "EXTERNAL",

        "result":
            "STRONG_SUPPORT",

        "detail":
            (
                "Andersson et al. 2021 report that deleting "
                "HSP104 did not abolish nucleolar ring-like "
                "Guk1-7 aggregate formation but made the "
                "structures more stable during recovery."
            ),

        "WT_value":
            np.nan,

        "hsp104D_value":
            np.nan,

        "WT_minus_hsp104D":
            np.nan,

        "interpretation":
            (
                "Direct published support for an Hsp104 role "
                "in recovery/resolution rather than initial "
                "formation."
            ),

        "inferential_scope":
            "published_direct",

        "source":
            "doi:10.1371/journal.pgen.1008951",
    },

    {
        "evidence_id":
            "HSP104_PUBLISHED_ASSAY_MATCH",

        "evidence_class":
            "published_methodological_support",

        "analysis_status":
            "EXTERNAL",

        "result":
            "SUPPORTIVE",

        "detail":
            (
                "Andersson et al. quantified Guk1-7 recovery "
                "using longitudinal focus:cytosol GFP intensity "
                "ratios in tracked cells."
            ),

        "WT_value":
            np.nan,

        "hsp104D_value":
            np.nan,

        "WT_minus_hsp104D":
            np.nan,

        "interpretation":
            (
                "Internal metric is biologically aligned with "
                "the published Guk1-7 recovery assay."
            ),

        "inferential_scope":
            "methodological",

        "source":
            "doi:10.1371/journal.pgen.1008951",
    },

    {
        "evidence_id":
            "HSP104_PUBLISHED_GENERAL_DISAGGREGATION",

        "evidence_class":
            "published_independent_system",

        "analysis_status":
            "EXTERNAL",

        "result":
            "SUPPORTIVE",

        "detail":
            (
                "Heat-induced Pab1-GFP stress granule "
                "resolution is dramatically slower in "
                "hsp104D cells."
            ),

        "WT_value":
            np.nan,

        "hsp104D_value":
            np.nan,

        "WT_minus_hsp104D":
            np.nan,

        "interpretation":
            (
                "Independent support for Hsp104 acting during "
                "post-stress aggregate/granule resolution."
            ),

        "inferential_scope":
            "mechanistic_context",

        "source":
            "doi:10.1371/journal.pgen.1011424",
    },

    {
        "evidence_id":
            "HSP104_190417_REPLICATION",

        "evidence_class":
            "internal_candidate_replication",

        "analysis_status":
            "NOT_ANALYZED",

        "result":
            "GENOTYPE_UNRESOLVED",

        "detail":
            (
                "190417 -05 contains four scenes but no "
                "defensible scene-to-genotype annotation. "
                "Published-movie spatial matching was not "
                "sufficiently discriminative."
            ),

        "WT_value":
            np.nan,

        "hsp104D_value":
            np.nan,

        "WT_minus_hsp104D":
            np.nan,

        "interpretation":
            (
                "Not used as genotype replication; assigning "
                "genotype from aggregate persistence would "
                "be circular."
            ),

        "inferential_scope":
            "excluded",

        "source":
            (
                "46r3_SYMMETRIC_REFINED_SCENE_MATCHES.tsv; "
                "190417 CZI metadata"
            ),
    },
]

ledger = pd.DataFrame(
    rows
)

LEDGER = (
    OUT /
    "47_HSP104_EVIDENCE_LEDGER.tsv"
)

ledger.to_csv(
    LEDGER,
    sep="\t",
    index=False,
)


                                                              
              
                                                              

claim = f"""
HSP104 EVIDENCE FREEZE — FINAL
==============================

STATUS
------
COMPLETE.

RESOLVED ACQUISITION
--------------------
190418 WT/hsp104D Guk1-7-GFP recovery timelapse -08.czi

Scene 0 = hsp104D
Scene 3 = WT
Scenes 1/2 = unresolved and excluded.

Scene identity was frozen from spatial identity to genotype-labelled
published movies and did not use aggregate persistence or recovery
kinetics.

PRIMARY FROZEN RESULT
---------------------
T1:T8 AUC of log2(focus/cytosol relative to T1)

WT:
    {wt_early:.6f}

hsp104D:
    {mut_early:.6f}

WT - hsp104D:
    {early_delta:.6f}

Interpretation:
The pre-existing early phenotype does NOT support the prediction that
Hsp104 deletion impairs the initial T1:T8 reduction in Guk1-7 focus
enrichment.

This negative/discordant result is retained without reinterpretation.

SECONDARY LATE-PHASE RESULT
---------------------------
Same cells selected before inspection of the late phenotype.
No reselection.
Same segmentation.
Same sequential tracking constraints.
Same radius-4 focus/cytosol measurement.

T8:T15 AUC relative to T8:

WT:
    {wt_late:.6f}

hsp104D:
    {mut_late:.6f}

WT - hsp104D:
    {late_delta:.6f}

T15 relative to T8:

WT:
    {wt_end:.6f}

hsp104D:
    {mut_end:.6f}

Interpretation:
WT cells exhibit substantially greater late-phase loss of Guk1-7
focus enrichment than hsp104D cells.

PUBLISHED DIRECT SUPPORT
------------------------
Andersson et al., PLOS Genetics 2021
DOI: 10.1371/journal.pgen.1008951

The study reports that deleting HSP104 did not abolish formation of
the nucleolar ring-like Guk1-7 aggregate structures, but made them
more stable during recovery.

The study's Guk1-7 recovery analysis also uses focus:cytosol GFP
intensity ratios in longitudinally tracked cells.

Published S4 Movie:
hsp104D guk1-7-GFP recovery after 42 C heat shock,
12 images acquired every 10 minutes from 10 to 120 minutes.

INDEPENDENT MECHANISTIC CONTEXT
-------------------------------
PLOS Genetics 2024
DOI: 10.1371/journal.pgen.1011424

Heat-induced Pab1-GFP stress-granule resolution is substantially
slower in hsp104D, independently supporting a post-stress
disaggregation/resolution function for Hsp104.

190417 ACQUISITION
------------------
NOT used as an independent genotype replicate.

Reason:
scene-to-genotype identity cannot be resolved independently from the
phenotype with sufficient confidence.

Inferring genotype from aggregate persistence would be circular.

FINAL SUPPORTED CLAIM
---------------------
Hsp104 contributes to the late resolution/stability of the Guk1-7
aggregate state during recovery from heat shock.

The data do not support a claim that Hsp104 is required for initial
aggregate formation or for the early T1:T8 decline in focus
enrichment.

INFERENCE LIMIT
---------------
The internal genotype comparison contains one resolved scene per
genotype.

The five tracked cells within each scene are cellular subsamples,
not biological replicates.

Therefore:
    - no genotype P-value is claimed;
    - the internal analysis is descriptive mechanistic support;
    - the direct published Hsp104/Guk1-7 experiment supplies the
      replicated external biological evidence.

PROHIBITED OVERCLAIMS
---------------------
Do NOT claim:

1. Hsp104 is required for Guk1-7 aggregate formation.

2. Hsp104 accelerates the entire recovery trajectory.

3. The T1:T8 primary phenotype supported the Hsp104 hypothesis.

4. The five cells per scene are biological replicates.

5. The internal WT/hsp104D comparison provides a valid genotype
   significance test.

6. 190417 independently replicated the genotype effect.

7. Hsp104 physically occupies the nucleolar ring-like Guk1-7
   structures; the published paper explicitly distinguishes
   aggregate stability from Hsp104 localization there.

FINAL PROJECT USE
-----------------
Evidence role:

    DIRECT EXTERNAL BIOLOGY
        +
    INDEPENDENT RAW-IMAGE QUANTITATIVE SUPPORT
        +
    MECHANISTIC CONSISTENCY

Use Hsp104 as evidence specifically for late aggregate resolution,
not aggregate initiation.
""".strip() + "\n"

FREEZE = (
    OUT /
    "47_HSP104_EVIDENCE_FREEZE.txt"
)

FREEZE.write_text(
    claim
)


                                                              
                        
                                                              

machine = {
    "status":
        "COMPLETE",

    "supported_claim":
        (
            "Hsp104 contributes to late resolution/stability "
            "of the Guk1-7 aggregate state during recovery "
            "from heat shock."
        ),

    "primary_early": {
        "WT_mean_AUC":
            wt_early,

        "hsp104D_mean_AUC":
            mut_early,

        "WT_minus_hsp104D":
            early_delta,

        "supports_predicted_effect":
            False,
    },

    "secondary_late": {
        "WT_mean_AUC":
            wt_late,

        "hsp104D_mean_AUC":
            mut_late,

        "WT_minus_hsp104D":
            late_delta,

        "WT_mean_log2_T15_vs_T8":
            wt_end,

        "hsp104D_mean_log2_T15_vs_T8":
            mut_end,

        "supports_late_resolution_effect":
            True,
    },

    "internal_replication":
        False,

    "internal_replication_reason":
        (
            "190417 scene genotypes unresolved "
            "independently of phenotype."
        ),

    "external_direct_reference":
        "10.1371/journal.pgen.1008951",

    "external_mechanistic_reference":
        "10.1371/journal.pgen.1011424",

    "genotype_inference_valid":
        False,

    "genotype_p_value_claimed":
        False,
}

MACHINE = (
    OUT /
    "47_HSP104_EVIDENCE_FREEZE.json"
)

MACHINE.write_text(
    json.dumps(
        machine,
        indent=2,
    )
    + "\n"
)


                                                              
                 
                                                              

def sha256(p):

    h = hashlib.sha256()

    with open(
        p,
        "rb",
    ) as f:

        while True:

            b = f.read(
                1024 * 1024
            )

            if not b:
                break

            h.update(b)

    return h.hexdigest()


manifest_files = (
    required
    + [
        LEDGER,
        FREEZE,
        MACHINE,
    ]
)

manifest_rows = []

for p in manifest_files:

    manifest_rows.append({
        "path":
            str(p),

        "sha256":
            sha256(p),

        "bytes":
            p.stat().st_size,
    })

manifest = pd.DataFrame(
    manifest_rows
)

MANIFEST = (
    OUT /
    "47_HSP104_EVIDENCE_MANIFEST.tsv"
)

manifest.to_csv(
    MANIFEST,
    sep="\t",
    index=False,
)


                                                              
                     
                                                              

complete = (
    OUT /
    "47_HSP104_COMPLETE.txt"
)

complete.write_text(
    "HSP104 evidence branch frozen COMPLETE.\n"
)


                                                              
       
                                                              

print()
print("=" * 80)
print("HSP104 FINAL EVIDENCE LEDGER")
print("=" * 80)

print(
    ledger[
        [
            "evidence_id",
            "analysis_status",
            "result",
            "WT_minus_hsp104D",
        ]
    ].to_string(
        index=False,
    )
)

print()
print("=" * 80)
print("FINAL CLAIM")
print("=" * 80)

print(
    machine[
        "supported_claim"
    ]
)

print()
print("=" * 80)
print("OUTPUT FILES")
print("=" * 80)

for p in sorted(
    OUT.iterdir()
):
    print(p.name)

print()
print("HSP104_BRANCH_STATUS=COMPLETE")
