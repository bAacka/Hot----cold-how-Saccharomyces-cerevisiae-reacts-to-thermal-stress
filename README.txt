STEFAN — FINAL RESULTS CODE AND FIGURES

SCOPE
=====

This package contains the analysis code supporting the results retained
in Sections 4.1–4.7, together with the final figure-generation scripts
and exported figures.

The package is intended primarily for external scientific/code review.
Large source datasets are not bundled.

All Python comments were removed from this review copy. Interpreter
shebang lines were retained. The analytical operations themselves were
not intentionally altered during comment removal.

The original local project root was removed from the review copy.
Where the original scripts used the project root, "/" now serves as a
neutral placeholder. To execute the scripts, ROOT must therefore be
changed to the appropriate local project/data root.


DIRECTORY STRUCTURE
===================

01_analysis/
    Analysis scripts directly supporting the retained Results.

02_upstream_processing/
    Locked image-processing steps required by the final 190321 WT
    microscopy analysis.

03_support_resource_acquisition/
    External mechanistic-resource inventory and download scripts used
    upstream of the structural analysis.

06_figure_scripts/
    Final Figure 4.1–4.7 generation scripts.

07_final_figures/
    Final PDF, PNG and SVG figure exports.

08_figure_legends/
    Final figure legends.

09_small_frozen_inputs/
    Small frozen project-level inputs required for exact provenance,
    including the final claim matrix and Figure 4.6 residual vectors.


RESULT-TO-CODE MAP
==================

4.1 Dynamic phosphoproteomic response
    01_build_kanshin_truth_set.py
    12_audit_heat_cold_response_classes.py
    13_freeze_response_layer.py

4.2 PKA-associated temperature programme
    16_networkin_temporal_differential.py
    29_pka_direct_site_trajectory_corrected.py
    29_export_corrected_direct_null.py

4.3 Independent PKA validation and replication
    28_uliana_context_interaction.py
    30_uliana_networkin_pka_validation.py
    31_phosphoatlas_pka_replication.py
    32_phosphoatlas_temperature_sensitivity.py

4.4 Structural and mechanistic context
    37_build_mechanistic_site_context.py
    38_complete_mechanistic_context.py
    39_structural_tractability_phase1.py
    40_experimental_structure_mapping_audit.py
    41_exact_experimental_residue_coverage.py
    42_cdc19_s22_native_environment.py
    43_disordered_phosphosite_screen.py

4.5 Guk1-7-GFP recovery phenotype
    Final focus-metric freeze/application/validation:
        65_freeze_focus_metric_v3.py
        66_apply_frozen_focus_metric_190321.py
        67_compare_frozen_automated_vs_andersson_S10.py

    Replicate phenotype analysis:
        46h_segment_190326_ssa12dd.py
        46k_track_focus_190326_T0_T8.py
        46l_track_focus_190319_T0_T8.py
        46m_compare_190319_190326_T1_T8.py
        46n_freeze_Guk1_recovery_phenotype.py

    Locked upstream 190321 WT image-processing scripts are retained
    separately under 02_upstream_processing/.

4.6 Yef3-T972 / Hsp104 / recovery evidence
    Yef3-T972:
        49b_test_T972_proteostasis_programs.py
        49c_test_global_axis_specificity.py
        export_exact_49c_SGRNP_residuals.py

    Guk1 / Hsp104:
        46s1_segment_190418_WT_hsp104D.py
        46s2_track_focus_190418_WT_hsp104D_T0_T8.py
        46r2_blinded_spatial_identity_match.py
        46t_finalize_hsp104_validation.py
        47_freeze_hsp104_evidence.py

    Public Pab1 / Hsp104 perturbation:
        52a_download_and_inventory_public_bridge.py
        52b_exact_schema_audit.py
        52c_quantify_public_perturbational_bridge.py

4.7 Integrated evidence map
    No new statistical analysis is introduced in 4.7.
    The figure is an evidence-synthesis schematic generated from the
    results established in Sections 4.1–4.6.


FIGURE 4.6C EXACT RESIDUAL PROVENANCE
=====================================

The Figure 4.6C residual vectors were originally exported from an exact
deterministic replay of the frozen 49c implementation.

The permanent script:
    export_exact_49c_SGRNP_residuals.py

reproduces the four authoritative 49c frozen outputs and exports the
107 Yef3-T972 / stress-granule-RNP residual pairs used in the figure.

The replay was validated against the frozen analysis output before this
package was finalized.


FIGURES
=======

The final Figure 4.1 output is the version with suffix "_v2".
Figures 4.2–4.7 use their final unsuffixed names.

For each figure, PDF, PNG and SVG exports are provided.


REPRODUCIBILITY NOTE
====================

This is a code-and-figure review package, not a complete redistribution
of all raw and public source datasets.

Scripts preserve their original relative project structure. A complete
rerun therefore additionally requires the corresponding source datasets
and external resources at the relative paths referenced by the code.
