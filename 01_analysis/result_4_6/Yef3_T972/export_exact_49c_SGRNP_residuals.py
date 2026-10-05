#!/usr/bin/env python3

from pathlib import Path
import argparse
import hashlib
import re
import shutil
import subprocess
import sys
import tempfile

import pandas as pd


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--project-root",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--output",
        type=Path,
        required=True,
    )

    args = ap.parse_args()

    root = args.project_root.resolve()
    output = args.output.resolve()

    source_dir = (
        root
        / "44_phenotype_bridge"
        / "49_T972_precausal_support"
    )

    original_script = (
        source_dir
        / "49c_test_global_axis_specificity.py"
    )

    if not original_script.is_file():
        raise FileNotFoundError(original_script)

    required_frozen = [
        "49c_GLOBAL_AXIS_DIAGNOSTICS.tsv",
        "49c_PRIMARY_GLOBAL_ADJUSTED_RESULTS.tsv",
        "49c_GLOBAL_ADJUSTED_RANDOM_GENESET_NULL.tsv",
        "49c_DISCOVERY_PERTURBATIONS_GLOBAL_CONTEXT.tsv",
    ]

    for name in required_frozen:
        path = source_dir / name

        if not path.is_file():
            raise FileNotFoundError(path)

    with tempfile.TemporaryDirectory(
        prefix="49c_exact_replay_"
    ) as tmp:
        tmp = Path(tmp)

        replay_dir = (
            tmp
            / "49_T972_precausal_support"
        )

                                                              
                                                          
                                             
        shutil.copytree(
            source_dir,
            replay_dir,
        )

        replay_script = (
            replay_dir
            / "49c_test_global_axis_specificity.py"
        )

        text = replay_script.read_text(
            encoding="utf-8"
        )

                            
                                                              
                                                                  
        out_pattern = re.compile(
            r'''
            OUT
            \s*=\s*
            \(
            \s*ROOT
            \s*/\s*"44_phenotype_bridge"
            \s*/\s*"49_T972_precausal_support"
            \s*
            \)
            ''',
            re.VERBOSE | re.DOTALL,
        )

        text, n = out_pattern.subn(
            f'OUT = Path(r"{replay_dir}")',
            text,
            count=1,
        )

        if n != 1:
            raise RuntimeError(
                "Could not uniquely redirect OUT "
                "in frozen 49c script."
            )

        anchor = (
            "    partial_rho, perm_p = (\n"
        )

        export_block = '''
    if module == "STRESS_GRANULE_RNP":
        residual_export = pd.DataFrame({
            "YEF3_T972_residual":
                target_resid,
            "STRESS_GRANULE_RNP_residual":
                module_resid,
        })

        residual_export.to_csv(
            OUT
            / "EXACT_STRESS_GRANULE_RNP_RESIDUAL_PAIRS.tsv",
            sep="\\t",
            index=False,
        )

'''

        if text.count(anchor) != 1:
            raise RuntimeError(
                "Expected exactly one primary "
                "partial-correlation anchor in 49c; "
                f"found {text.count(anchor)}."
            )

        text = text.replace(
            anchor,
            export_block + anchor,
            1,
        )

        replay_script.write_text(
            text,
            encoding="utf-8",
        )

        result = subprocess.run(
            [
                sys.executable,
                str(replay_script),
            ],
            cwd=str(root),
            text=True,
            capture_output=True,
        )

        if result.returncode != 0:
            print(
                "===== REPLAY STDOUT ====="
            )
            print(result.stdout)

            print(
                "===== REPLAY STDERR =====",
                file=sys.stderr,
            )
            print(
                result.stderr,
                file=sys.stderr,
            )

            raise RuntimeError(
                "Frozen 49c replay failed with "
                f"exit code {result.returncode}"
            )

                                                              
                                                              
                         
                                                              

        for name in required_frozen:
            original = pd.read_csv(
                source_dir / name,
                sep="\t",
            )

            replayed = pd.read_csv(
                replay_dir / name,
                sep="\t",
            )

            pd.testing.assert_frame_equal(
                original,
                replayed,
                check_dtype=False,
                check_exact=False,
                rtol=0,
                atol=1e-12,
            )

            print(
                "REPRODUCED:",
                name,
            )

        residual_file = (
            replay_dir
            / "EXACT_STRESS_GRANULE_RNP_RESIDUAL_PAIRS.tsv"
        )

        if not residual_file.is_file():
            raise RuntimeError(
                "Residual export was not created."
            )

        residual = pd.read_csv(
            residual_file,
            sep="\t",
        )

        expected_columns = [
            "YEF3_T972_residual",
            "STRESS_GRANULE_RNP_residual",
        ]

        if list(residual.columns) != expected_columns:
            raise RuntimeError(
                "Unexpected residual schema: "
                f"{list(residual.columns)}"
            )

        if len(residual) != 107:
            raise RuntimeError(
                "Expected 107 independent-deletion "
                f"residual pairs; found {len(residual)}"
            )

                                                              
                                                          
                             
                                                              

        authoritative = (
            root
            / "54_figures"
            / "02_curated_inputs"
            / "Fig4_6"
            / "Fig4_6_STRESS_GRANULE_RNP_residual_pairs.tsv"
        )

        if not authoritative.is_file():
            raise FileNotFoundError(
                authoritative
            )

        frozen = pd.read_csv(
            authoritative,
            sep="\t",
        )

        pd.testing.assert_frame_equal(
            frozen,
            residual,
            check_dtype=False,
            check_exact=False,
            rtol=0,
            atol=1e-12,
        )

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copy2(
            residual_file,
            output,
        )

        print()
        print(
            "EXACT 49c REPLAY SUCCESS"
        )
        print(
            "Residual pairs:",
            len(residual),
        )
        print(
            "SHA256:",
            sha256(output),
        )
        print(
            "Output:",
            output,
        )


if __name__ == "__main__":
    main()
