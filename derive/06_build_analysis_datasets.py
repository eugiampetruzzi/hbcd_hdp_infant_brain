#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import read_table, require_unique, write_python_provenance, write_tsv


def choose(generated: Path, frozen: Path) -> Path:
    return generated if generated.is_file() else frozen


def restricted(frame: pd.DataFrame, complications: pd.DataFrame) -> pd.DataFrame:
    result = frame.merge(complications, on="participant_id", how="inner", validate="one_to_one")
    keep = result["primary_hdp"].eq(1) | (
        result["primary_hdp"].eq(0) & result["none_of_these"].eq(1)
    )
    return result.loc[keep].copy()


def main() -> None:
    s = get_settings()
    md, out = s.manuscript_data, s.output
    complications_path = md / "covariates/v02_pregnancy_complications.tsv"
    complications = read_table(complications_path)[[
        "participant_id", "other_complication", "none_of_these",
    ]]
    require_unique(complications, ["participant_id"], complications_path)
    inputs = [complications_path]

    dwi_cov_path = choose(md / "analysis_datasets/v02_dwi_primary_covariates.tsv", md / "covariates/v02_dwi_primary_covariates.tsv")
    dwi_out_path = choose(out / "outcomes/diffusion/v02_global_cerebral_wm.tsv", md / "outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv")
    dwi = read_table(dwi_cov_path)
    dwi = dwi.loc[dwi["primary_model_complete"].astype(bool)].copy()
    dwi_out = read_table(dwi_out_path)[["participant_id", "global_median_md", "global_median_fa", "global_median_mkt"]]
    dwi = restricted(dwi.drop(columns=["global_median_mkt"], errors="ignore").merge(dwi_out, on="participant_id", how="inner", validate="one_to_one"), complications)
    write_tsv(dwi, out / "analysis_datasets/v02_diffusion_primary.tsv")
    inputs += [dwi_cov_path, dwi_out_path]

    structural_path = md / "analysis_datasets/v02_global_cortical_morphometry_model_data.tsv"
    structural = restricted(read_table(structural_path), complications)
    write_tsv(structural, out / "analysis_datasets/v02_structural_primary.tsv")
    inputs.append(structural_path)

    functional_path = md / "analysis_datasets/v02_global_functional_hdp_data.tsv"
    functional = restricted(read_table(functional_path), complications)
    functional = functional.loc[functional["retained_duration_minutes"].ge(5)].copy()
    write_tsv(functional, out / "analysis_datasets/v02_functional_primary.tsv")
    inputs.append(functional_path)

    observed = {
        "diffusion": [len(dwi), int(dwi.primary_hdp.sum())],
        "structural": [len(structural), int(structural.primary_hdp.sum())],
        "functional": [len(functional), int(functional.primary_hdp.sum())],
    }
    expected = {"diffusion": [1229, 344], "structural": [1100, 302], "functional": [864, 225]}
    if observed != expected:
        raise ValueError(f"Analysis sample mismatch: observed {observed}, expected {expected}")
    write_python_provenance(
        out / "provenance/06_build_analysis_datasets.json", inputs, {"counts": observed}
    )
    print(f"Validated analysis samples: {observed}")


if __name__ == "__main__":
    main()
