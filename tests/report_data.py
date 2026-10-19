"""Invented participant tables for Python reporting, independent of the R pipeline."""

import numpy as np
import pandas as pd

from synthetic import write_tsv


OUTCOMES = ["global_median_md", "global_median_fa", "global_median_mkt",
            "global_area_weighted_thickness", "global_area_weighted_absolute_sulc",
            "modularity", "mean_participation_coefficient"]


def write_report_inputs(module):
    rng = np.random.default_rng(17)
    n = 32
    ids = [f"sub-synthetic{i:03}" for i in range(n)]
    core = pd.DataFrame({"participant_id": ids, "primary_hdp": np.tile([0, 1], n // 2),
                         "ga_delivery": rng.uniform(34, 41, n), "child_sex": rng.integers(0, 2, n),
                         "candidate_age_weeks": rng.uniform(1, 12, n),
                         "chronic_hypertension": rng.integers(0, 2, n)})
    values = pd.DataFrame({"participant_id": ids, **{c: rng.normal(2, .3, n) for c in OUTCOMES}})
    values["raw_neighbor_corr"] = rng.uniform(.8, .99, n)
    values["brainswipes_summary_qc"] = rng.uniform(1, 3, n)
    values["mean_network_connectivity"] = rng.normal(.1, .02, n)
    values["mean_fd_post_censoring"] = rng.uniform(.05, .2, n)
    values["retained_duration_minutes"] = rng.uniform(5, 15, n)
    all_data = core.merge(values, on="participant_id")
    write_tsv(module.DATA / "covariates/v02_dwi_primary_covariates.tsv",
              core.rename(columns={"child_sex": "sex_m"}).assign(raw_neighbor_corr=values.raw_neighbor_corr))
    write_tsv(module.DATA / "analysis_datasets/v02_global_cortical_morphometry_model_data.tsv", all_data)
    write_tsv(module.DATA / "analysis_datasets/v02_global_functional_hdp_data.tsv", all_data)
    for rel in ["outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv",
                "outcomes/functional/v02_functional_network_metrics.tsv"]:
        write_tsv(module.DATA / rel, values.drop(columns=["raw_neighbor_corr", "mean_fd_post_censoring",
                                                        "retained_duration_minutes"]))
    write_tsv(module.DATA / "covariates/v02_functional_qc.tsv",
              values[["participant_id", "mean_fd_post_censoring", "retained_duration_minutes"]])
    person = pd.DataFrame({"participant_id": ids, "site": np.tile(["A", "B"], n // 2),
                           "mother_age": rng.uniform(20, 40, n), "education": np.resize([1, 2, 3, 4], n),
                           "income": np.resize(np.arange(1, 8), n), "race_eth": np.resize(np.arange(5), n)})
    person.loc[0, ["education", "income"]] = np.nan
    person.to_csv(module.DATA / "covariates/hdp_missingness_person_level.csv", index=False)
    write_tsv(module.DATA / "covariates/v02_pregnancy_complications.tsv",
              {"participant_id": ids, "gestational_hypertension": core.primary_hdp,
               "preeclampsia_eclampsia": core.primary_hdp,
               "other_complication": np.resize([0, 0, 1, 1], n)})
    sample = values[["participant_id"] + OUTCOMES].copy()
    for modality in ["diffusion", "structural", "functional"]:
        sample[f"in_{modality}"] = True
    # Each modality has a different excluded infant; all still enter the union.
    for i, modality in enumerate(["diffusion", "structural", "functional"]):
        sample.loc[i, f"in_{modality}"] = False
    module.TABLES.mkdir(parents=True, exist_ok=True)
    sample.to_csv(module.TABLES / "analytic_sample.csv", index=False)
    return sample, core, values


def model_results():
    return pd.DataFrame([
        {"modality": modality, "outcome": outcome, "model_set": model, "beta": .2,
         "ci_lo": .1, "ci_hi": .3, "p": .001, "q": .007, "n": 32, "n_hdp": 16, "n_control": 16}
        for modality, outcome in [("Diffusion MRI", "Mean diffusivity"),
                                  ("Structural MRI", "Cortical thickness"),
                                  ("Functional MRI", "Network modularity")]
        for model in ["primary", "ga_adjusted"]])


def sensitivity_results():
    return pd.DataFrame({"group": ["Imaging-related factors"],
                         "analysis": ["No harmonization, scanner fixed effects"],
                         "n": [32], "n_hdp": [16], "beta": [.2], "ci_lo": [.1], "ci_hi": [.3], "p": [.001]})
