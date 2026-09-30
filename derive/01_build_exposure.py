#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import read_table, require_unique, write_python_provenance, write_tsv


PREG_PREFIX = "pex_bm_healthv2_preg__compl_001___"
HISTORY_FIELD = "pex_bm_health_preg__chroncond_001___1"


def main() -> None:
    s = get_settings()
    pregnancy_path = s.root / "pregnancy_exposure/pex_bm_healthv2_preg.tsv"
    history_path = s.root / "pregnancy_exposure/pex_bm_health_preg__chroncond.tsv"
    pregnancy = read_table(pregnancy_path)
    history = read_table(history_path)
    questionnaire_fields = [PREG_PREFIX + str(index) for index in range(1, 13)]
    for field in questionnaire_fields:
        pregnancy[field] = pd.to_numeric(pregnancy[field], errors="coerce")
    history[HISTORY_FIELD] = pd.to_numeric(history[HISTORY_FIELD], errors="coerce")
    pregnancy = pregnancy.loc[pregnancy["session_id"].eq("ses-V02")].copy()
    history = history.loc[history["session_id"].eq("ses-V01")].copy()
    require_unique(pregnancy, ["participant_id"], pregnancy_path)
    require_unique(history, ["participant_id"], history_path)
    history = history[["participant_id", HISTORY_FIELD]].rename(
        columns={HISTORY_FIELD: "history_hypertension"}
    )
    data = pregnancy.merge(history, on="participant_id", how="left")
    gestational = data[PREG_PREFIX + "2"]
    preeclampsia = data[PREG_PREFIX + "3"]
    other_items = [1, 4, 5, 6, 7, 8, 9, 10, 11]
    data["pregnancy_hdp"] = gestational.eq(1) | preeclampsia.eq(1)
    data["other_complication"] = data[
        [PREG_PREFIX + str(index) for index in other_items]
    ].eq(1).any(axis=1)
    data["none_of_these"] = data[PREG_PREFIX + "12"].eq(1)
    data["hypertension_data_complete"] = (
        data["history_hypertension"].isin([0, 1])
        & gestational.isin([0, 1])
        & preeclampsia.isin([0, 1])
    )
    data["group"] = "outside_health_cohort"
    eligible = data["hypertension_data_complete"]
    data.loc[eligible & data["pregnancy_hdp"], "group"] = "HDP"
    control = (
        eligible
        & data["history_hypertension"].eq(0)
        & ~data["pregnancy_hdp"]
        & ~data["other_complication"]
        & data["none_of_these"]
    )
    data.loc[control, "group"] = "control"
    cohort = data.loc[data["group"].isin(["HDP", "control"]), [
        "participant_id", "group", "pregnancy_hdp", "history_hypertension",
        "other_complication", "none_of_these", "hypertension_data_complete",
    ]].copy()
    counts = cohort["group"].value_counts().to_dict()
    expected = {"HDP": 495, "control": 1258}
    if counts != expected or len(cohort) != 1753:
        raise ValueError(f"Cohort mismatch: observed {counts}, expected {expected}")
    output = s.output / "cohort/v02_hdp_control_cohort.tsv"
    write_tsv(cohort, output)
    write_python_provenance(
        s.output / "provenance/01_build_exposure.json",
        [pregnancy_path, history_path],
        {"counts": counts, "total": len(cohort)},
    )
    print(f"Wrote {output}: {len(cohort)} participants ({counts})")


if __name__ == "__main__":
    main()
