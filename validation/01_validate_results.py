#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import write_tsv


REPO = Path(__file__).resolve().parents[1]


def record(rows: list[dict], name: str, observed, expected, passed: bool) -> None:
    rows.append({
        "check": name,
        "observed": observed,
        "expected": expected,
        "status": "PASS" if passed else "FAIL",
    })


def main() -> None:
    settings = get_settings()
    expected_counts = json.loads((REPO / "config/expected_counts.json").read_text())
    expected_results = json.loads((REPO / "config/expected_results.json").read_text())
    rows: list[dict] = []

    for modality in ("diffusion", "structural", "functional"):
        path = settings.output / "analysis_datasets" / f"v02_{modality}_primary.tsv"
        if path.exists():
            data = pd.read_csv(path, sep="\t")
            expected = expected_counts[modality]
            record(rows, f"{modality} total", len(data), expected["total"], len(data) == expected["total"])
            n_hdp = int(data["primary_hdp"].sum())
            record(rows, f"{modality} HDP", n_hdp, expected["HDP"], n_hdp == expected["HDP"])
            record(rows, f"{modality} controls", len(data) - n_hdp, expected["control"], len(data) - n_hdp == expected["control"])

    primary_path = settings.output / "results/table_primary_global.tsv"
    primary = pd.read_csv(primary_path, sep="\t")
    record(rows, "primary outcome count", len(primary), 7, len(primary) == 7)
    record(rows, "primary estimates finite", bool(np.isfinite(primary[["beta", "se", "p", "q"]]).all().all()), True,
           bool(np.isfinite(primary[["beta", "se", "p", "q"]]).all().all()))

    sensitivity_path = settings.output / "results/table_no_complications.tsv"
    sensitivity = pd.read_csv(sensitivity_path, sep="\t").set_index("outcome")
    tolerance = expected_results["tolerance"]
    for outcome, expected in expected_results["no_complications"].items():
        for statistic in ("beta", "p"):
            observed = float(sensitivity.loc[outcome, statistic])
            limit = tolerance[statistic]
            record(rows, f"sensitivity {outcome} {statistic}", observed, expected[statistic],
                   abs(observed - expected[statistic]) <= limit)

    report = pd.DataFrame(rows)
    destination = settings.output / "validation/validation_report.tsv"
    write_tsv(report, destination)
    failures = report[report.status == "FAIL"]
    if len(failures):
        raise SystemExit("Validation failed:\n" + failures.to_string(index=False))
    print(f"All {len(report)} validation checks passed; wrote {destination}")


if __name__ == "__main__":
    main()
