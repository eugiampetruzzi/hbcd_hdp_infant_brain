#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import write_python_provenance, write_tsv


def p_text(value: float) -> str:
    if not np.isfinite(value):
        return ""
    return "<.001" if value < 0.001 else f"{value:.3f}".lstrip("0")


def format_model_table(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for column in ("beta", "se", "ci_low", "ci_high"):
        if column in result:
            result[column] = result[column].map(lambda x: "" if pd.isna(x) else f"{x:.3f}")
    for column in ("p", "q", "indirect_p", "indirect_q"):
        if column in result:
            result[column] = result[column].map(p_text)
    return result


def main() -> None:
    settings = get_settings()
    source = settings.output / "results"
    destination = settings.output / "tables"
    destination.mkdir(parents=True, exist_ok=True)
    inputs: list[Path] = []
    for path in sorted(source.glob("table_*.tsv")):
        frame = pd.read_csv(path, sep="\t")
        output = destination / path.name
        write_tsv(format_model_table(frame), output)
        inputs.append(path)
    if not inputs:
        raise SystemExit(f"No analysis tables found beneath {source}")
    write_python_provenance(
        settings.output / "provenance/01_build_tables.json",
        inputs,
        {"formatted_tables": len(inputs)},
    )
    print(f"Formatted {len(inputs)} tables in {destination}")


if __name__ == "__main__":
    main()
