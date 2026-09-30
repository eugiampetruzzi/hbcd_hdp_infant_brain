from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="," if path.suffix == ".csv" else "\t", low_memory=False)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def require_columns(frame: pd.DataFrame, columns: list[str], source: Path) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"Missing columns in {source}: {missing}")


def require_unique(frame: pd.DataFrame, columns: list[str], source: Path) -> None:
    duplicated = frame.duplicated(columns, keep=False)
    if duplicated.any():
        example = frame.loc[duplicated, columns].head().to_dict("records")
        raise ValueError(f"Duplicate keys in {source}: {example}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_python_provenance(path: Path, inputs: list[Path], extra: dict | None = None) -> None:
    payload = {
        "python": sys.version,
        "platform": platform.platform(),
        "inputs": [
            {"path": str(item), "sha256": sha256(item), "bytes": item.stat().st_size}
            for item in inputs
        ],
    }
    if extra:
        payload.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
