from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import load_dotenv


@dataclass(frozen=True)
class Settings:
    root: Path
    manuscript_data: Path
    derivatives: Path
    output: Path
    systems: Path
    bootstrap_seed: int
    bootstrap_reps: int


def get_settings() -> Settings:
    load_dotenv()
    root = Path(os.environ.get("HBCD_ROOT", "/projects/iang/eanlp/eu_hbcd-release-2.1"))
    manuscript = Path(
        os.environ.get("HBCD_MANUSCRIPT_DATA", root / "imaging/manuscript_data")
    )
    derivatives = Path(os.environ.get("HBCD_DERIVATIVES", root / "imaging"))
    output = Path(
        os.environ.get(
            "HBCD_OUTPUT", manuscript / "reproducible_hdp_results"
        )
    )
    systems = Path(
        os.environ.get(
            "HBCD_MYERSLABONTE_LABELS",
            manuscript / "config/myerslabonte_parcel_systems.tsv",
        )
    )
    return Settings(
        root=root,
        manuscript_data=manuscript,
        derivatives=derivatives,
        output=output,
        systems=systems,
        bootstrap_seed=int(os.environ.get("HBCD_BOOTSTRAP_SEED", "20260930")),
        bootstrap_reps=int(os.environ.get("HBCD_BOOTSTRAP_REPS", "5000")),
    )
