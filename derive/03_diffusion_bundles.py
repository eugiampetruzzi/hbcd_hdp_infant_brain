#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import write_python_provenance, write_tsv


def main() -> None:
    s = get_settings()
    root = s.derivatives / "dwi_qsirecon"
    paths = sorted(root.rglob("*_ses-V02_space-ACPC_bundles-DSIStudio_scalarstats.tsv"))
    if not paths:
        raise RuntimeError(f"No bundle scalar-statistics files beneath {root}")
    for metric in ("fa", "md", "mkt"):
        rows = []
        for path in paths:
            frame = pd.read_csv(path, sep="\t")
            frame = frame.loc[frame["variable_name"].eq(f"dki_{metric}")].copy()
            if frame.empty:
                continue
            participant = str(frame["subject_id"].iloc[0])
            wide = frame.set_index("bundle")["masked_median"].replace([np.inf, -np.inf], np.nan)
            row = {"participant_id": participant, "session_id": "ses-V02"}
            row.update({f"{metric}_{bundle}": value for bundle, value in wide.items()})
            rows.append(row)
        result = pd.DataFrame(rows).sort_values("participant_id")
        bundle_columns = [name for name in result if name.startswith(metric + "_")]
        if len(bundle_columns) != 67:
            raise ValueError(f"Expected 67 {metric} bundles; found {len(bundle_columns)}")
        output = s.output / f"outcomes/diffusion/v02_bundle_{metric}_masked_median.tsv"
        write_tsv(result, output)
        print(f"Wrote {output}: {len(result)} participants, 67 bundles")
    write_python_provenance(
        s.output / "provenance/03_diffusion_bundles.json",
        paths,
        {"source_files": len(paths), "metrics": ["fa", "md", "mkt"]},
    )


if __name__ == "__main__":
    main()
