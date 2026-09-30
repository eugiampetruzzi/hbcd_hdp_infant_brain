#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.imaging import fixed_community_modularity, signed_participation
from hbcd_hdp.io import require_columns, write_python_provenance, write_tsv


HASH = "xcp_d-0f306a2f+0ef9c88a"
MATRIX_PATTERN = "*_ses-V02_*task-rest_space-fsLR_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv"


def read_matrix(path: Path) -> tuple[np.ndarray, list[str]]:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    matrix = frame.to_numpy(dtype=float)
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Non-square correlation matrix: {path}")
    if not np.allclose(matrix, matrix.T, atol=1e-8, equal_nan=True):
        raise ValueError(f"Asymmetric correlation matrix: {path}")
    return matrix, list(map(str, frame.columns))


def main() -> None:
    s = get_settings()
    systems = pd.read_csv(s.systems, sep="\t")
    require_columns(systems, ["parcel", "system"], s.systems)
    system_by_parcel = systems.set_index("parcel")["system"]
    root = s.derivatives / "func_xcpd/hbcd/derivatives" / HASH
    paths = sorted(path for path in root.rglob(MATRIX_PATTERN) if "_run-" not in path.name and "_dir-" not in path.name)
    if not paths:
        raise RuntimeError(f"No run-free V02 functional matrices beneath {root}")
    global_rows, system_rows = [], []
    for path in paths:
        matrix, parcels = read_matrix(path)
        missing = sorted(set(parcels) - set(system_by_parcel.index))
        if missing:
            raise ValueError(f"Missing system labels for {len(missing)} parcels; first: {missing[:5]}")
        communities = system_by_parcel.loc[parcels].to_numpy()
        positive, negative = signed_participation(matrix, communities)
        mean_participation = np.nanmean(np.column_stack([positive, negative]), axis=1)
        participant = re.search(r"(sub-\d+)", path.name).group(1)
        global_rows.append({
            "participant_id": participant,
            "session_id": "ses-V02",
            "parcel_count": len(parcels),
            "system_count": len(np.unique(communities)),
            "modularity": fixed_community_modularity(matrix, communities),
            "mean_participation_coefficient": float(np.nanmean(mean_participation)),
            "mean_positive_participation": float(np.nanmean(positive)),
            "mean_negative_participation": float(np.nanmean(negative)),
            "mean_network_connectivity": float(np.nanmean(matrix[np.triu_indices_from(matrix, 1)])),
            "matrix_path": str(path),
        })
        for system in np.unique(communities):
            keep = communities == system
            system_rows.append({
                "participant_id": participant,
                "session_id": "ses-V02",
                "myerslabonte_system": system,
                "participation_coefficient": float(np.nanmean(mean_participation[keep])),
                "positive_participation": float(np.nanmean(positive[keep])),
                "negative_participation": float(np.nanmean(negative[keep])),
            })
    output_root = s.output / "outcomes/functional"
    global_output = output_root / "v02_functional_network_metrics.tsv"
    system_output = output_root / "v02_myerslabonte_system_participation.tsv"
    write_tsv(pd.DataFrame(global_rows), global_output)
    write_tsv(pd.DataFrame(system_rows), system_output)
    write_python_provenance(
        s.output / "provenance/05_functional_networks.json",
        [s.systems, *paths],
        {"participants": len(global_rows), "processing_hash": HASH},
    )
    print(f"Wrote {global_output}: {len(global_rows)} participants")
    print(f"Wrote {system_output}: {len(system_rows)} participant-system rows")


if __name__ == "__main__":
    main()
