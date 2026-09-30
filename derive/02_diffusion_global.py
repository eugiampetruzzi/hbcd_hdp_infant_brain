#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import write_python_provenance, write_tsv


PATTERNS = {
    "fa": "*_ses-V02_space-ACPC_model-tensor_param-fa_dwimap.nii.gz",
    "md": "*_ses-V02_space-ACPC_model-dki_param-md_dwimap.nii.gz",
    "mkt": "*_ses-V02_space-ACPC_*_param-mkt_dwimap.nii.gz",
    "dseg": "*_ses-V02_space-ACPC_desc-aseg_dseg.nii.gz",
}


def key(path: Path) -> tuple[str, str]:
    match = re.search(r"(sub-\d+)_(ses-V\d+)", path.name)
    if not match:
        raise ValueError(f"Cannot parse participant/session: {path}")
    return match.group(1), match.group(2)


def index(root: Path, pattern: str) -> dict[tuple[str, str], Path]:
    result: dict[tuple[str, str], Path] = {}
    for path in root.rglob(pattern):
        item = key(path)
        if item in result:
            raise ValueError(f"Duplicate {pattern} for {item}")
        result[item] = path
    return result


def main() -> None:
    s = get_settings()
    recon = s.derivatives / "dwi_qsirecon"
    prep = s.derivatives / "dwi_qsiprep"
    maps = {metric: index(recon, pattern) for metric, pattern in PATTERNS.items() if metric != "dseg"}
    maps["dseg"] = index(prep, PATTERNS["dseg"])
    common = set.intersection(*(set(item) for item in maps.values()))
    if not common:
        raise RuntimeError("No common V02 diffusion maps and segmentations found")
    rows = []
    for participant_id, session_id in sorted(common):
        images = {name: nib.load(path) for name, path in ((name, paths[(participant_id, session_id)]) for name, paths in maps.items())}
        reference = images["mkt"]
        dseg = resample_from_to(images["dseg"], reference, order=0)
        white_matter = np.isin(np.asarray(dseg.dataobj), [2, 41])
        row = {"participant_id": participant_id, "session_id": session_id}
        for metric in ("md", "fa", "mkt"):
            values = np.asarray(images[metric].dataobj, dtype=float)[white_matter]
            usable = np.isfinite(values) & (values != 0)
            row[f"global_median_{metric}"] = float(np.median(values[usable])) if usable.any() else np.nan
            row[f"valid_{metric}_voxel_count"] = int(usable.sum())
        row["cerebral_wm_voxel_count"] = int(white_matter.sum())
        rows.append(row)
    output = s.output / "outcomes/diffusion/v02_global_cerebral_wm.tsv"
    write_tsv(pd.DataFrame(rows), output)
    write_python_provenance(
        s.output / "provenance/02_diffusion_global.json",
        [maps[name][item] for name in maps for item in common],
        {"participants": len(rows), "white_matter_labels": [2, 41]},
    )
    print(f"Wrote {output}: {len(rows)} participants")


if __name__ == "__main__":
    main()
