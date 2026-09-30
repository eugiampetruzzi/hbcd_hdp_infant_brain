#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from hbcd_hdp.config import get_settings
from hbcd_hdp.imaging import vertex_areas, weighted_mean
from hbcd_hdp.io import write_python_provenance, write_tsv


HASH = "xcp_d-0f306a2f+0ef9c88a"


def one(paths: list[Path], description: str) -> Path:
    if len(paths) != 1:
        raise ValueError(f"Expected one {description}; found {len(paths)}")
    return paths[0]


def cifti_surface_values(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    image = nib.load(path)
    data = np.asarray(image.dataobj).squeeze()
    axis = image.header.get_axis(1)
    result = {}
    for structure, slc, model in axis.iter_structures():
        if "CORTEX_LEFT" in structure:
            result["L"] = (data[slc], np.asarray(model.vertex, dtype=int))
        elif "CORTEX_RIGHT" in structure:
            result["R"] = (data[slc], np.asarray(model.vertex, dtype=int))
    if set(result) != {"L", "R"}:
        raise ValueError(f"Missing cortical CIFTI structures: {path}")
    return result


def surface_area(path: Path) -> np.ndarray:
    image = nib.load(path)
    coordinates = np.asarray(image.darrays[0].data)
    triangles = np.asarray(image.darrays[1].data)
    return vertex_areas(coordinates, triangles)


def export_parcels(root: Path, output_root: Path) -> list[Path]:
    outputs = []
    for metric in ("thickness", "sulc"):
        paths = sorted(root.glob(
            f"sub-*/ses-V02/anat/*space-fsLR_seg-MyersLabonte_stat-mean_desc-{metric}_morph.tsv"
        ))
        rows = []
        expected = None
        for path in paths:
            frame = pd.read_csv(path, sep="\t")
            if len(frame) != 1:
                raise ValueError(f"Expected one row: {path}")
            columns = [re.sub(r"_hash-.*$", "", item) for item in frame.columns]
            if expected is None:
                expected = columns
            elif columns != expected:
                raise ValueError(f"Parcel order mismatch: {path}")
            participant = next(item for item in path.parts if item.startswith("sub-"))
            run = re.search(r"_run-([^_]+)_", path.name)
            rows.append([participant, "ses-V02", run.group(1) if run else "", *frame.iloc[0]])
        result = pd.DataFrame(rows, columns=["participant_id", "session_id", "source_run", *(expected or [])])
        output = output_root / f"v02_myerslabonte_{metric}.tsv"
        write_tsv(result, output)
        outputs.append(output)
    return outputs


def main() -> None:
    s = get_settings()
    root = s.derivatives / "structure_function/hbcd/derivatives" / HASH
    if not root.exists():
        alternate = s.derivatives / "structure-function coupling/hbcd/derivatives" / HASH
        root = alternate if alternate.exists() else root
    output_root = s.output / "outcomes/structural"
    parcel_outputs = export_parcels(root, output_root)
    rows, inputs = [], []
    for subject_dir in sorted(root.glob("sub-*")):
        anat = subject_dir / "ses-V02/anat"
        if not anat.exists():
            continue
        thickness = list(anat.glob("*thickness.dscalar.nii"))
        sulc = list(anat.glob("*sulc.dscalar.nii"))
        left = list(anat.glob("*hemi-L*midthickness.surf.gii"))
        right = list(anat.glob("*hemi-R*midthickness.surf.gii"))
        if not (thickness and sulc and left and right):
            continue
        paths = [one(thickness, "thickness map"), one(sulc, "sulc map"), one(left, "left surface"), one(right, "right surface")]
        inputs.extend(paths)
        areas = {"L": surface_area(paths[2]), "R": surface_area(paths[3])}
        measures = {"thickness": cifti_surface_values(paths[0]), "sulc": cifti_surface_values(paths[1])}
        row = {"participant_id": subject_dir.name, "session_id": "ses-V02"}
        for metric in ("thickness", "sulc"):
            values, weights = [], []
            for hemisphere in ("L", "R"):
                data, vertices = measures[metric][hemisphere]
                values.append(data)
                weights.append(areas[hemisphere][vertices])
            values = np.concatenate(values)
            weights = np.concatenate(weights)
            if metric == "sulc":
                row["global_area_weighted_sulc"] = weighted_mean(values, weights)
                values = np.abs(values)
                key = "global_area_weighted_absolute_sulc"
            else:
                key = "global_area_weighted_thickness"
            row[key] = weighted_mean(values, weights)
            row[f"{metric}_valid_vertex_count"] = int(np.isfinite(values).sum())
        row["cortical_surface_area_mm2"] = float(sum(item.sum() for item in areas.values()))
        rows.append(row)
    output = output_root / "v02_global_cortical_morphometry.tsv"
    write_tsv(pd.DataFrame(rows), output)
    write_python_provenance(
        s.output / "provenance/04_structural_morphometry.json",
        inputs + parcel_outputs,
        {"participants": len(rows), "processing_hash": HASH},
    )
    print(f"Wrote {output}: {len(rows)} participants")


if __name__ == "__main__":
    main()
