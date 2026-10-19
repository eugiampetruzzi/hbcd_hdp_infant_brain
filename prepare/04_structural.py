"""Global area-weighted cortical thickness and sulcal depth, and Myers-Labonte parcel means (XCP-D)."""

import re

import nibabel as nib
import numpy as np
import pandas as pd

from config import DATA, XCPD, release

ANAT = "sub-*/ses-V02/anat/"


def vertex_areas(path):
    surface = nib.load(path)
    coords, faces = surface.darrays[0].data.astype(float), surface.darrays[
        1
    ].data.astype(int)
    p = coords[faces]
    face_area = 0.5 * np.linalg.norm(
        np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1
    )
    areas = np.zeros(len(coords))
    for corner in range(3):
        np.add.at(areas, faces[:, corner], face_area / 3)
    return areas


def cortex(path):
    img = nib.load(path)
    data = np.asarray(img.dataobj)[0]
    hemis = {
        "CIFTI_STRUCTURE_CORTEX_LEFT": "LEFT",
        "CIFTI_STRUCTURE_CORTEX_RIGHT": "RIGHT",
    }
    return {
        hemis[name]: (data[slc], model.vertex)
        for name, slc, model in img.header.get_axis(1).iter_structures()
        if name in hemis
    }


def global_morphometry(thickness_path):
    prefix = thickness_path.name.replace(
        "_space-fsLR_den-91k_thickness.dscalar.nii", ""
    )
    anat = thickness_path.parent
    thickness = cortex(thickness_path)
    sulc = cortex(anat / f"{prefix}_space-fsLR_den-91k_sulc.dscalar.nii")
    sums = np.zeros(4)
    for hemi in ("LEFT", "RIGHT"):
        t, vertices = thickness[hemi]
        s = sulc[hemi][0]
        surface = (
            anat
            / f"{prefix}_hemi-{hemi[0]}_space-fsLR_den-32k_desc-hcp_midthickness.surf.gii"
        )
        w = vertex_areas(surface)[vertices]
        # Vertices with zero or nonfinite values (medial wall) are left out.
        t_ok = np.isfinite(t) & (t != 0) & (w > 0)
        s_ok = np.isfinite(s) & (s != 0) & (w > 0)
        sums += [
            np.sum(t[t_ok] * w[t_ok]),
            np.sum(w[t_ok]),
            np.sum(np.abs(s[s_ok]) * w[s_ok]),
            np.sum(w[s_ok]),
        ]
    return {
        "participant_id": prefix.split("_")[0],
        "session_id": "ses-V02",
        "source_run": re.search(r"_run-([^_]+)", prefix).group(1),
        "global_area_weighted_thickness": sums[0] / sums[1],
        "global_area_weighted_absolute_sulc": sums[2] / sums[3],
    }


def main(expected_count=1442):
    """Prepare a release, validating its expected number of structural scans."""
    root = release(XCPD)
    paths = sorted(root.glob(ANAT + "*_space-fsLR_den-91k_thickness.dscalar.nii"))
    assert len(paths) == expected_count, len(paths)
    out = DATA / "outcomes/structural"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([global_morphometry(p) for p in paths]).to_csv(
        out / "v02_global_cortical_morphometry.tsv", sep="\t", index=False
    )

    # Parcel columns keep the release order, which matches the atlas label order (figure 4).
    for metric in ("thickness", "sulc"):
        rows = []
        for path in sorted(
            root.glob(ANAT + f"*_seg-MyersLabonte_stat-mean_desc-{metric}_morph.tsv")
        ):
            row = pd.read_csv(path, sep="\t", na_values="n/a")
            row.columns = [re.sub(r"_hash-.*$", "", c) for c in row.columns]
            run = re.search(r"_run-([^_]+)_", path.name)
            row.insert(0, "source_run", run.group(1) if run else "")
            row.insert(0, "session_id", "ses-V02")
            row.insert(0, "participant_id", path.name.split("_")[0])
            rows.append(row)
        table = pd.concat(rows).replace([np.inf, -np.inf], np.nan)
        table.to_csv(out / f"v02_myerslabonte_{metric}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
