"""Small, entirely invented files with analytically checkable imaging values."""

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.cifti2 import cifti2_axes as axes


def write_tsv(path, data, index=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data).to_csv(path, sep="\t", index=index, na_rep="n/a")
    return path


def volume(path, data, affine=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    image = nib.Nifti1Image(
        np.asarray(data, dtype=np.float32), np.eye(4) if affine is None else affine
    )
    nib.save(image, path)
    return path


def surface(path, scale=1):
    # Two unequal triangles (areas 1 and 2), plus one isolated vertex.
    coords = np.array(
        [[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 2], [3, 3, 3]], np.float32
    )
    faces = np.array([[0, 1, 2], [0, 1, 3]], np.int32)
    image = nib.gifti.GiftiImage(
        darrays=[
            nib.gifti.GiftiDataArray(coords * scale, intent="NIFTI_INTENT_POINTSET"),
            nib.gifti.GiftiDataArray(faces, intent="NIFTI_INTENT_TRIANGLE"),
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(image, path)
    return path


def brain_axis(nvertices=5):
    left = axes.BrainModelAxis.from_surface([0, 1, 2, 3, 4], nvertices, "CortexLeft")
    right = axes.BrainModelAxis.from_surface([0, 1, 2, 3, 4], nvertices, "CortexRight")
    return left + right


def scalar(path, values, axis=None):
    header = nib.Cifti2Header.from_axes(
        (axes.ScalarAxis(["synthetic"]), brain_axis() if axis is None else axis)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Cifti2Image(np.array([values], np.float32), header), path)
    return path


def structural_files(root):
    anat = root / "sub-test/ses-V02/anat"
    prefix = "sub-test_ses-V02_run-1"
    thickness = scalar(
        anat / f"{prefix}_space-fsLR_den-91k_thickness.dscalar.nii",
        [1, 3, 0, np.nan, 99, 2, 4, 0, np.inf, 99],
    )
    scalar(
        anat / f"{prefix}_space-fsLR_den-91k_sulc.dscalar.nii",
        [-2, 0, 6, np.nan, 99, -4, 0, 12, np.inf, 99],
    )
    for hemi in ("L", "R"):
        surface(
            anat
            / f"{prefix}_hemi-{hemi}_space-fsLR_den-32k_desc-hcp_midthickness.surf.gii",
            scale=1 if hemi == "L" else 2,
        )
    for metric in ("thickness", "sulc"):
        write_tsv(
            anat / f"{prefix}_seg-MyersLabonte_stat-mean_desc-{metric}_morph.tsv",
            {"L_Test_1_hash-abc": [2.0], "R_Test_2_hash-abc": [np.inf]},
        )
    return thickness


def diffusion_files(module):
    dwi = module.QSIRECON / "a/qsirecon-DIPYDKI/sub-test/ses-V02/dwi"
    values = np.array([1, 3, 0, np.nan, np.inf, 5, 100, 100]).reshape(2, 2, 2)
    maps = {
        metric: volume(
            dwi / f"sub-test_ses-V02_space-ACPC_{entity}_dwimap.nii.gz", values * scale
        )
        for (metric, entity), scale in zip(module.MAPS.items(), [1, 0.1, 0.001])
    }
    segmentation = volume(
        module.QSIPREP
        / "sub-test/ses-V02/anat/sub-test_ses-V02_space-ACPC_desc-aseg_dseg.nii.gz",
        np.array([2, 41, 2, 41, 2, 41, 3, 3]).reshape(2, 2, 2),
    )
    write_tsv(
        dwi / "sub-test_ses-V02_space-ACPC_bundles-DSIStudio_scalarstats.tsv",
        {
            "subject_id": ["sub-test"] * 6,
            "bundle": ["A", "B"] * 3,
            "variable_name": [f"dki_{m}" for m in module.MAPS for _ in range(2)],
            "masked_median": [1, np.inf, 0.2, 0.3, 0.001, 0.002],
        },
    )
    return maps, segmentation


def functional_files(root, participant="sub-test", incomplete=False):
    folder = root / participant / "ses-V02/func"
    prefix = f"{participant}_ses-V02_task-rest"
    path = (
        folder
        / f"{prefix}_space-fsLR_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv"
    )
    # Positive within-system edges; equal negative between-system edges.
    matrix = np.array(
        [
            [1, 0.5, -0.25, -0.25],
            [0.5, 1, -0.25, -0.25],
            [-0.25, -0.25, 1, 0.5],
            [-0.25, -0.25, 0.5, 1],
        ]
    )
    if incomplete:
        matrix[0, 2] = matrix[2, 0] = np.nan
    labels = ["L_A_1_hash-abc", "R_A_2_hash-abc", "L_B_1_hash-abc", "R_B_2_hash-abc"]
    write_tsv(path, pd.DataFrame(matrix, index=labels, columns=labels), index=True)
    write_tsv(
        folder / f"{prefix}_motion.tsv",
        {"framewise_displacement": [0.1, 10, 0.3, np.nan]},
    )
    write_tsv(
        folder / f"{prefix}_outliers.tsv", {"framewise_displacement": [0, 1, 0, 0]}
    )
    return path
