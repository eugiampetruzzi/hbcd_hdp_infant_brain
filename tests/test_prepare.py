import runpy

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from synthetic import (
    brain_axis,
    diffusion_files,
    functional_files,
    scalar,
    structural_files,
    surface,
    volume,
    write_tsv,
)


@pytest.mark.parametrize(
    "field,value,expected",
    [
        ("scanner_manufacturer", " Siemens Healthineers ", "Siemens"),
        ("scanner_manufacturer", "GE MEDICAL SYSTEMS", "GE"),
        ("scanner_manufacturer", "Philips Medical Systems", "Philips"),
        ("scanner_manufacturer", "Other", "Other"),
        ("scanner_serial_number", "123.0", "123"),
        ("scanner_serial_number", "ABC.0", "ABC.0"),
        ("scanner_software", "v1\\\\v2", "v1"),
        ("scanner_software", None, ""),
    ],
)
def test_scanner_cleaning(modules, field, value, expected):
    assert modules("prepare/02_scanners.py").clean(field, value) == expected


def test_scanner_pipeline_counts_distinct_acquisitions(modules):
    m = modules("prepare/02_scanners.py")
    path = m.release(
        "imaging/qc_mriqc/hbcd/rawdata/sub-test/ses-V02/sub-test_ses-V02_scans.tsv"
    )
    write_tsv(
        path,
        {
            "filename": [
                "dwi/a_dwi.nii.gz",
                "dwi/b_dwi.nii.gz",
                "anat/a_T2w.nii.gz",
                "anat/b_T2w.nii.gz",
                "func/a_bold.nii.gz",
                "anat/a_T1w.nii.gz",
            ],
            "scanner_serial_number": ["123.0", "123", "1", "2", "", "3"],
            "scanner_software": ["v1"] * 4 + ["", "v1"],
        },
    )
    m.main()
    row = pd.read_csv(m.DATA / "covariates/v02_scanner_metadata.tsv", sep="\t").iloc[0]
    assert row.participant_id == "sub-test"
    assert row.dwi_unique_scanners == 1
    assert row.dwi_scanner_serial_number == 123
    assert row.t2w_unique_scanners == 2
    assert pd.isna(row.t2w_scanner_serial_number)
    assert row.bold_unique_scanners == 0


def test_scanner_script_entry_point_writes_metadata(modules):
    m = modules("prepare/02_scanners.py")
    path = m.release(
        "imaging/qc_mriqc/hbcd/rawdata/sub-cli/ses-V02/sub-cli_ses-V02_scans.tsv"
    )
    write_tsv(
        path,
        {
            "filename": ["func/a_bold.nii.gz"],
            "scanner_serial_number": [456],
            "scanner_software": ["v2"],
        },
    )
    runpy.run_path(m.__file__, run_name="__main__")
    data = pd.read_csv(m.DATA / "covariates/v02_scanner_metadata.tsv", sep="\t")
    assert data.participant_id.tolist() == ["sub-cli"]
    assert data.bold_scanner_serial_number.tolist() == [456]


@pytest.fixture
def diffusion(modules):
    m = modules("prepare/03_diffusion.py")
    maps, segmentation = diffusion_files(m)
    return m, maps, segmentation


def test_white_matter_medians_exclude_background_zero_and_nonfinite(diffusion):
    m, maps, _ = diffusion
    assert m.white_matter_medians("sub-test", maps) == pytest.approx(
        {"global_median_mkt": 3, "global_median_fa": 0.3, "global_median_md": 0.003}
    )


def test_diffusion_requires_both_hemispheres(diffusion):
    m, maps, segmentation = diffusion
    volume(segmentation, np.full((2, 2, 2), 2))
    assert all(np.isnan(v) for v in m.white_matter_medians("sub-test", maps).values())


def test_diffusion_resamples_labels_in_world_coordinates(diffusion):
    m, maps, segmentation = diffusion
    affine = np.eye(4)
    affine[0, 3] = -1
    labels = np.zeros((3, 2, 2))
    labels[1, 0, 0] = 2
    labels[2, 0, 1] = 41
    volume(segmentation, labels, affine)
    assert m.white_matter_medians("sub-test", maps)["global_median_mkt"] == 3


def test_diffusion_no_valid_values_returns_nan(diffusion):
    m, maps, _ = diffusion
    volume(maps["fa"], np.zeros((2, 2, 2)))
    assert np.isnan(m.white_matter_medians("sub-test", maps)["global_median_fa"])


def test_diffusion_duplicate_paths_choose_last(diffusion):
    m, maps, _ = diffusion
    path = m.QSIRECON / "z/qsirecon-DIPYDKI/sub-test/ses-V02/dwi" / maps["fa"].name
    volume(path, np.ones((2, 2, 2)))
    assert m.by_participant("model-tensor_param-fa_dwimap.nii.gz") == {"sub-test": path}


def test_diffusion_pipeline_small_release(diffusion):
    m, _, _ = diffusion
    volume(
        m.QSIRECON
        / "a/qsirecon-DIPYDKI/sub-incomplete/ses-V02/dwi/sub-incomplete_ses-V02_space-ACPC_model-tensor_param-fa_dwimap.nii.gz",
        np.ones((2, 2, 2)),
    )
    m.main(expected_count=1)
    out = m.DATA / "outcomes/diffusion"
    row = pd.read_csv(out / "v02_global_cerebral_wm_fa_md.tsv", sep="\t").iloc[0]
    assert row.global_median_mkt == 3
    assert row.session_id == "ses-V02"
    bundles = pd.read_csv(out / "v02_bundle_mkt_masked_median.tsv", sep="\t")
    assert bundles.mkt_A.tolist() == [1]
    assert bundles.mkt_B.isna().all()


def test_surface_triangle_area_distribution(modules, tmp_path):
    m = modules("prepare/04_structural.py")
    np.testing.assert_allclose(
        m.vertex_areas(surface(tmp_path / "test.surf.gii")), [1, 1, 1 / 3, 2 / 3, 0]
    )


@pytest.fixture
def structural(modules):
    m = modules("prepare/04_structural.py")
    return m, structural_files(m.release(m.XCPD))


def test_cortex_preserves_vertex_indices(structural):
    m, path = structural
    result = m.cortex(path)
    assert set(result) == {"LEFT", "RIGHT"}
    np.testing.assert_array_equal(result["RIGHT"][1], [0, 1, 2, 3, 4])
    np.testing.assert_allclose(result["LEFT"][0], [1, 3, 0, np.nan, 99])


def test_cortex_excludes_subcortical_structures(modules, tmp_path):
    m = modules("prepare/04_structural.py")
    subcortex = nib.cifti2.cifti2_axes.BrainModelAxis.from_mask(
        np.ones((1, 1, 1), bool), "ThalamusLeft", np.eye(4)
    )
    path = scalar(
        tmp_path / "with_volume.dscalar.nii",
        list(range(10)) + [999],
        brain_axis() + subcortex,
    )
    data = m.cortex(path)
    assert set(data) == {"LEFT", "RIGHT"}
    np.testing.assert_array_equal(data["RIGHT"][0], [5, 6, 7, 8, 9])


def test_global_morphometry_weights_hemispheres_and_valid_vertices(structural):
    m, path = structural
    row = m.global_morphometry(path)
    # Left contributes 4/2; right contributes 24/8. Sulc: 4/(4/3), 32/(16/3).
    assert row["global_area_weighted_thickness"] == pytest.approx(2.8)
    assert row["global_area_weighted_absolute_sulc"] == pytest.approx(5.4)
    assert row["participant_id"] == "sub-test"
    assert row["source_run"] == "1"


def test_structural_pipeline_small_release(structural):
    m, _ = structural
    m.main(expected_count=1)
    out = m.DATA / "outcomes/structural"
    global_data = pd.read_csv(out / "v02_global_cortical_morphometry.tsv", sep="\t")
    assert global_data.global_area_weighted_thickness.tolist() == [2.8]
    parcels = pd.read_csv(out / "v02_myerslabonte_thickness.tsv", sep="\t")
    assert list(parcels) == [
        "participant_id",
        "session_id",
        "source_run",
        "L_Test_1",
        "R_Test_2",
    ]
    assert parcels.L_Test_1.tolist() == [2]
    assert parcels.R_Test_2.isna().all()


@pytest.fixture
def functional(modules):
    m = modules("prepare/05_functional.py")
    return m, functional_files(m.release(m.XCPD))


def test_connectivity_reader_removes_hashes_and_self_connections(functional):
    m, path = functional
    matrix, labels = m.read_matrix(path)
    assert labels == ["L_A_1", "R_A_2", "L_B_1", "R_B_2"]
    assert np.isnan(matrix.diagonal()).all()
    assert matrix[0, 1] == 0.5
    assert matrix[0, 2] == -0.25


def test_motion_censors_outliers_but_duration_counts_retained_frames(functional):
    m, path = functional
    assert m.motion(path) == pytest.approx(
        {"mean_fd_post_censoring": 0.2, "retained_duration_minutes": 0.08625}
    )


def test_fd_reader_accepts_hashed_column(modules, tmp_path):
    m = modules("prepare/05_functional.py")
    path = write_tsv(
        tmp_path / "motion.tsv",
        {"irrelevant": [9, 9], "framewise_displacement_hash-abc": [0.2, np.nan]},
    )
    np.testing.assert_allclose(m.framewise_displacement(path), [0.2, np.nan])


def test_signed_network_metrics_have_hand_calculated_values(modules):
    m = modules("prepare/05_functional.py")
    z = np.array([[0, 2, -1, -1], [2, 0, -1, -1], [-1, -1, 0, 2], [-1, -1, 2, 0]])
    result = m.network_metrics(z, np.array(["A", "A", "B", "B"]))
    assert result == pytest.approx(
        {
            "parcel_count": 4,
            "system_count": 2,
            "modularity": 0.75,
            "mean_participation_coefficient": 0,
            "mean_positive_participation": 0,
            "mean_negative_participation": 0,
            "mean_network_connectivity": 0,
        }
    )


def test_participation_measures_connections_across_systems(modules):
    m = modules("prepare/05_functional.py")
    z = np.array([[0, 2, 1, -1], [2, 0, -1, 1], [1, -1, 0, 2], [-1, 1, 2, 0]])
    result = m.network_metrics(z, np.array(["A", "A", "B", "B"]))
    assert result["mean_positive_participation"] == pytest.approx(4 / 9)
    assert result["mean_negative_participation"] == 0
    assert result["mean_participation_coefficient"] == pytest.approx(2 / 9)
    assert result["modularity"] == pytest.approx(7 / 24)
    assert result["mean_network_connectivity"] == pytest.approx(2 / 3)


def test_functional_pipeline_excludes_incomplete_matrices(functional):
    m, path = functional
    functional_files(m.release(m.XCPD), "sub-incomplete", incomplete=True)
    m.main(expected_count=2, expected_parcels=4)
    data = pd.read_csv(
        m.DATA / "outcomes/functional/v02_functional_network_metrics.tsv", sep="\t"
    )
    assert data.participant_id.tolist() == ["sub-test"]
    assert data.parcel_count.tolist() == [4]
    assert data.mean_network_connectivity.iloc[0] == pytest.approx(
        (np.arctanh(0.5) - 2 * np.arctanh(0.25)) / 3
    )
    qc = pd.read_csv(m.DATA / "covariates/v02_functional_qc.tsv", sep="\t")
    assert len(qc) == 2


def test_functional_coverage_boundary_unassigned_and_run_specific_exclusions(modules):
    m = modules("prepare/05_functional.py")
    for i in range(5):
        path = functional_files(m.release(m.XCPD), f"sub-test{i}")
        frame = pd.read_csv(path, sep="\t", index_col=0)
        labels = list(frame) + ["None", "L_Sparse_1"]
        matrix = pd.DataFrame(0.1, index=labels, columns=labels)
        matrix.loc[frame.index, frame.columns] = frame
        if i == 0:
            matrix.iloc[0, :] = np.nan
            matrix.iloc[:, 0] = np.nan
        if i < 2:
            matrix.loc["L_Sparse_1", :] = np.nan
            matrix.loc[:, "L_Sparse_1"] = np.nan
        write_tsv(path, matrix, index=True)
        # These deliberately lack motion files: selecting them would fail.
        for entity in ["run-1", "dir-AP"]:
            extra = path.with_name(
                path.name.replace("_task-rest", f"_{entity}_task-rest")
            )
            extra.write_bytes(path.read_bytes())
    m.main(expected_count=5, expected_parcels=4)
    result = pd.read_csv(
        m.DATA / "outcomes/functional/v02_functional_network_metrics.tsv", sep="\t"
    )
    assert result.participant_id.tolist() == [f"sub-test{i}" for i in range(1, 5)]
    assert result.parcel_count.tolist() == [4] * 4


@pytest.mark.parametrize("script", ["03_diffusion", "04_structural", "05_functional"])
def test_default_release_count_is_still_checked(modules, script):
    with pytest.raises(AssertionError):
        modules(f"prepare/{script}.py").main()
