"""Competing release columns must not change measurement or subject identity."""

import shutil

import pandas as pd
import pytest

from header_fixtures import release_frame
from synthetic import diffusion_files, write_tsv


def test_bundle_statistics_use_masked_median_by_subject_bundle_and_metric(modules):
    m = modules("prepare/03_diffusion.py")
    maps, segmentation = diffusion_files(m)
    # A second complete synthetic scan, with reversed bundle measurements.
    for path in [*maps.values(), segmentation]:
        other = type(path)(str(path).replace("sub-test", "sub-other"))
        other.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, other)
    expected = {"fa": (.2, .3), "md": (.001, .002), "mkt": (.5, .6)}
    for subject in ["sub-test", "sub-other"]:
        rows = []
        for metric, values in expected.items():
            for bundle, value in zip(["A", "B"], values if subject == "sub-test" else values[::-1]):
                rows.append(dict(subject_id=subject, session_id="ses-V02", bundle=bundle,
                                 variable_name=f"dki_{metric}", masked_median=value,
                                 median=.7, mean=.8, masked_mean=.85,
                                 weighted_mean=.9, masked_weighted_mean=.95,
                                 zero_proportion=.99, source_file="synthetic"))
        rows.append(dict(subject_id=subject, bundle="A", variable_name="unselected", masked_median=999))
        frame = release_frame("bundle_statistics", rows)
        path = maps["fa"].parent / "sub-test_ses-V02_space-ACPC_bundles-DSIStudio_scalarstats.tsv"
        path = type(path)(str(path).replace("sub-test", subject))
        write_tsv(path, frame.iloc[::-1, ::-1])
    m.main(expected_count=2)
    for metric, values in expected.items():
        result = pd.read_csv(m.DATA / f"outcomes/diffusion/v02_bundle_{metric}_masked_median.tsv", sep="\t").set_index("participant_id")
        assert set(result.index) == {"sub-test", "sub-other"}
        for subject, reference in [("sub-test", values), ("sub-other", values[::-1])]:
            assert result.loc[subject, [f"{metric}_A", f"{metric}_B"]].astype(float).tolist() == pytest.approx(reference)


def test_scanner_variants_preserve_identity_without_substituting_site(modules):
    m = modules("prepare/02_scanners.py")
    for variant, subject in enumerate(["sub-a", "sub-b", "sub-unknown"]):
        rows = dict(filename=["dwi/a_dwi.nii.gz", "anat/a_T2w.nii.gz", "func/a_bold.nii.gz"])
        if variant < 2:
            rows.update(scanner_serial_number=["S1"] * 3, scanner_software=["v1"] * 3,
                        scanner_manufacturer=["Siemens"] * 3, scanner_model=["synthetic"] * 3)
        frame = release_frame("scanner_metadata", rows, variant)
        # Same site must not manufacture scanner identity for the unidentified variant.
        frame["site"] = "same-site"
        path = m.release(f"imaging/qc_mriqc/hbcd/rawdata/{subject}/ses-V02/{subject}_ses-V02_scans.tsv")
        write_tsv(path, frame)
    m.main()
    result = pd.read_csv(m.DATA / "covariates/v02_scanner_metadata.tsv", sep="\t").set_index("participant_id")
    for prefix in ["dwi", "t2w", "bold"]:
        assert result.loc[["sub-a", "sub-b"], f"{prefix}_scanner_serial_number"].tolist() == ["S1", "S1"]
        assert result.loc[["sub-a", "sub-b"], f"{prefix}_scanner_software"].tolist() == ["v1", "v1"]
        assert result.loc[["sub-a", "sub-b"], f"{prefix}_unique_scanners"].tolist() == [1, 1]
    if "sub-unknown" in result.index:
        assert pd.isna(result.loc["sub-unknown", "dwi_scanner_serial_number"])
