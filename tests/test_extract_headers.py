"""Exercise the Bash entry point against invented release files."""

import json
import os
from pathlib import Path
import subprocess

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/extract_headers.sh"
ANAT = "imaging/structure_function/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"
FUNC = "imaging/func_xcpd/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"


@pytest.fixture
def release_headers(tmp_path):
    root = tmp_path / "release with spaces"
    root.mkdir()
    return root


def run_extractor(root):
    env = os.environ.copy()
    env.pop("HBCD_RELEASE", None)
    env["datadir"] = "/unused-legacy-root"
    env["HBCD_HDP_DATA"] = "/unused-derived-root"
    if root is not None:
        env["HBCD_RELEASE"] = str(root)
    return subprocess.run(
        ["bash", str(SCRIPT)], env=env, capture_output=True, text=True
    )


def put(root, relative, content):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def test_extracts_ordered_distinct_headers_without_records_or_paths(release_headers):
    root = release_headers
    for participant, header in [
        ("sub-private001", "R_A_2_hash-ab\tL_A_1_hash-cd"),
        ("sub-private002", "R_A_2_hash-ab\tL_A_1_hash-cd"),
        ("sub-private003", "L_A_1_hash-cd\tR_A_2_hash-ab"),
    ]:
        put(
            root,
            f"{ANAT}/{participant}/ses-V02/anat/{participant}_run-1_seg-MyersLabonte_stat-mean_desc-thickness_morph.tsv",
            header + "\nSECRET_VALUE\tSECRET_VALUE\n",
        )
    put(
        root,
        "demographics_sed/sed_basic_demographics.tsv",
        "\ufeffparticipant_id\tsex\tga\r\nsub-private001\tSECRET_VALUE\tSECRET_VALUE\r\n",
    )
    result = run_extractor(root)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    schemas = report["files"]["parcel_thickness"]["schemas"]
    assert len(schemas) == 2
    assert {tuple(x["columns"]) for x in schemas} == {
        ("R_A_2_hash-ab", "L_A_1_hash-cd"),
        ("L_A_1_hash-cd", "R_A_2_hash-ab"),
    }
    assert report["files"]["demographics"]["schemas"][0]["columns"] == [
        "participant_id",
        "sex",
        "ga",
    ]
    assert report["files"]["medications"]["status"] == "missing"
    for secret in ["SECRET_VALUE", "sub-private", str(root)]:
        assert secret not in result.stdout + result.stderr
    assert result.stdout == run_extractor(root).stdout


def test_extracts_matrix_row_labels_but_never_connectivity_values(release_headers):
    put(
        release_headers,
        f"{FUNC}/sub-private/ses-V02/func/sub-private_task-rest_space-fsLR_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv",
        "\tL_A_1_hash-a\tNone\nNone\tSECRET_VALUE\tSECRET_VALUE\nL_A_1_hash-a\tSECRET_VALUE\tSECRET_VALUE\n",
    )
    result = run_extractor(release_headers)
    assert result.returncode == 0, result.stderr
    schema = json.loads(result.stdout)["files"]["connectivity"]["schemas"][0]
    assert schema["columns"] == ["", "L_A_1_hash-a", "None"]
    assert schema["row_labels"] == ["None", "L_A_1_hash-a"]
    assert "SECRET_VALUE" not in result.stdout


@pytest.mark.parametrize("location", ["unset", "nonexistent", "file"])
def test_invalid_release_fails_without_partial_report(tmp_path, location):
    root = None if location == "unset" else tmp_path / "release"
    if location == "file":
        root.write_text("not a directory")
    result = run_extractor(root)
    assert result.returncode != 0
    assert not result.stdout
    assert "HBCD_RELEASE" in result.stderr


@pytest.mark.parametrize(
    "content", ["", "participant_id\tsub-private001\n", "\tL_A_1\nsub-private001\t.5\n"]
)
def test_empty_or_identifier_bearing_matrix_is_not_emitted(release_headers, content):
    put(
        release_headers,
        f"{FUNC}/sub-private/ses-V02/func/sub-private_task-rest_space-fsLR_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv",
        content,
    )
    result = run_extractor(release_headers)
    assert result.returncode != 0
    assert not result.stdout
    assert "sub-private" not in result.stderr
    assert "connectivity" in result.stderr
