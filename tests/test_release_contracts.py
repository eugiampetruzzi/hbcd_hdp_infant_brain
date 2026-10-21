"""Release-path and numerical contracts independent of production constants."""

import numpy as np
import pandas as pd
import pytest

from synthetic import functional_files, scalar, structural_files, write_tsv

ANATOMICAL = "imaging/structure_function/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"
FUNCTIONAL = "imaging/func_xcpd/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"


def test_structural_reads_anatomical_release_not_functional_distractor(modules):
    m = modules("prepare/04_structural.py")
    structural_files(m.release(ANATOMICAL))
    distractor = structural_files(m.release(FUNCTIONAL))
    scalar(distractor, np.full(10, 99))
    m.main(expected_count=1)
    result = pd.read_csv(
        m.DATA / "outcomes/structural/v02_global_cortical_morphometry.tsv", sep="\t"
    )
    assert result.participant_id.tolist() == ["sub-test"]
    assert result.global_area_weighted_thickness.tolist() == pytest.approx([2.8])


@pytest.fixture
def pearson_release(modules):
    m = modules("prepare/05_functional.py")
    # Positive strengths are .6, negative strengths .2 at each node.
    # Within-system positive weight is 1.6; the null expectation is 1.2.
    # Q+ = 1/6, Q- contribution = -1/8, hence signed Q = 7/24.
    # Positive participation is 1 - (2/3)^2 - (1/3)^2 = 4/9.
    matrix = np.array(
        [
            [1, 0.4, 0.2, -0.2],
            [0.4, 1, -0.2, 0.2],
            [0.2, -0.2, 1, 0.4],
            [-0.2, 0.2, 0.4, 1],
        ]
    )
    for participant, reverse in [("sub-ordered", False), ("sub-permuted", True)]:
        path = functional_files(m.release(FUNCTIONAL), participant)
        frame = pd.read_csv(path, sep="\t", index_col=0)
        frame.iloc[:, :] = matrix
        if reverse:
            frame = frame.iloc[[2, 0, 3, 1], [3, 1, 0, 2]]
        write_tsv(path, frame, index=True)
    functional_files(m.release(ANATOMICAL), "sub-wrong-folder")
    return m


def test_pearson_pipeline_metrics_and_subject_alignment(pearson_release):
    m = pearson_release
    m.main(expected_count=2, expected_parcels=4)
    result = pd.read_csv(
        m.DATA / "outcomes/functional/v02_functional_network_metrics.tsv", sep="\t"
    )
    assert result.participant_id.tolist() == ["sub-ordered", "sub-permuted"]
    expected = {
        "parcel_count": 4,
        "system_count": 2,
        "modularity": 7 / 24,
        "mean_positive_participation": 4 / 9,
        "mean_negative_participation": 0,
        "mean_participation_coefficient": 2 / 9,
        "mean_network_connectivity": 2 / 15,
    }
    for column, value in expected.items():
        np.testing.assert_allclose(result[column], [value, value], atol=1e-12)


@pytest.mark.parametrize(
    "kind,modularity,positive,negative,connectivity",
    [
        # Q+ = 1/2; negative within-system residual = -4, total strength = 12.
        ("perfect", 5 / 6, 0, 0, -1 / 3),
        ("positive", -1 / 6, 4 / 9, 0, 0.2),
        ("negative", 1 / 6, 0, 4 / 9, -0.2),
        ("zero", 0, 0, 0, 0),
    ],
)
def test_pearson_pipeline_retains_boundary_and_single_sign_graphs(
    modules, kind, modularity, positive, negative, connectivity
):
    m = modules("prepare/05_functional.py")
    path = functional_files(m.release(FUNCTIONAL))
    frame = pd.read_csv(path, sep="\t", index_col=0)
    if kind == "perfect":
        # A valid rank-one correlation matrix, including off-diagonal +/-1.
        matrix = np.outer([1, 1, -1, -1], [1, 1, -1, -1])
    else:
        matrix = np.full((4, 4), {"positive": 0.2, "negative": -0.2, "zero": 0}[kind])
        np.fill_diagonal(matrix, 1)
    frame.iloc[:, :] = matrix
    write_tsv(path, frame, index=True)
    with np.errstate(all="raise"):
        m.main(expected_count=1, expected_parcels=4)
    result = pd.read_csv(
        m.DATA / "outcomes/functional/v02_functional_network_metrics.tsv", sep="\t"
    )
    assert result.participant_id.tolist() == ["sub-test"]
    assert result.modularity.iloc[0] == pytest.approx(modularity)
    assert result.mean_positive_participation.iloc[0] == pytest.approx(positive)
    assert result.mean_negative_participation.iloc[0] == pytest.approx(negative)
    assert result.mean_participation_coefficient.iloc[0] == pytest.approx(
        (positive + negative) / 2
    )
    assert result.mean_network_connectivity.iloc[0] == pytest.approx(connectivity)
