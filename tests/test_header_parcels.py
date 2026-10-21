"""Release parcel names identify networks and anatomy independently of table order."""

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from nibabel.cifti2 import cifti2_axes as axes

from header_fixtures import release_columns
from synthetic import brain_axis, functional_files, structural_files, write_tsv


def test_full_parcel_networks_exclude_unassigned_and_survive_axis_permutations(modules):
    m = modules("prepare/05_functional.py")
    columns = release_columns("parcel_thickness")
    names = [c.split("_hash-")[0] for c in columns]
    active = [names.index(n) for n in ["PostDefault", "L_PostDefault_1", "L_MedialVisual_1", "L_MedialVisual_2"]]
    matrix = np.eye(len(names))
    matrix[np.ix_(active, active)] = np.array([[1, .4, .2, -.2], [.4, 1, -.2, .2], [.2, -.2, 1, .4], [-.2, .2, .4, 1]])
    unassigned = [i for i, n in enumerate(names) if "_None_" in n]
    matrix[np.ix_(unassigned, unassigned)] = .9
    np.fill_diagonal(matrix, 1)
    frame = pd.DataFrame(matrix, index=columns, columns=columns)
    path = functional_files(m.release(m.XCPD))
    results = []
    for permutation in [False, True]:
        data = frame.sample(frac=1, random_state=7).sample(frac=1, axis=1, random_state=11) if permutation else frame
        write_tsv(path, data.rename_axis("Node"), index=True)
        m.main(expected_count=1, expected_parcels=279)
        results.append(pd.read_csv(m.DATA / "outcomes/functional/v02_functional_network_metrics.tsv", sep="\t"))
    pd.testing.assert_frame_equal(*results)
    row = results[0].iloc[0]
    assert row.parcel_count == 279
    assert row.system_count == 23
    assert row.modularity == pytest.approx(7 / 24)
    assert row.mean_positive_participation == pytest.approx(16 / (9 * 279))
    assert row.mean_negative_participation == 0
    assert row.mean_participation_coefficient == pytest.approx(8 / (9 * 279))
    assert row.mean_network_connectivity == pytest.approx(.8 / (279 * 278 / 2))


def test_structural_full_headers_preserve_named_values_after_column_permutation(modules):
    m = modules("prepare/04_structural.py")
    path = structural_files(m.release(m.XCPD_ANAT))
    for metric in ["thickness", "sulc"]:
        columns = release_columns(f"parcel_{metric}")
        frame = pd.DataFrame([np.arange(len(columns)) + .25], columns=columns)
        target = next(path.parent.glob(f"*desc-{metric}_morph.tsv"))
        write_tsv(target, frame.iloc[:, ::-1])
    m.main(expected_count=1)
    for metric in ["thickness", "sulc"]:
        row = pd.read_csv(m.DATA / f"outcomes/structural/v02_myerslabonte_{metric}.tsv", sep="\t").iloc[0]
        for i, column in enumerate(release_columns(f"parcel_{metric}")):
            assert row[column.split("_hash-")[0]] == i + .25


def test_same_hemisphere_parcel_swaps_cannot_relocate_significant_thickening(modules):
    m = modules("python/figure4_map.py")
    names = [c.split("_hash-")[0] for c in release_columns("parcel_thickness")]
    positive, negative = "L_PostDefault_1", "L_PostDefault_2"
    m.ATLAS.mkdir(parents=True)
    m.TABLES.mkdir(parents=True)
    # Label numbers deliberately differ from TSV positions and are noncontiguous.
    labels = {0: ("background", (0, 0, 0, 0)), 7: (positive, (1, 0, 0, 1)),
              11: (negative, (0, 1, 0, 1))}
    header = nib.Cifti2Header.from_axes((axes.LabelAxis(["synthetic"], [labels]), brain_axis(32492)))
    nib.save(nib.Cifti2Image(np.array([[7, 11, 0, 7, 11, 0, 0, 0, 0, 0]], np.int32), header), m.ATLAS / "parcels.dlabel.nii")
    pd.DataFrame(dict(metric=["Parcel cortical thickness"] * 2, feature=[positive, negative],
                      beta=[.3, -.2], q=[.01, .01])).to_csv(m.TABLES / "regional.csv", index=False)
    for ordered in [names, names[::-1]]:
        write_tsv(m.DATA / "outcomes/structural/v02_myerslabonte_thickness.tsv",
                  pd.DataFrame(columns=["participant_id", "session_id", "source_run", *ordered]))
        maps = m.surface_maps(*m.parcel_values())
        np.testing.assert_allclose(maps["L"][:5], [.3, 0, 0, .3, 0])
        assert not maps["R"].any()
