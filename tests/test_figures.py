import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from nibabel.cifti2 import cifti2_axes as axes

from synthetic import brain_axis, write_tsv


@pytest.fixture
def figure_data(modules):
    m = modules("python/figure4_map.py")
    m.ATLAS.mkdir(parents=True)
    m.TABLES.mkdir(parents=True)
    m.FIGURES.mkdir(parents=True)
    write_tsv(m.DATA / "outcomes/structural/v02_myerslabonte_thickness.tsv",
              {"participant_id": [], "session_id": [], "source_run": [],
               "L_A_1": [], "R_A_2": [], "L_B_1": []})
    pd.DataFrame({"metric": ["Parcel cortical thickness"] * 3 + ["Parcel sulcal depth"],
                  "feature": ["L_A_1", "R_A_2", "L_B_1", "L_A_1"],
                  "q": [.01, .05, .2, .001], "beta": [.2, .3, .4, 9]}).to_csv(
        m.TABLES / "regional.csv", index=False)
    labels = {i: (f"parcel{i}", (i / 4, 0., 0., 1.)) for i in range(4)}
    header = nib.Cifti2Header.from_axes((axes.LabelAxis(["parcels"], [labels]), brain_axis(32492)))
    nib.save(nib.Cifti2Image(np.array([[1, 2, 0, 3, 1, 2, 1, 3, 0, 2]], np.int32), header),
             m.ATLAS / "parcels.dlabel.nii")
    for hemi in ["L", "R"]:
        coords = np.zeros((32492, 3), np.float32)
        coords[:5] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]
        faces = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 4], [2, 3, 4]], np.int32)
        image = nib.gifti.GiftiImage(darrays=[
            nib.gifti.GiftiDataArray(coords, intent="NIFTI_INTENT_POINTSET"),
            nib.gifti.GiftiDataArray(faces, intent="NIFTI_INTENT_TRIANGLE")])
        nib.save(image, m.ATLAS / f"fsLR_{hemi}_inflated.surf.gii")
    bars = np.ones((100, 120, 3))
    bars[10:90, 10:100] = [.1, .3, .7]
    plt.imsave(m.FIGURES / "figure4b_networks.png", bars)
    yield m
    plt.close("all")


def test_parcel_values_use_atlas_order_and_strict_fdr(figure_data):
    np.testing.assert_array_equal(figure_data.parcel_values(), [0, .2, 0, 0])


def test_surface_maps_place_values_at_correct_vertices(figure_data):
    maps = figure_data.surface_maps(np.array([0, .1, .2, .3]))
    assert set(maps) == {"L", "R"}
    assert maps["L"].shape == (32492,)
    np.testing.assert_allclose(maps["L"][:5], [.1, .2, 0, .3, .1])
    np.testing.assert_allclose(maps["R"][:5], [.2, .1, .3, 0, .2])
    assert not maps["L"][5:].any()


def test_crop_preserves_colored_region_and_clips_padding(modules):
    m = modules("python/figure4_map.py")
    image = np.ones((12, 14, 4))
    image[3:8, 4:9, :3] = .2
    result = m.crop(image, pad=2)
    assert result.shape == (8, 8, 4)
    assert (result[..., :3].min(axis=2) < .97).sum() == 25
    np.testing.assert_array_equal(m.crop(image, pad=30), image)


def test_figure_main_renders_real_surfaces_and_composes_outputs(figure_data):
    m = figure_data
    m.main()
    for filename in ["figure4a_map.png", "figure4_cortical_thickness.png"]:
        image = plt.imread(m.FIGURES / filename)
        assert image.shape[0] > 100
        assert image.shape[1] > 100
        assert np.any(image[..., :3] < .5)
    assert (m.FIGURES / "figure4_cortical_thickness.pdf").read_bytes().startswith(b"%PDF")
