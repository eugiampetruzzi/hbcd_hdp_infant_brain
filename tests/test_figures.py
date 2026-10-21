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
    write_tsv(
        m.DATA / "outcomes/structural/v02_myerslabonte_thickness.tsv",
        {
            "participant_id": [],
            "session_id": [],
            "source_run": [],
            "L_A_1": [],
            "R_A_2": [],
            "L_B_1": [],
        },
    )
    pd.DataFrame(
        {
            "metric": ["Parcel cortical thickness"] * 3 + ["Parcel sulcal depth"],
            "feature": ["L_A_1", "R_A_2", "L_B_1", "L_A_1"],
            "q": [0.01, 0.05, 0.2, 0.001],
            "beta": [0.2, 0.3, 0.4, 9],
        }
    ).to_csv(m.TABLES / "regional.csv", index=False)
    labels = {i: (name, (i / 4, 0.0, 0.0, 1.0))
              for i, name in enumerate(["background", "L_A_1", "R_A_2", "L_B_1"])}
    header = nib.Cifti2Header.from_axes(
        (axes.LabelAxis(["parcels"], [labels]), brain_axis(32492))
    )
    nib.save(
        nib.Cifti2Image(np.array([[1, 3, 0, 3, 1, 2, 2, 2, 0, 2]], np.int32), header),
        m.ATLAS / "parcels.dlabel.nii",
    )
    for hemi in ["L", "R"]:
        coords = np.zeros((32492, 3), np.float32)
        coords[:5] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]
        faces = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 4], [2, 3, 4]], np.int32)
        image = nib.gifti.GiftiImage(
            darrays=[
                nib.gifti.GiftiDataArray(coords, intent="NIFTI_INTENT_POINTSET"),
                nib.gifti.GiftiDataArray(faces, intent="NIFTI_INTENT_TRIANGLE"),
            ]
        )
        nib.save(image, m.ATLAS / f"fsLR_{hemi}_inflated.surf.gii")
    bars = np.ones((100, 120, 3))
    bars[10:90, 10:100] = [0.1, 0.3, 0.7]
    plt.imsave(m.FIGURES / "figure4b_networks.png", bars)
    yield m
    plt.close("all")


def test_parcel_values_use_atlas_order_and_strict_fdr(figure_data):
    values, parcels = figure_data.parcel_values()
    np.testing.assert_array_equal(values, [0, 0.2, 0, 0])
    assert parcels == ["L_A_1", "R_A_2", "L_B_1"]


def test_parcel_values_exclude_significant_thinning(figure_data):
    path = figure_data.TABLES / "regional.csv"
    regional = pd.read_csv(path)
    regional.loc[0, "beta"] = -0.2
    regional.to_csv(path, index=False)
    values, _ = figure_data.parcel_values()
    np.testing.assert_array_equal(values, [0, 0, 0, 0])


def test_surface_maps_place_values_at_correct_vertices(figure_data):
    maps = figure_data.surface_maps(
        np.array([0, 0.1, 0.2, 0.3]), ["L_A_1", "R_A_2", "L_B_1"]
    )
    assert set(maps) == {"L", "R"}
    assert maps["L"].shape == (32492,)
    np.testing.assert_allclose(maps["L"][:5], [0.1, 0.3, 0, 0.3, 0.1])
    np.testing.assert_allclose(maps["R"][:5], [0.2, 0.2, 0.2, 0, 0.2])
    assert not maps["L"][5:].any()


def test_surface_maps_follow_names_when_hemisphere_order_changes(figure_data):
    maps = figure_data.surface_maps(
        np.array([0, 0.1, 0.2, 0.3]), ["R_A_2", "L_A_1", "L_B_1"]
    )
    np.testing.assert_allclose(maps["L"][:5], [.2, .3, 0, .3, .2])
    np.testing.assert_allclose(maps["R"][:5], [.1, .1, .1, 0, .1])


def test_crop_preserves_colored_region_and_clips_padding(modules):
    m = modules("python/figure4_map.py")
    image = np.ones((12, 14, 4))
    image[3:8, 4:9, :3] = 0.2
    result = m.crop(image, pad=2)
    assert result.shape == (8, 8, 4)
    assert (result[..., :3].min(axis=2) < 0.97).sum() == 25
    np.testing.assert_array_equal(m.crop(image, pad=30), image)


def test_figure_main_renders_real_surfaces_and_composes_outputs(figure_data):
    m = figure_data
    m.main()
    for filename in ["figure4a_map.png", "figure4_cortical_thickness.png"]:
        image = plt.imread(m.FIGURES / filename)
        assert image.shape[0] > 100
        assert image.shape[1] > 100
        assert np.any(image[..., :3] < 0.5)
    assert (
        (m.FIGURES / "figure4_cortical_thickness.pdf").read_bytes().startswith(b"%PDF")
    )


def test_tile_centers_unequal_images_without_losing_pixels(modules):
    m = modules("python/figure4_map.py")
    red = np.full((2, 4, 3), [1, 0, 0])
    green = np.full((4, 2, 3), [0, 1, 0])
    blue = np.full((1, 1, 3), [0, 0, 1])
    black = np.zeros((4, 4, 3))
    tiled = m.tile([[red, green], [blue, black]], gap=2)
    expected = np.ones((10, 10, 3))
    expected[1:3, 0:4] = red
    expected[0:4, 7:9] = green
    expected[7:8, 1:2] = blue
    expected[6:10, 6:10] = black
    np.testing.assert_array_equal(tiled, expected)


@pytest.mark.parametrize(
    "hemi,view", [("L", "lateral"), ("R", "lateral"), ("L", "medial"), ("R", "medial")]
)
def test_render_view_produces_rgb_and_closes_its_figure(figure_data, hemi, view):
    m = figure_data
    maps = m.surface_maps(np.array([0, 0.1, 0.2, 0.3]), ["L_A_1", "R_A_2", "L_B_1"])
    before = plt.get_fignums()
    rendered = m.render_view(maps, hemi, view)
    assert rendered.ndim == 3 and rendered.shape[2] == 3
    assert min(rendered.shape[:2]) > 10
    assert np.isfinite(rendered).all()
    assert ((rendered >= 0) & (rendered <= 1)).all()
    assert np.any(rendered < 0.5)
    assert plt.get_fignums() == before


def test_draw_map_keeps_view_order_and_colorbar_below_brains(figure_data, monkeypatch):
    m = figure_data
    colors = {
        ("L", "lateral"): [1, 0, 0],
        ("R", "lateral"): [0, 1, 0],
        ("L", "medial"): [0, 0, 1],
        ("R", "medial"): [0, 0, 0],
    }
    # Isolate expensive surface rendering; exercise actual tiling, layout, and PNG output.
    monkeypatch.setattr(
        m,
        "render_view",
        lambda maps, hemi, view: np.full((20, 30, 3), colors[hemi, view], dtype=float),
    )
    path = m.FIGURES / "layout.png"
    m.draw_map({}, path)
    image = plt.imread(path)[..., :3]
    centers = []
    for color in colors.values():
        mask = np.all(np.isclose(image, color, atol=0.01), axis=2)
        # The black tile dominates any black colorbar lettering.
        ys, xs = np.where(mask[: int(0.8 * len(image))])
        assert len(xs) > 100
        centers.append((np.median(xs), np.median(ys)))
    assert centers[0][0] < centers[1][0]
    assert centers[2][0] < centers[3][0]
    assert centers[0][1] < centers[2][1]
    assert centers[1][1] < centers[3][1]
    # Colorbar pixels remain in the region below the tiles.
    bottom = image[int(0.85 * len(image)) :]
    assert np.any(bottom[..., 2] - bottom[..., 0] > 0.15)


def test_figure_renders_when_no_parcel_is_significant(figure_data):
    m = figure_data
    path = m.TABLES / "regional.csv"
    regional = pd.read_csv(path)
    regional["q"] = 1.0
    regional.to_csv(path, index=False)
    values, parcels = m.parcel_values()
    assert not values.any()
    maps = m.surface_maps(values, parcels)
    assert not maps["L"].any() and not maps["R"].any()
    m.main()
    image = plt.imread(m.FIGURES / "figure4a_map.png")
    assert np.isfinite(image).all()
    assert min(image.shape[:2]) > 100
