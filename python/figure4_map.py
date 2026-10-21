"""Figure 4: parcels with thicker cortex in HDP infants (a) next to the network bars (b)."""

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from nilearn import plotting

from config import ATLAS, DATA, FIGURES, TABLES

plt.rcParams["font.family"] = "Arial"
VMAX = 0.4
CMAP = LinearSegmentedColormap.from_list("thick", ["#C9DCF0", "#2E6DB4", "#0B2E5C"])


def parcel_values():
    # Values stay paired with parcel names; the CIFTI label table supplies anatomy.
    thickness = pd.read_csv(
        DATA / "outcomes/structural/v02_myerslabonte_thickness.tsv", sep="\t", nrows=0
    )
    parcels = [
        c
        for c in thickness.columns
        if c not in ("participant_id", "session_id", "source_run")
    ]
    regional = pd.read_csv(TABLES / "regional.csv")
    sig = regional[
        (regional.metric == "Parcel cortical thickness")
        & (regional.q < 0.05)
        & (regional.beta > 0)
    ]
    beta = dict(zip(sig.feature, sig.beta))
    return np.array([0.0] + [beta.get(p, 0.0) for p in parcels]), parcels


def surface_maps(values, parcels):
    img = nib.load(ATLAS / "parcels.dlabel.nii")
    labels = np.asarray(img.get_fdata()).ravel().astype(int)
    atlas_names = img.header.get_axis(0).label[0]
    by_name = dict(zip(parcels, values[1:], strict=True))
    lookup = np.zeros(labels.max() + 1)
    for label in np.unique(labels):
        if label > 0:
            parcel = atlas_names[label][0].split("_hash-")[0]
            # Unverified positional matching can silently move an effect within a hemisphere.
            if parcel not in by_name:
                raise ValueError(f"Atlas label {label} names unknown parcel {parcel!r}")
            lookup[label] = by_name[parcel]
    maps = {}
    for name, slc, model in img.header.get_axis(1).iter_structures():
        hemi = "L" if name.endswith("LEFT") else "R"
        maps[hemi] = np.zeros(32492)
        maps[hemi][model.vertex] = lookup[labels[slc]]
    return maps


def render_view(maps, hemi, view):
    fig = plt.figure(figsize=(4, 3.2))
    ax = fig.add_axes([0, 0, 1, 1], projection="3d")
    surf = nib.load(ATLAS / f"fsLR_{hemi}_inflated.surf.gii")
    plotting.plot_surf_stat_map(
        [surf.darrays[0].data, surf.darrays[1].data],
        maps[hemi],
        hemi="left" if hemi == "L" else "right",
        view=view,
        bg_map=None,
        cmap=CMAP,
        vmin=0,
        vmax=VMAX,
        threshold=1e-6,
        symmetric_cbar=False,
        colorbar=False,
        axes=ax,
        figure=fig,
        engine="matplotlib",
    )
    buf = io.BytesIO()
    fig.savefig(buf, dpi=250, facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return crop(mpimg.imread(buf), pad=4)[..., :3]


def tile(rows, gap=12):
    """Place cropped views edge to edge with a small white gap, centering each in its cell."""
    h = max(img.shape[0] for row in rows for img in row)
    w = max(img.shape[1] for row in rows for img in row)
    out = np.ones(
        (
            len(rows) * h + (len(rows) - 1) * gap,
            len(rows[0]) * w + (len(rows[0]) - 1) * gap,
            3,
        )
    )
    for r, row in enumerate(rows):
        for c, img in enumerate(row):
            y = r * (h + gap) + (h - img.shape[0]) // 2
            x = c * (w + gap) + (w - img.shape[1]) // 2
            out[y : y + img.shape[0], x : x + img.shape[1]] = img
    return out


def draw_map(maps, out):
    views = [[("L", "lateral"), ("R", "lateral")], [("L", "medial"), ("R", "medial")]]
    brains = tile(
        [[render_view(maps, hemi, view) for hemi, view in row] for row in views]
    )
    width = 7.2
    brain_h = width * brains.shape[0] / brains.shape[1]
    bar_h = 1.0
    fig = plt.figure(figsize=(width, brain_h + bar_h))
    total = brain_h + bar_h
    ax = fig.add_axes([0, bar_h / total, 1, brain_h / total])
    ax.imshow(brains)
    ax.axis("off")
    cax = fig.add_axes([0.3, 0.62 * bar_h / total, 0.4, 0.2 * bar_h / total])
    bar = fig.colorbar(
        plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, VMAX)),
        cax=cax,
        orientation="horizontal",
        ticks=[0, 0.2, 0.4],
    )
    bar.ax.tick_params(labelsize=13)
    bar.set_label("Thicker cortex in HDP infants (SD)", fontsize=13)
    fig.savefig(out, dpi=250, facecolor="white")
    plt.close(fig)


def crop(img, pad=10):
    rows, cols = np.where(img[..., :3].min(axis=2) < 0.97)
    return img[
        max(rows.min() - pad, 0) : rows.max() + pad,
        max(cols.min() - pad, 0) : cols.max() + pad,
    ]


def compose(map_png, bars_png):
    a, b = crop(mpimg.imread(map_png)), crop(mpimg.imread(bars_png))
    width_in = 183 / 25.4
    wa, wb = 0.52, 0.46
    height = max(wa * a.shape[0] / a.shape[1], wb * b.shape[0] / b.shape[1]) + 0.02
    fig = plt.figure(figsize=(width_in, width_in * height))
    for img, left, w in [(a, 0, wa), (b, 1 - wb, wb)]:
        h = w * img.shape[0] / img.shape[1]
        ax = fig.add_axes([left, (height - h) / 2 / height, w, h / height])
        ax.imshow(img)
        ax.axis("off")
    fig.text(0.005, 0.985, "a", fontsize=12, fontweight="bold", va="top")
    fig.text(1 - wb - 0.01, 0.985, "b", fontsize=12, fontweight="bold", va="top")
    for ext in ("png", "pdf"):
        fig.savefig(
            FIGURES / f"figure4_cortical_thickness.{ext}", dpi=300, facecolor="white"
        )


def main():
    map_png = FIGURES / "figure4a_map.png"
    draw_map(surface_maps(*parcel_values()), map_png)
    compose(map_png, FIGURES / "figure4b_networks.png")


if __name__ == "__main__":
    main()
