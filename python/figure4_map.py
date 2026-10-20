"""Figure 4: parcels with thicker cortex in HDP infants (a) next to the network bars (b)."""

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
    # Atlas label i is the i-th parcel column of the HBCD Myers-Labonte table.
    thickness = pd.read_csv(DATA / "outcomes/structural/v02_myerslabonte_thickness.tsv", sep="\t", nrows=0)
    parcels = [c for c in thickness.columns if c not in ("participant_id", "session_id", "source_run")]
    regional = pd.read_csv(TABLES / "regional.csv")
    sig = regional[(regional.metric == "Parcel cortical thickness") & (regional.q < .05) & (regional.beta > 0)]
    beta = dict(zip(sig.feature, sig.beta))
    return np.array([0.0] + [beta.get(p, 0.0) for p in parcels]), parcels


def surface_maps(values, parcels):
    img = nib.load(ATLAS / "parcels.dlabel.nii")
    labels = np.asarray(img.get_fdata()).ravel().astype(int)
    maps = {}
    for name, slc, model in img.header.get_axis(1).iter_structures():
        hemi = "L" if name.endswith("LEFT") else "R"
        # Column order must match atlas label order: every parcel's hemisphere prefix must
        # agree with the hemisphere its atlas label is drawn on (parcel 1 has no prefix, left).
        for label in np.unique(labels[slc]):
            if label > 0:
                name_hemi = "R" if parcels[label - 1].startswith("R_") else "L"
                assert name_hemi == hemi, (label, parcels[label - 1], hemi)
        maps[hemi] = np.zeros(32492)
        maps[hemi][model.vertex] = values[labels[slc]]
    return maps


def draw_map(maps, out):
    views = [("L", "lateral"), ("R", "lateral"), ("L", "medial"), ("R", "medial")]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), subplot_kw={"projection": "3d"})
    for ax, (hemi, view) in zip(axes.flat, views):
        surf = nib.load(ATLAS / f"fsLR_{hemi}_inflated.surf.gii")
        plotting.plot_surf_stat_map([surf.darrays[0].data, surf.darrays[1].data], maps[hemi],
                                    hemi="left" if hemi == "L" else "right", view=view, bg_map=None,
                                    cmap=CMAP, vmin=0, vmax=VMAX, threshold=1e-6, symmetric_cbar=False,
                                    colorbar=False, axes=ax, figure=fig, engine="matplotlib")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0.17, wspace=-0.28, hspace=-0.15)
    cax = fig.add_axes([0.33, 0.11, 0.34, 0.022])
    bar = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, VMAX)), cax=cax,
                       orientation="horizontal", ticks=[0, 0.2, 0.4])
    bar.ax.tick_params(labelsize=10)
    bar.set_label("Thicker cortex in HDP infants (SD)", fontsize=10)
    fig.savefig(out, dpi=250, facecolor="white")
    plt.close(fig)


def crop(img, pad=10):
    rows, cols = np.where(img[..., :3].min(axis=2) < 0.97)
    return img[max(rows.min() - pad, 0):rows.max() + pad, max(cols.min() - pad, 0):cols.max() + pad]


def compose(map_png, bars_png):
    a, b = crop(mpimg.imread(map_png)), crop(mpimg.imread(bars_png))
    width_in = 183 / 25.4
    wa, wb = 0.60, 0.38
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
        fig.savefig(FIGURES / f"figure4_cortical_thickness.{ext}", dpi=300, facecolor="white")


def main():
    map_png = FIGURES / "figure4a_map.png"
    draw_map(surface_maps(*parcel_values()), map_png)
    compose(map_png, FIGURES / "figure4b_networks.png")


if __name__ == "__main__":
    main()
