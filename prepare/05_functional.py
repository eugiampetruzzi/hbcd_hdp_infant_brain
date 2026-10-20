"""Motion summaries and Myers-Labonte network metrics from the concatenated XCP-D resting-state runs."""

import re

import numpy as np
import pandas as pd

from config import DATA, XCPD, release

TR = 1.725
COVERAGE = 0.80


def read_matrix(path, order=None):
    m = pd.read_csv(path, sep="\t", index_col=0, na_values="n/a")
    m.columns = [re.sub(r"_hash-.*$", "", c) for c in m.columns]
    m.index = [re.sub(r"_hash-.*$", "", str(c)) for c in m.index]
    # Rows and columns must name the same parcels; every matrix is put in one parcel order.
    assert m.columns.is_unique and set(m.index) == set(m.columns), path
    m = m.loc[m.columns, m.columns] if order is None else m.loc[order, order]
    labels = list(m.columns)
    m = m.to_numpy(float)
    np.fill_diagonal(m, np.nan)
    return m, labels


def framewise_displacement(path):
    table = pd.read_csv(path, sep="\t", na_values="n/a")
    return table[
        next(c for c in table if c.startswith("framewise_displacement"))
    ].to_numpy(float)


def motion(path):
    prefix = path.parent / (path.name.split("_task-rest_")[0] + "_task-rest")
    fd = framewise_displacement(f"{prefix}_motion.tsv")
    retained = framewise_displacement(f"{prefix}_outliers.tsv") != 1
    return {
        "mean_fd_post_censoring": np.nanmean(fd[retained]),
        "retained_duration_minutes": retained.sum() * TR / 60,
    }


# Signed weighted modularity of the fixed system partition (positive minus scaled negative
# term; Rubinov & Sporns, 2011) and signed participation coefficients over the same systems.
def network_metrics(z, systems):
    modules = (systems[:, None] == np.unique(systems)).astype(float)
    same = modules @ modules.T
    q, p, s = [], [], []
    for w in (np.clip(z, 0, None), np.clip(-z, 0, None)):
        k = w.sum(1)
        s.append(k.sum())
        # A sign with no edges contributes 0; nodes without edges of a sign have participation 0
        # (Brain Connectivity Toolbox convention).
        with np.errstate(invalid="ignore", divide="ignore"):
            q.append(np.sum((w - np.outer(k, k) / k.sum()) * same) if k.sum() else 0.0)
            p.append(np.where(k > 0, 1 - ((w @ modules) ** 2).sum(1) / k**2, 0.0))
    return {
        "parcel_count": len(z),
        "system_count": modules.shape[1],
        "modularity": (q[0] / s[0] if s[0] else 0.0)
        - (q[1] / (s[0] + s[1]) if s[1] else 0.0),
        "mean_participation_coefficient": np.mean((p[0] + p[1]) / 2),
        "mean_positive_participation": p[0].mean(),
        "mean_negative_participation": p[1].mean(),
        "mean_network_connectivity": z[np.triu_indices(len(z), 1)].mean(),
    }


def main(expected_count=1442, expected_parcels=276):
    """Prepare a release, validating scan and retained parcel counts."""
    # Run-free files hold the runs concatenated within the session.
    pattern = "sub-*/ses-V02/func/*task-rest_space-fsLR_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv"
    paths = sorted(
        p
        for p in release(XCPD).glob(pattern)
        if "_run-" not in p.name and "_dir-" not in p.name
    )
    assert len(paths) == expected_count, len(paths)

    order = read_matrix(paths[0])[1]
    present = 0
    for path in paths:
        m, labels = read_matrix(path, order)
        present = present + np.isfinite(m).any(0)
    # Myers-Labonte labels are <hemisphere>_<system>_<index>; parcel 1 has neither prefix nor index.
    systems = np.array(
        [re.sub(r"_\d+$", "", re.sub(r"^[LR]_", "", label)) for label in labels]
    )
    keep = (present / len(paths) >= COVERAGE) & (systems != "None")
    assert keep.sum() == expected_parcels, keep.sum()

    qc, metrics = [], []
    for path in paths:
        ids = {"participant_id": path.name.split("_")[0], "session_id": "ses-V02"}
        qc.append({**ids, **motion(path)})
        z = np.arctanh(read_matrix(path, order)[0][np.ix_(keep, keep)])
        np.fill_diagonal(z, 0)
        # Only complete matrices over the retained parcels are summarized.
        if np.isfinite(z).all():
            metrics.append(
                {**ids, **network_metrics(z, systems[keep]), "matrix_path": str(path)}
            )

    for rows, rel in (
        (qc, "covariates/v02_functional_qc.tsv"),
        (metrics, "outcomes/functional/v02_functional_network_metrics.tsv"),
    ):
        (DATA / rel).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(DATA / rel, sep="\t", index=False)


if __name__ == "__main__":
    main()
