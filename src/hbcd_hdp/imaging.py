from __future__ import annotations

import numpy as np


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not keep.any():
        return float("nan")
    return float(np.average(values[keep], weights=weights[keep]))


def vertex_areas(coordinates: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(coordinates, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    a, b, c = (coordinates[triangles[:, i]] for i in range(3))
    triangle_area = np.linalg.norm(np.cross(b - a, c - a), axis=1) / 2.0
    areas = np.zeros(len(coordinates), dtype=float)
    for position in range(3):
        np.add.at(areas, triangles[:, position], triangle_area / 3.0)
    return areas


def fixed_community_modularity(matrix: np.ndarray, communities: np.ndarray) -> float:
    weights = np.maximum(np.asarray(matrix, dtype=float), 0.0).copy()
    np.fill_diagonal(weights, 0.0)
    strength = weights.sum(axis=1)
    total = strength.sum()
    if total <= 0:
        return float("nan")
    expected = np.outer(strength, strength) / total
    same = communities[:, None] == communities[None, :]
    return float(((weights - expected) * same).sum() / total)


def participation_coefficient(matrix: np.ndarray, communities: np.ndarray) -> np.ndarray:
    weights = np.asarray(matrix, dtype=float).copy()
    np.fill_diagonal(weights, 0.0)
    strength = weights.sum(axis=1)
    result = np.full(len(weights), np.nan)
    valid = strength != 0
    fractions = np.zeros(len(weights), dtype=float)
    for community in np.unique(communities):
        within = weights[:, communities == community].sum(axis=1)
        fractions += np.square(np.divide(within, strength, out=np.zeros_like(within), where=valid))
    result[valid] = 1.0 - fractions[valid]
    return result


def signed_participation(matrix: np.ndarray, communities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(matrix, dtype=float)
    positive = participation_coefficient(np.maximum(matrix, 0.0), communities)
    negative = participation_coefficient(np.maximum(-matrix, 0.0), communities)
    return positive, negative
