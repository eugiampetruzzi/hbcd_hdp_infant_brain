import numpy as np

from hbcd_hdp.imaging import (
    fixed_community_modularity,
    participation_coefficient,
    vertex_areas,
    weighted_mean,
)


def test_weighted_mean_ignores_invalid_and_zero_weights():
    assert weighted_mean(np.array([1, 3, np.nan]), np.array([1, 1, 10])) == 2
    assert np.isnan(weighted_mean(np.array([1]), np.array([0])))


def test_vertex_area_is_shared_across_triangle_vertices():
    coordinates = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    areas = vertex_areas(coordinates, np.array([[0, 1, 2]]))
    np.testing.assert_allclose(areas, np.repeat(1 / 6, 3))


def test_network_metrics_on_two_disconnected_communities():
    matrix = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=float)
    communities = np.array([0, 0, 1, 1])
    assert np.isclose(fixed_community_modularity(matrix, communities), 0.5)
    np.testing.assert_allclose(participation_coefficient(matrix, communities), 0)
