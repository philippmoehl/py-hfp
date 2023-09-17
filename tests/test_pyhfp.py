import numpy as np
import pytest

from pyhfp import segment


def test_segment_valid_input():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    faces = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    num_segments = 2
    indices = segment(vertices, faces, num_segments)
    assert isinstance(indices, np.ndarray)
    assert len(indices) == len(vertices)

def test_segment_invalid_input():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    faces = np.array([[0, 1, 2]])
    with pytest.raises(ValueError):
        segment(vertices, faces, 5)  # num_segments > number of vertices

    with pytest.raises(ValueError):
        segment([0, 0, 0], [1, 1, 1], 1)  # not a numpy array

    with pytest.raises(ValueError):
        segment(vertices, faces, 0)  # non-positive num_segments
