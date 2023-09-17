import numpy as np
import trimesh

from pyhfp import hfp


def segment(vertices, faces, num_segments):
    """
    Wraps the HFP C++ extension into Numpy and segments the data.

    Args:
        vertices (numpy.ndarray): Input array representing vertices.
        faces (numpy.ndarray): Input array representing faces.
        num_segments (int): Number of desired segments.

    Returns:
        numpy.ndarray: An array of segment indices per vertex.
    """
    if not isinstance(vertices, np.ndarray) or not isinstance(faces, np.ndarray):
        raise ValueError("Both 'vertices' and 'faces' must be numpy arrays.")
    if not isinstance(num_segments, int) or num_segments <= 0:
        raise ValueError("'num_segments' must be a positive integer.")
    if num_segments > vertices.shape[0]:
        raise ValueError("'num_segments' cannot be larger than the amount of 'vertices'.")
    if not is_watertight(vertices, faces):
        raise ValueError("Mesh needs to be watertight.")
        
    try:
        out = hfp.run_hfp(vertices, faces, num_segments)
        idx = np.array(out)
    except Exception as e:
        raise RuntimeError(f"Error during segmentation: {str(e)}")

    return np.flip(idx).astype(int)


def is_watertight(vertices, faces):
    """
    Check if mesh is waatertight.

    Args:
        vertices (numpy.ndarray): Input array representing vertices.
        faces (numpy.ndarray): Input array representing faces.
        num_segments (int): Number of desired segments.

    Returns:
        bool: True if watertight.
    """
    mesh = trimesh.Trimesh(vertices, faces)
    return mesh.is_watertight