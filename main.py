import hfp
import matplotlib.pyplot as plt
import trimesh
import numpy as np
import random
import scipy

def wrapper_hfp(vs, fs, num_clusters):
    """
    Wraps HFP C++ extension.
    """
    out = hfp.run_hfp(vs, fs, num_clusters)
    idx = np.array(out)

    return np.flip(idx).astype(int)


def load_vf(file_name):
    """
    Loads file with vedo.
    """
    mesh = trimesh.load_mesh(file_name)
    vs = np.array(mesh.vertices)
    fs = np.array(mesh.faces)
    return vs, fs

def map2range(arr):
    arr_uni, _, counts = np.unique(arr, return_counts=True, return_index=True)
    arr_range = np.arange(len(arr_uni))
    arr_new = np.repeat(arr_range, counts)
    
    idx = np.argsort(arr)
    idx_rev = np.argsort(idx)
    return arr_new[idx_rev]

def segments_split(vs, fs, idx):
    segments = {}
    for segment in np.unique(idx):
        segment_idx = np.isin(idx, segment).astype(int)
        segment_idx = np.where(segment_idx == 1)[0]
        fs_part = fs[segment_idx]
        vs_split = vs[np.unique(fs_part)]
        fs_split = map2range(fs_part.reshape(-1)).reshape(-1, 3)
        mesh = trimesh.Trimesh(vs_split, fs_split)
        segments[segment] = (mesh)
    return segments

def visualize_select(meshes):
    colors = np.array([
        plt.get_cmap("hsv")(int(c*256/len(meshes))) for c in range(len(meshes))])

    random.shuffle(colors)
    for mesh, color in zip(meshes, colors):
        mesh.visual.face_colors = np.repeat(color.reshape(1,-1), len(mesh.faces), axis=0)
    
    trimesh.Scene(meshes).show()

vs, fs = load_vf("data/cup.obj")
mesh = trimesh.Trimesh(vs, fs)

idx = wrapper_hfp(vs, fs, 2)

segments = segments_split(vs, fs, idx)

visualize_select(segments.values())