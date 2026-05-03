"""Coplanar-adjacent face clustering — the only piece of the original
cluster_sig.py the preauthorizer needs. The OBJ reader and signature
helpers were dropped; this module is OBJ-free.
"""
import numpy as np


def face_clusters(verts, faces, eps=1e-5):
    """Return list of clusters; each cluster is a list of face indices.
    Two adjacent faces are in the same cluster iff they're coplanar to eps."""
    nf = len(faces)
    norms = np.zeros((nf, 3))
    offs = np.zeros(nf)
    for i, f in enumerate(faces):
        v1 = verts[f[1]] - verts[f[0]]
        v2 = verts[f[2]] - verts[f[0]]
        n = np.cross(v1, v2)
        nm = np.linalg.norm(n)
        if nm > 1e-12:
            n /= nm
        norms[i] = n
        offs[i] = -np.dot(n, verts[f[0]])

    edge_to_faces = {}
    for i, f in enumerate(faces):
        for u, v in [(f[0], f[1]), (f[1], f[2]), (f[0], f[2])]:
            e = (min(u, v), max(u, v))
            edge_to_faces.setdefault(e, []).append(i)

    parent = list(range(nf))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb

    for e, fs in edge_to_faces.items():
        if len(fs) != 2:
            continue
        a, b = fs
        opp_b = [v for v in faces[b] if v not in (e[0], e[1])][0]
        d = abs(norms[a] @ verts[opp_b] + offs[a])
        if d < eps:
            union(a, b)

    clusters = {}
    for i in range(nf):
        clusters.setdefault(find(i), []).append(i)
    return list(clusters.values())
