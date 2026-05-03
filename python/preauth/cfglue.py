"""cfglue.py — bigface-based donor/receiver records and gluing.

See DONORS.md for the design. A bigface is a maximal coplanar-adjacent
patch of ≥ 2 triangles in the flop; it's a triangulated disk. Donor
and receiver offers are keyed by the CLERS encoding of the bigface
(encode with boundary_cycle), so string equality = glue-compatible
shape. Corner-degree check applies position-by-position.
"""
from __future__ import annotations
import sys, os
from dataclasses import dataclass, field
from typing import Optional

from . import _CLERS_SRC  # noqa: F401  (canonical clers on path)
from .cluster_sig import face_clusters
from clers import encode as _encode_closed, reverse_poly, official_unoriented_name, build_edge_map


def encode(poly, start=None, boundary_cycle=None):
    """Encode a triangulated surface as a CLERS string.

    Without `boundary_cycle`: closed-sphere encoding — defers to canonical
    `clers.encode`.

    With `boundary_cycle`: open-disk encoding (used here for bigface CLERS).
    Pass the boundary as a list of vertex indices in CCW order starting
    with (a, b, …). Same main loop as the closed case; only initial seed
    differs. Local copy of the loop because canonical clers doesn't expose
    a boundary-seeded variant.
    """
    if boundary_cycle is None:
        return _encode_closed(poly, start=start)

    em = build_edge_map(poly)
    a0, b0 = (start if start is not None else (poly[0][0], poly[0][1]))
    assert boundary_cycle[0] == a0 and boundary_cycle[1] == b0, \
        "boundary_cycle must start with (a0, b0)"

    vset = set(boundary_cycle)
    tree = set()
    n = len(boundary_cycle)
    for i in range(1, n):
        tree.add((boundary_cycle[i], boundary_cycle[(i + 1) % n]))

    code = []
    stack = [(a0, b0)]
    while stack:
        x, y = stack.pop()
        z = em[(x, y)]
        if z not in vset:
            code.append("C")
            tree.add((z, y))
            vset.add(z)
            stack.append((x, z))
        elif (y, z) in tree and (z, x) in tree:
            code.append("E")
        elif (z, x) in tree:
            code.append("A")
            stack.append((z, y))
        elif (y, z) in tree:
            code.append("B")
            stack.append((x, z))
        else:
            code.append("D")
            stack.append((z, y))
            stack.append((x, z))
    return "".join(code)


# -----------------------------------------------------------------
# Bigface extraction
# -----------------------------------------------------------------

@dataclass
class Bigface:
    source: str                  # parent flop name
    faces: tuple                 # bigface triangles as (u, v, w) tuples (parent's indices)
    ccw_boundary: tuple          # vertex indices walked CCW around the boundary
    ambient_deg: dict            # {v -> parent flop edge-degree} (covers bigface verts)
    cluster_deg: dict            # {v -> count of bigface triangles at v}


def _boundary_cycle_ccw(cluster, flop_faces):
    """Walk boundary CCW using face winding. Returns tuple of vertex indices, or None."""
    edge_count = {}
    for fi in cluster:
        f = flop_faces[fi]
        for i in range(3):
            u, v = f[i], f[(i + 1) % 3]
            e = (min(u, v), max(u, v))
            edge_count[e] = edge_count.get(e, 0) + 1
    boundary_edges = {e for e, c in edge_count.items() if c == 1}
    nxt = {}
    for fi in cluster:
        f = flop_faces[fi]
        for i in range(3):
            u, v = f[i], f[(i + 1) % 3]
            if (min(u, v), max(u, v)) in boundary_edges:
                if u in nxt:
                    return None
                nxt[u] = v
    if not nxt:
        return None
    start = next(iter(nxt))
    cycle = [start]
    cur = nxt[start]
    while cur != start:
        cycle.append(cur)
        if cur not in nxt:
            return None
        cur = nxt[cur]
    if len(cycle) != len(nxt):
        return None
    return tuple(cycle)


def _ambient_degrees(flop_faces):
    nbrs = {}
    for f in flop_faces:
        for u, v in [(f[0], f[1]), (f[1], f[2]), (f[0], f[2])]:
            nbrs.setdefault(u, set()).add(v)
            nbrs.setdefault(v, set()).add(u)
    return {v: len(ns) for v, ns in nbrs.items()}


def _per_vertex_count(faces):
    cd = {}
    for f in faces:
        for v in f:
            cd[v] = cd.get(v, 0) + 1
    return cd


def extract_bigfaces(flop_verts, flop_faces, source='?'):
    """Return list of Bigface objects for the flop (drops single-face clusters)."""
    clusters = face_clusters(flop_verts, flop_faces)
    clusters = [c for c in clusters if len(c) > 1]
    amb = _ambient_degrees(flop_faces)

    bigfaces = []
    for cluster in clusters:
        cycle = _boundary_cycle_ccw(cluster, flop_faces)
        if cycle is None:
            continue
        faces = tuple(tuple(flop_faces[fi]) for fi in cluster)
        cdeg = _per_vertex_count(faces)
        bf_amb = {v: amb[v] for v in cdeg}
        bigfaces.append(Bigface(
            source=source,
            faces=faces,
            ccw_boundary=cycle,
            ambient_deg=bf_amb,
            cluster_deg=cdeg,
        ))
    return bigfaces


# -----------------------------------------------------------------
# Donor and receiver records
# -----------------------------------------------------------------

@dataclass
class Donor:
    bigface_code: str            # CLERS of bigface disk from (a, b) CCW
    ambient_degs: tuple          # parent-flop edge-degrees along CCW walk from a
    cluster_degs: tuple          # bigface triangle counts at each boundary vertex (same walk)
    iso_key: str                 # full-polyhedron CLERS starting at (a, b) — dedup within flop
    bigface: Bigface
    start_vertex: int


@dataclass
class Receiver:
    bigface_code: str            # CLERS of reverse-bigface disk from (a, b) CW
    ambient_degs: tuple          # parent-flop edge-degrees along CW walk from a
    cluster_degs: tuple
    iso_key: str                 # full-polyhedron CLERS of reverse_poly starting at (a, b)
    bigface: Bigface
    start_vertex: int


def donor_offers(bigface: Bigface, flop_faces):
    """Yield Donor per distinct boundary vertex, deduplicated by iso_key
    (whole-polyhedron CLERS), so that symmetric equivalents collapse."""
    seen_iso = set()
    ccw = bigface.ccw_boundary
    n = len(ccw)
    for start_idx in range(n):
        rotated = ccw[start_idx:] + ccw[:start_idx]
        a, b = rotated[0], rotated[1]
        iso_key = encode(flop_faces, start=(a, b))
        if iso_key in seen_iso:
            continue
        seen_iso.add(iso_key)
        bf_code = encode(bigface.faces, start=(a, b), boundary_cycle=list(rotated))
        amb = tuple(bigface.ambient_deg[v] for v in rotated)
        clust = tuple(bigface.cluster_deg.get(v, 0) for v in rotated)
        yield Donor(
            bigface_code=bf_code,
            ambient_degs=amb,
            cluster_degs=clust,
            iso_key=iso_key,
            bigface=bigface,
            start_vertex=a,
        )


def receiver_offers(bigface: Bigface, flop_faces):
    """Yield Receiver per distinct boundary vertex (CW walk), deduplicated by
    iso_key of the reverse-oriented polyhedron."""
    seen_iso = set()
    ccw = bigface.ccw_boundary
    cw = tuple(reversed(ccw))
    rev_flop = reverse_poly(list(flop_faces))
    rev_bf_faces = reverse_poly(list(bigface.faces))
    n = len(cw)
    for start_idx in range(n):
        rotated = cw[start_idx:] + cw[:start_idx]
        a, b = rotated[0], rotated[1]
        iso_key = encode(rev_flop, start=(a, b))
        if iso_key in seen_iso:
            continue
        seen_iso.add(iso_key)
        bf_code = encode(rev_bf_faces, start=(a, b), boundary_cycle=list(rotated))
        amb = tuple(bigface.ambient_deg[v] for v in rotated)
        clust = tuple(bigface.cluster_deg.get(v, 0) for v in rotated)
        yield Receiver(
            bigface_code=bf_code,
            ambient_degs=amb,
            cluster_degs=clust,
            iso_key=iso_key,
            bigface=bigface,
            start_vertex=a,
        )


# -----------------------------------------------------------------
# Matching
# -----------------------------------------------------------------

def corner_degree_test(donor: Donor, receiver: Receiver) -> bool:
    if len(donor.ambient_degs) != len(receiver.ambient_degs):
        return False
    for i in range(len(donor.ambient_degs)):
        final = donor.ambient_degs[i] + receiver.ambient_degs[i] - 2 * donor.cluster_degs[i]
        if final not in (4, 5, 6):
            return False
    return True


def find_matches(inventory):
    """inventory: list of (source, flop_faces, [bigface, ...]).

    Returns list of (donor, receiver, donor_flop_faces, recv_flop_faces)."""
    recv_by_code = {}
    for source, flop_faces, bfs in inventory:
        for bf in bfs:
            for r in receiver_offers(bf, flop_faces):
                recv_by_code.setdefault(r.bigface_code, []).append((r, flop_faces))

    matches = []
    for source, flop_faces, bfs in inventory:
        for bf in bfs:
            for d in donor_offers(bf, flop_faces):
                for r, r_flop_faces in recv_by_code.get(d.bigface_code, []):
                    if corner_degree_test(d, r):
                        matches.append((d, r, flop_faces, r_flop_faces))
    return matches


# -----------------------------------------------------------------
# Gluing
# -----------------------------------------------------------------

def glue(donor: Donor, receiver: Receiver, d_flop_faces, r_flop_faces):
    """Combine donor's flop + receiver's flop along matched bigfaces.

    Returns 1-indexed face list of the glued polyhedron."""
    # Donor's CCW boundary cycle from its start
    d_ccw = donor.bigface.ccw_boundary
    d_idx = d_ccw.index(donor.start_vertex)
    d_cycle = d_ccw[d_idx:] + d_ccw[:d_idx]

    # Receiver's CW boundary cycle from its start
    r_cw = tuple(reversed(receiver.bigface.ccw_boundary))
    r_idx = r_cw.index(receiver.start_vertex)
    r_cycle = r_cw[r_idx:] + r_cw[:r_idx]

    assert len(d_cycle) == len(r_cycle)

    # Boundary and interior vertex sets
    d_bdry = set(d_ccw)
    d_big_verts = {v for f in donor.bigface.faces for v in f}
    d_interior = d_big_verts - d_bdry  # vanishes on glue

    r_bdry = set(receiver.bigface.ccw_boundary)
    r_big_verts = {v for f in receiver.bigface.faces for v in f}
    r_interior = r_big_verts - r_bdry  # vanishes on glue

    # Receiver → donor vertex map: identified boundary vertices
    r_to_d = {r_cycle[i]: d_cycle[i] for i in range(len(d_cycle))}

    # Non-bigface-face triangles of donor flop — keep with donor's labels
    d_big_face_set = {tuple(f) for f in donor.bigface.faces}
    new_faces = []
    for f in d_flop_faces:
        if tuple(f) in d_big_face_set:
            continue
        new_faces.append(tuple(f))

    # Non-bigface-face triangles of receiver flop — relabel
    r_big_face_set = {tuple(f) for f in receiver.bigface.faces}
    # Offset for receiver-only vertices (neither matched nor interior-bigface)
    d_used = {v for f in new_faces for v in f}
    offset = (max(d_used) if d_used else 0) + 1
    for f in r_flop_faces:
        if tuple(f) in r_big_face_set:
            continue
        if any(v in r_interior for v in f):
            # Shouldn't happen: interior bigface verts are only in bigface faces
            continue
        relabeled = []
        for v in f:
            if v in r_to_d:
                relabeled.append(r_to_d[v])
            else:
                # Fresh receiver-only vertex
                if v not in r_to_d:
                    r_to_d[v] = offset
                    offset += 1
                relabeled.append(r_to_d[v])
        new_faces.append(tuple(relabeled))

    # Compact 1-indexed relabeling
    used = sorted({v for f in new_faces for v in f})
    lbl = {v: i + 1 for i, v in enumerate(used)}
    return [tuple(lbl[v] for v in f) for f in new_faces]


# -----------------------------------------------------------------
# Pancake / valid-flop utilities
# -----------------------------------------------------------------

def is_pancake(verts, faces):
    """True if the flop is a pancake: a single cluster covering the entire
    surface (every cluster edge shared by 2 cluster triangles)."""
    cs = face_clusters(verts, faces)
    cs = [c for c in cs if len(c) > 1]
    if len(cs) != 1 or len(cs[0]) != len(faces):
        return False
    edge_count = {}
    for fi in cs[0]:
        f = faces[fi]
        for u, v in [(f[0], f[1]), (f[1], f[2]), (f[0], f[2])]:
            e = (min(u, v), max(u, v))
            edge_count[e] = edge_count.get(e, 0) + 1
    return all(c == 2 for c in edge_count.values())


def degree_set(faces):
    nbrs = {}
    for f in faces:
        for u, v in [(f[0], f[1]), (f[1], f[2]), (f[0], f[2])]:
            nbrs.setdefault(u, set()).add(v)
            nbrs.setdefault(v, set()).add(u)
    return sorted(len(s) for s in nbrs.values())


def valid_flop(faces):
    V = len({v for f in faces for v in f})
    F = len(faces)
    if V < 4 or F != 2 * V - 4:
        return False
    degs = degree_set(faces)
    return degs and degs[0] >= 4 and degs[-1] <= 6




if __name__ == '__main__':
    main(sys.argv)
