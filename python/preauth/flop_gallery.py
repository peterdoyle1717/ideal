"""Donor enumerators for the preauthorizer.

Extracted from undented/floppy/src/gen_flop_gallery.py — only the
combinatorial donor enumerators (pancakes, octahedra, octahedron donors,
and the compound-from-donors gluing). HTML rendering, OBJ-based
loser_donors, and the CLI driver are dropped because the preauthorizer
is OBJ-free and runs purely on synthetic constructions.
"""
from __future__ import annotations
from collections import deque
import numpy as np

from . import _CLERS_SRC  # noqa: F401  (canonical clers on path)
from .cfglue import extract_bigfaces, find_matches, glue, valid_flop
from .gen_pancakes import (enumerate_pancakes, pancake_facelist,
                           edge_lengths, canonical_sextuple)
from .boundary_tri import (canonical_state, neighbors, boundary_triangulation,
                           vertex_set_of_triangulation, has_no_degree3,
                           has_face_interior_degree6)
from clers import official_unoriented_name


def compute_pancakes(vmax):
    sigs = {}
    for k, l, r, q, v in enumerate_pancakes(vmax):
        sig = canonical_sextuple(edge_lengths(k, l, r, q))
        if sig not in sigs:
            sigs[sig] = (k, l, r, q)
    codes = set()
    for (k, l, r, q) in sigs.values():
        res = pancake_facelist(k, l, r, q)
        if res is None:
            continue
        V, faces = res
        codes.add(official_unoriented_name(faces))
    return codes


def _enumerate_octahedra(vmax, require_no_deg3):
    """Yield (clers, verts_np, faces_1idx) for each canonical octahedron <= vmax."""
    fmax = 2 * vmax - 4
    start = canonical_state(1, 1, 1, 0, 1)
    dq = deque([start]); seen = {start}
    emitted_codes = set()
    while dq:
        s = dq.popleft()
        tris = boundary_triangulation(*s)
        if len(tris) > fmax:
            continue
        verts = vertex_set_of_triangulation(*s)
        V = len(verts)
        ok = (len(tris) == 2 * V - 4 and V <= vmax
              and has_face_interior_degree6(*s)
              and (not require_no_deg3 or has_no_degree3(*s)))
        if ok:
            vmap1 = {w: i + 1 for i, w in enumerate(verts)}
            vmap0 = {w: i for i, w in enumerate(verts)}
            poly1 = [(vmap1[a], vmap1[b], vmap1[c]) for a, b, c in tris]
            poly0 = [(vmap0[a], vmap0[b], vmap0[c]) for a, b, c in tris]
            code = official_unoriented_name(poly1)
            if code not in emitted_codes:
                emitted_codes.add(code)
                verts_np = np.array([list(v) for v in verts], dtype=float)
                yield code, verts_np, poly0
        for nxt in neighbors(s):
            if nxt not in seen:
                seen.add(nxt); dq.append(nxt)


def compute_octahedra(vmax):
    return {code for code, _, _ in _enumerate_octahedra(vmax, require_no_deg3=True)}


def octahedron_donors(vmax):
    """Return [(clers, verts, faces), ...] for deg-3-allowed octahedra <= vmax.

    Octahedra with a deg-3 tip never appear as solver-rejects, but they can
    serve as donors in gluing — when the bigface covers the deg-3 tip, the
    glued compound can be a flop the prover wouldn't otherwise see."""
    return list(_enumerate_octahedra(vmax, require_no_deg3=False))


def compute_compounds_from_donors(donors):
    """Donors: iterable of (code, verts, faces). Glue every donor pair via
    matching bigfaces, return the set of resulting valid-flop CLERSes."""
    inv = []
    for code, verts, faces in donors:
        bfs = extract_bigfaces(verts, faces, source=code)
        if bfs:
            inv.append((code, faces, bfs))
    matches = find_matches(inv)
    codes = set()
    for d, r, df, rf in matches:
        glued = glue(d, r, df, rf)
        if not valid_flop(glued):
            continue
        codes.add(official_unoriented_name(glued))
    return codes
