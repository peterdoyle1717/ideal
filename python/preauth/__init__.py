"""Preauthorization — by-construction proofs of net realizability.

Identifies CLERS strings whose existence is established by explicit
geometric construction (pancakes, octahedra, thick pancakes, the v=14
hex antiprism, and compounds glued from these via matching bigfaces).

Doesn't require OBJs from a solver — pure combinatorial enumeration
plus a synthetic 3D embedding good enough for face_clusters to detect
coplanar bigfaces.
"""
import sys
from pathlib import Path

# Use the canonical clers Python module from clers/src/.
_CLERS_SRC = Path(__file__).resolve().parent.parent.parent.parent / 'clers' / 'src'
if str(_CLERS_SRC) not in sys.path:
    sys.path.insert(0, str(_CLERS_SRC))
