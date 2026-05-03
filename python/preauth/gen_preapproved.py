#!/usr/bin/env python3
"""gen_preapproved.py — emit per-v CLERS lists of preauthorized floppers.

Inductively builds the preauthorized set up to --vmax:
  pancake          double cover of lattice polygon
  octahedron       cannonball region (no deg-3)
  thick pancake    planar polygon × [0,1] hip-roof
  misc             v=14 hex antiprism (CCCCACCACACACACACAACAAAE)
  compound         glued from smaller flops along matching bigfaces

Output: <outdir>/<v>.txt, one CLERS per line.

Usage:
  python3 -m preauth.gen_preapproved --vmax 50 --outdir /tmp/preapp
"""
from __future__ import annotations
import sys, os, argparse
from collections import defaultdict

from . import _CLERS_SRC  # noqa: F401
from .gen_thick_pancakes import enumerate_polygons, build_thick_pancake
from .flop_gallery import (compute_pancakes, compute_octahedra,
                           octahedron_donors, compute_compounds_from_donors)
from clers import official_unoriented_name

MISC = {'CCCCACCACACACACACAACAAAE'}  # v=14 hex antiprism


def vof(clers):
    return (len(clers) + 4) // 2


def gen_thick_pancake_codes(vmax):
    seen = set()
    for k, l, r, q in enumerate_polygons(vmax):
        res = build_thick_pancake(k, l, r, q)
        if res is None:
            continue
        V, faces = res
        if V > vmax or len(faces) != 2 * V - 4:
            continue
        seen.add(official_unoriented_name(faces))
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vmax', type=int, default=50)
    ap.add_argument('--outdir', default='/tmp/preapp_ideal')
    ap.add_argument('--no-compound', action='store_true')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f'pancakes vmax={args.vmax} ...', file=sys.stderr)
    pancakes = compute_pancakes(args.vmax)
    print(f'  {len(pancakes)} pancakes', file=sys.stderr)

    print(f'octahedra (no deg-3) vmax={args.vmax} ...', file=sys.stderr)
    octas = compute_octahedra(args.vmax)
    print(f'  {len(octas)} octahedra', file=sys.stderr)

    print(f'thick pancakes vmax={args.vmax} ...', file=sys.stderr)
    thicks = gen_thick_pancake_codes(args.vmax)
    print(f'  {len(thicks)} thick pancakes', file=sys.stderr)

    print(f'misc: {len(MISC)} (v=14 hex antiprism)', file=sys.stderr)

    compounds = set()
    if not args.no_compound:
        print('compounds (donor pool: deg-3-allowed octahedra) ...', file=sys.stderr)
        donors = octahedron_donors(args.vmax)
        print(f'  donor pool: {len(donors)}', file=sys.stderr)
        compounds = compute_compounds_from_donors(donors)
        compounds = {c for c in compounds if vof(c) <= args.vmax}
        print(f'  {len(compounds)} compounds', file=sys.stderr)

    all_codes = pancakes | octas | thicks | MISC | compounds
    by_v = defaultdict(set)
    for c in all_codes:
        by_v[vof(c)].add(c)

    print(f'\nper-v counts (total / pancake / oct / thick / misc / compound):',
          file=sys.stderr)
    total = 0
    for v in sorted(by_v):
        codes = sorted(by_v[v])
        with open(os.path.join(args.outdir, f'{v}.txt'), 'w') as f:
            for c in codes:
                f.write(c + '\n')
        np_ = sum(1 for c in codes if c in pancakes)
        no = sum(1 for c in codes if c in octas)
        nt = sum(1 for c in codes if c in thicks)
        nm = sum(1 for c in codes if c in MISC)
        nc = sum(1 for c in codes if c in compounds)
        print(f'  v={v:3d}  total={len(codes):6d}  '
              f'p={np_:5d}  o={no:5d}  t={nt:5d}  m={nm}  c={nc:5d}',
              file=sys.stderr)
        total += len(codes)
    print(f'\ntotal CLERSes: {total}', file=sys.stderr)
    print(f'wrote per-v files to {args.outdir}', file=sys.stderr)


if __name__ == '__main__':
    main()
