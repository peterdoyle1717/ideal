# ideal

Ideal horoball packings of triangulated spheres: solvers, placements, and bracket proofs.

This repository contains small C and Python tools for computing ideal horoball packings of triangulated spheres and for certifying existence by a sub/supersolution bracket. It is a project tool rather than a polished package.

Input is a face list

    a,b,c;d,e,f;...

one triangulation per line.

## Scope

The underlying geometric model is the usual product metric: adjacent vertices `i` and `j` have edge length `u[i] u[j]`, with one vertex sent to infinity and its neighbors pinned on the boundary. The unknowns are the remaining positive vertex weights.

This repository has two main jobs:

- compute the weights and associated flat placement
- certify existence by a bracket argument with explicit slack checks

The motivating combinatorial inputs in the broader project are prime 6-nets, but the basic tools take triangulated-sphere face lists as input.

## Layout

    ideal/
    ├── Makefile
    ├── data/
    │   └── eps_needed.txt
    ├── python/
    │   ├── horou.py
    │   ├── horoz.py
    │   └── proof.py
    ├── scripts/
    │   ├── run_proof_doob.sh
    │   └── find_eps.sh
    └── src/
        ├── horou_c.c
        ├── horoz_c.c
        └── proof_c.c

## Build

    make
    make clean

Requires a C compiler and `libm` (standard on macOS and Linux).

## Main tools

### horou

Solve for the horoball weights `u[v]` of a triangulation.

C version:

    ../clers/bin/clers decode < 20.txt | src/horou_c > horou/20.bin

Python version:

    echo "CCAE" | ../clers/bin/clers decode | python3 python/horou.py

### horoz

Solve for the weights and place vertices in the upper half-plane.

C version:

    ../clers/bin/clers decode < 20.txt | src/horoz_c > horoz/20.bin

Python usage:

    from horou import horou
    from horoz import horoz
    u = horou(poly)
    pos = horoz(poly, u)

### proof

Certify existence by bracketing the solution between a sub-solution and a super-solution.

The proof code checks five positive-slack conditions:

- mono
- excess
- triangle
- convex
- boundary

Python examples:

    echo "CCAE" | ../clers/bin/clers decode | python3 python/proof.py --verbose

    ../clers/bin/clers decode < 20.txt | python3 python/proof.py

C example:

    ../clers/bin/clers decode < 20.txt | src/proof_c > proof_20.bin

## Results currently checked in

The checked-in table `data/eps_needed.txt` records the smallest tested `eps` proving all prime 6-nets for each `v` through `v = 80`.

In particular:

- the checked-in results currently go through `v = 80`
- `1/8000` suffices from `v = 60` through `v = 80`

## Parallel scripts

For server-side runs:

    ./scripts/run_proof_doob.sh 81
    ./scripts/run_proof_doob.sh 81 100
    ./scripts/find_eps.sh 81

These are project scripts for batch runs, not general installation machinery.

## Notes

This repository sits in the circle-packing / discrete-conformal / ideal-polyhedral orbit of ideas, but it is not meant to be a survey or exposition. The point is to keep the computational core small, explicit, and inspectable.

Where the separate `clers` repository handles combinatorial naming and decoding of triangulations, this repository handles the ideal horoball-packing side.

## Provenance

The code and documentation in this repository were drafted primarily with Claude Code under the author's direction, with additional advice, review, and supervision from ChatGPT.

## License

See `LICENSE`.
