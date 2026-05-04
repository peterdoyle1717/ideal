#!/bin/bash
# euclid_prover.sh — gated entry point for the Euclidean existence prover.
#
# Accepts either a single OBJ file or a directory of OBJs. Each OBJ goes
# through two stages:
#   1. src/embed_check (CGAL exact predicates) — certifies non-self-
#      intersection across all triangle pairs, including pairs sharing
#      a vertex. Failures are recorded as REJECT-not-embedded and never
#      reach the prover.
#   2. src/euclid_prover — runs only on OBJs that passed embed_check.
#      The prover's certificate-of-existence-near-input math depends on
#      this certified-embedded precondition.
#
# Usage:
#   scripts/euclid_prover.sh OBJ.obj             # single OBJ -> stdout tsv line
#   scripts/euclid_prover.sh OBJDIR [OUT.tsv]    # parallel sweep, default OUT: <OBJDIR>.tsv
#
# Env: JOBS (default 80, directory mode only)

set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EMBED="$ROOT/src/embed_check"
PROVE="$ROOT/src/euclid_prover"
[ -x "$EMBED" ] || { echo "FAIL embed_check_unavailable: build with 'make src/embed_check' (needs CGAL+gmp+mpfr)" >&2; exit 2; }
[ -x "$PROVE" ] || { echo "FAIL euclid_prover_unavailable: build with 'make src/euclid_prover'" >&2; exit 2; }
[ $# -ge 1 ] || { echo "usage: $(basename "$0") OBJ.obj | OBJDIR [OUT.tsv]" >&2; exit 1; }

INPUT="$1"

# Single-OBJ mode: gate, prove, emit one tsv line on stdout. Mirrors the
# old wrapper's single-file interface so callers that pass an .obj path
# directly still work.
if [ -f "$INPUT" ]; then
    name=$(basename "$INPUT" .obj)
    embed_flag=$("$EMBED" "$INPUT" 2>/dev/null | head -1 | awk '{print $1}')
    if [ "$embed_flag" != "1" ]; then
        printf "%s\tREJECT-not-embedded\n" "$name"
        exit 0
    fi
    verdict=$("$PROVE" "$INPUT" 2>&1 | grep -E "^final:" | tail -1 | awk '{print $2}')
    [ -z "$verdict" ] && verdict=ERROR
    printf "%s\t%s\n" "$name" "$verdict"
    exit 0
fi

[ -d "$INPUT" ] || { echo "not a file or directory: $INPUT" >&2; exit 1; }

OBJDIR="$INPUT"
OUT="${2:-$(basename "$OBJDIR").tsv}"
JOBS="${JOBS:-80}"

: > "$OUT"
ls "$OBJDIR"/*.obj 2>/dev/null | xargs -P "$JOBS" -I {} bash -c '
    obj="$1"; embed="$2"; prove="$3"; out="$4"
    name=$(basename "$obj" .obj)
    # Stage 1: CGAL embed_check gate. Output line is "1 stem" or "0 stem".
    embed_flag=$("$embed" "$obj" 2>/dev/null | head -1 | awk "{print \$1}")
    if [ "$embed_flag" != "1" ]; then
        printf "%s\tREJECT-not-embedded\n" "$name" >> "$out"
        exit 0
    fi
    # Stage 2: existence prover on the certified-embedded OBJ.
    verdict=$("$prove" "$obj" 2>&1 | grep -E "^final:" | tail -1 | awk "{print \$2}")
    [ -z "$verdict" ] && verdict=ERROR
    printf "%s\t%s\n" "$name" "$verdict" >> "$out"
' _ {} "$EMBED" "$PROVE" "$OUT"

echo "wrote $OUT  ($(wc -l < "$OUT" | tr -d ' ') OBJs)"
awk -F'\t' '{c[$2]++} END {for (k in c) printf "  %s=%d\n", k, c[k]}' "$OUT"
