#!/bin/bash
# prove.sh — gated entry point for the Euclidean existence prover.
#
# Wraps src/prove_c with a mandatory CGAL embed_check on input. OBJs
# that fail embed_check are reported as FAIL not_embedded and never
# reach the Kantorovich math; OBJs that pass go on to prove_c.
#
# This is the "lady and bathtub" boundary check — even when input
# OBJs come from puffup_c (which has its own dent gate), the prover
# treats every OBJ as untrusted.
#
# Usage:
#   scripts/prove.sh OBJDIR
#   scripts/prove.sh OBJ.obj

set -eu

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
EMBED="$ROOT/src/embed_check"
PROVE="$ROOT/src/prove_c"

if [ ! -x "$EMBED" ]; then
    echo "FAIL embed_check_unavailable: build with 'make src/embed_check' (needs CGAL+gmp+mpfr)" >&2
    exit 2
fi
if [ ! -x "$PROVE" ]; then
    echo "FAIL prove_c_unavailable: build with 'make src/prove_c'" >&2
    exit 2
fi
if [ $# -lt 1 ]; then
    echo "usage: prove.sh OBJDIR | OBJ.obj" >&2
    exit 1
fi

INPUT="$1"
if [ -d "$INPUT" ]; then
    OBJS=( "$INPUT"/*.obj )
elif [ -f "$INPUT" ]; then
    OBJS=( "$INPUT" )
else
    echo "not found: $INPUT" >&2
    exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/passed"

# embed_check prints one line per input: "1 stem" embedded, "0 stem" not.
# Output is in input order, so we can zip OBJS[] with embed.txt line-by-line.
"$EMBED" "${OBJS[@]}" > "$TMP/embed.txt" 2> "$TMP/embed.err" || true

n_emb=0
n_si=0
fail_lines="$TMP/fail_lines.txt"
: > "$fail_lines"

i=0
while read -r v _stem; do
    [ -z "${v:-}" ] && continue
    obj="${OBJS[$i]}"
    name="$(basename "$obj" .obj)"
    if [ "$v" = "1" ]; then
        ln -s "$(cd "$(dirname "$obj")" && pwd)/$(basename "$obj")" "$TMP/passed/$name.obj"
        n_emb=$((n_emb + 1))
    else
        printf "  %-40s   ?   ?   ?  FAIL  not_embedded (CGAL)\n" "$name" >> "$fail_lines"
        n_si=$((n_si + 1))
    fi
    i=$((i + 1))
done < "$TMP/embed.txt"

# Print FAIL not_embedded lines first, then run prove_c on the survivors.
[ -s "$fail_lines" ] && cat "$fail_lines"

if [ "$n_emb" -gt 0 ]; then
    "$PROVE" "$TMP/passed"
else
    printf "# pass=0  fail=%d\n" "$n_si"
fi
