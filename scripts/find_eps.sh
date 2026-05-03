#!/bin/bash
# find_eps.sh V
#
# Find the smallest eps (from 1/500, 1/1000, ...) that proves all prime 6-nets
# at vertex count V.  Prints one line: "v=V: eps=1/N (M nets, all proved)".
#
# Usage:
#   ./scripts/find_eps.sh 25
#   PRIME_DIR=~/neo/clers/fuller/prime ./scripts/find_eps.sh 60
#
# Requires: src/ideal_proof compiled (run: make), ../clers/python/clers_decode.py

set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
PROOF="$ROOT/src/ideal_proof"
DECODE="$(cd "$ROOT/../clers/python" && pwd)/clers_decode.py"
PRIME="${PRIME_DIR:-$(cd "$ROOT/../clers/fuller/prime" && pwd)}"
TMP="$ROOT/tmp_findeps"
JOBS=80

[ -x "$PROOF" ]  || { echo "ERROR: $PROOF not found — run: make"; exit 1; }
[ -f "$DECODE" ] || { echo "ERROR: $DECODE not found"; exit 1; }
[ -n "$1" ]      || { echo "Usage: $0 V"; exit 1; }

V=$1

local_src zcat
if   [ -f "$PRIME/${V}.txt.gz" ]; then local_src="$PRIME/${V}.txt.gz"; zcat="zcat"
elif [ -f "$PRIME/${V}.txt"    ]; then local_src="$PRIME/${V}.txt";    zcat="cat"
else echo "ERROR: $PRIME/${V}.txt[.gz] not found"; exit 1; fi

nlines=$($zcat "$local_src" | wc -l)

if [ "$nlines" -eq 0 ]; then
    echo "v=$V: 0 nets — trivially proved"
    exit 0
fi

mkdir -p "$TMP"

shards=$JOBS
[ "$nlines" -lt "$shards" ] && shards=$nlines

count_nan() {
    python3 -c "
import sys, struct, math
d = sys.stdin.buffer.read()
n = len(d) // 4
print(sum(1 for x in struct.unpack(f'{n}f', d) if math.isnan(x)))
"
}

for denom in 500 1000 2000 4000 8000 16000 32000; do
    eps=$(python3 -c "print(f'{1/$denom:.10g}')")

    $zcat "$local_src" | python3 "$DECODE" | split -n "l/$shards" - "$TMP/fe_${V}_"

    nice -n 19 parallel -j "$JOBS" \
        "$PROOF $eps < {} > {}.out" \
        ::: "$TMP"/fe_${V}_* 2>/dev/null

    nfail=$(cat "$TMP"/fe_${V}_*.out | count_nan)

    rm -f "$TMP"/fe_${V}_* "$TMP"/fe_${V}_*.out

    if [ "$nfail" -eq 0 ]; then
        echo "v=$V: eps=1/$denom ($nlines nets, all proved)"
        rmdir "$TMP" 2>/dev/null || true
        exit 0
    else
        echo "v=$V: eps=1/$denom — $nfail/$nlines failed"
    fi
done

echo "v=$V: UNPROVED — no eps in sequence worked"
rmdir "$TMP" 2>/dev/null || true
exit 1
