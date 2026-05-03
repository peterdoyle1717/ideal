#!/bin/bash
# run_proof_doob.sh [V1 [V2]]
#
# Run sub/supersolution proof checker for prime 6-nets on doob.
# Reads  $PRIME_DIR/N.txt (or $PRIME_DIR/N.txt.gz)
# Writes ideal/proof_N.bin  (1 float32 per net, ordered)
#
# Each record: 1 float32 — slack/eps (>0 = proved, NaN = solver failed).
# Record i matches line i of prime/N.txt exactly.
#
# Usage:
#   ./scripts/run_proof_doob.sh           # all v=4..80
#   ./scripts/run_proof_doob.sh 80        # single v
#   ./scripts/run_proof_doob.sh 81 100    # extend range
#
# Requires:
#   src/ideal_proof compiled  (run: make)
#   ../clers/python/clers_decode.py
#   $PRIME_DIR/N.txt[.gz]  (default: ../clers/fuller/prime)

set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
PROOF="$ROOT/src/ideal_proof"
DECODE="$(cd "$ROOT/../clers/python" && pwd)/clers_decode.py"
PRIME="${PRIME_DIR:-$(cd "$ROOT/../clers/fuller/prime" && pwd)}"
TMP="$ROOT/tmp_proof"
JOBS=80

[ -x "$PROOF" ]    || { echo "ERROR: $PROOF not found — run: make"; exit 1; }
[ -f "$DECODE" ]   || { echo "ERROR: $DECODE not found"; exit 1; }

mkdir -p "$TMP"

proof_one() {
    local v=$1
    local dst="$ROOT/proof_${v}.bin"

    [ -f "$dst" ] && { echo "v=$v: $dst already exists, skipping"; return 0; }

    # Accept plain .txt or gzipped .txt.gz
    local src zcat
    if   [ -f "$PRIME/${v}.txt.gz" ]; then src="$PRIME/${v}.txt.gz"; zcat="zcat"
    elif [ -f "$PRIME/${v}.txt"    ]; then src="$PRIME/${v}.txt";    zcat="cat"
    else echo "ERROR: $PRIME/${v}.txt[.gz] not found"; return 1; fi

    local nlines
    nlines=$($zcat "$src" | wc -l)
    [ "$nlines" -eq 0 ] && { echo "v=$v: 0 nets — skipping"; return 0; }

    local shards=$JOBS
    [ "$nlines" -lt "$shards" ] && shards=$nlines

    echo "v=$v: $nlines nets, $shards shards → $dst"

    $zcat "$src" | python3 "$DECODE" | split -n "l/$shards" - "$TMP/pf_${v}_"

    nice -n 19 parallel -j "$JOBS" \
        "$PROOF < {} > {}.bin" \
        ::: "$TMP"/pf_${v}_*

    cat "$TMP"/pf_${v}_*.bin > "$dst"

    local nbytes nbytes_expected
    nbytes=$(wc -c < "$dst")
    nbytes_expected=$(( nlines * 4 ))
    if [ "$nbytes" -ne "$nbytes_expected" ]; then
        echo "  ERROR: expected $nbytes_expected bytes, got $nbytes — removing $dst"
        rm -f "$dst" "$TMP"/pf_${v}_* "$TMP"/pf_${v}_*.bin
        return 1
    fi

    local nnan
    nnan=$(python3 -c "
import struct, math, sys
data = open('$dst','rb').read()
n = len(data)//4
vals = struct.unpack(f'{n}f', data)
print(sum(1 for v in vals if math.isnan(v)))
")
    echo "  → proof_${v}.bin: $nlines records, $nnan NaN (unproved)"
    rm -f "$TMP"/pf_${v}_* "$TMP"/pf_${v}_*.bin
}

case $# in
    0) for v in $(seq 4 80); do proof_one "$v"; done ;;
    1) proof_one "$1" ;;
    *) for v in $(seq "$1" "$2"); do proof_one "$v"; done ;;
esac

rmdir "$TMP" 2>/dev/null || true
echo "=== done ==="
