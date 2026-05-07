#!/usr/bin/env python3
"""Compare ideal vs vertexwish starts via puffup_c_lm at multiple α.

For each prime CLERS in v_min..v_max and each α in alphas_deg, runs the C
LM solver twice (--start ideal and --start vertexwish), parses bend
output, and reports per-(v, α): how many succeeded under each, and how
many had max|Δ| of bend vectors below SAME_THRESH.
"""
import argparse
import math
import re
import subprocess
import sys
from pathlib import Path

PUFFUP = Path("/Users/doyle/Dropbox/neo/ideal/src/puffup_c_lm")
CLERS_BIN = "/Users/doyle/Dropbox/neo/clers/bin/clers"
PRIMES = Path("/Users/doyle/Dropbox/neo/data/primes")

SAME_THRESH = 1e-6


def decode_netcode(clers: str) -> str:
    return subprocess.run(
        [CLERS_BIN, "decode"], input=clers + "\n",
        capture_output=True, text=True, check=True).stdout.strip()


def run_one(netcode: str, alpha_deg: float, start: str):
    p = subprocess.run(
        [str(PUFFUP), "--alpha-deg", str(alpha_deg), "--start", start, "--print-bends"],
        input=netcode, capture_output=True, text=True,
    )
    out = {"success": False, "iters": None, "resid": None, "bends": {},
           "msg": "", "base_face": None}
    for line in p.stdout.splitlines():
        if line.startswith("success: "):
            out["success"] = (line.split()[1] == "true")
        elif line.startswith("iters: "):
            out["iters"] = int(line.split()[1])
        elif line.startswith("final_resid: "):
            out["resid"] = float(line.split()[1])
        elif line.startswith("message: "):
            out["msg"] = line[len("message: "):]
        elif line.startswith("base_face: "):
            out["base_face"] = tuple(int(x) for x in line.split()[1:])
        elif line.startswith("bend "):
            parts = line.split()
            edge = parts[1]
            val = float(parts[2])
            out["bends"][edge] = val
    return out


def base_edges(base_face):
    a, b, c = base_face
    def k(u, v): return f"{min(u,v)}-{max(u,v)}"
    return {k(a, b), k(b, c), k(c, a)}


def max_abs_diff(bw: dict, bi: dict, exclude: set = None) -> float:
    """Max |Δ| over edges in BOTH dicts, optionally excluding `exclude`."""
    exclude = exclude or set()
    keys = (bw.keys() | bi.keys()) - exclude
    if not keys:
        return float("inf")
    m = 0.0
    for k in keys:
        if k not in bw or k not in bi:
            return float("inf")
        d = abs(bw[k] - bi[k])
        if d > m:
            m = d
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v-min", type=int, default=20)
    ap.add_argument("--v-max", type=int, default=26)
    ap.add_argument("--alphas-deg", type=float, nargs="+", default=[1.0, 2.5, 5.0])
    ap.add_argument("--limit", type=int, default=0,
                    help="if >0, cap CLERSes per v to this count")
    args = ap.parse_args()

    print(f"v={args.v_min}..{args.v_max}  alphas={args.alphas_deg}  limit={args.limit}", flush=True)
    print(flush=True)
    print(f"  {'v':>3}  {'α°':>5}  "
          f"{'tot':>5}  {'idl_ok':>6}  {'wsh_ok':>6}  "
          f"{'both_ok':>7}  {'same':>5}  {'diff':>5}  "
          f"{'idl_iters':>9}  {'wsh_iters':>9}", flush=True)

    for v in range(args.v_min, args.v_max + 1):
        path = PRIMES / f"{v}.txt"
        if not path.exists():
            continue
        clersts = [c.strip() for c in path.read_text().splitlines() if c.strip()]
        if args.limit > 0:
            clersts = clersts[: args.limit]

        netcodes = []
        for c in clersts:
            try:
                netcodes.append(decode_netcode(c))
            except Exception:
                netcodes.append(None)

        for alpha_deg in args.alphas_deg:
            n_total = 0
            n_idl_ok = 0
            n_wsh_ok = 0
            n_both = 0
            n_same = 0
            n_diff = 0
            sum_idl_iters = 0
            sum_wsh_iters = 0
            for nc in netcodes:
                if nc is None:
                    continue
                n_total += 1
                ri = run_one(nc, alpha_deg, "ideal")
                rw = run_one(nc, alpha_deg, "vertexwish")
                if ri["success"]:
                    n_idl_ok += 1
                    sum_idl_iters += ri["iters"] or 0
                if rw["success"]:
                    n_wsh_ok += 1
                    sum_wsh_iters += rw["iters"] or 0
                if ri["success"] and rw["success"]:
                    n_both += 1
                    excl = base_edges(ri["base_face"]) if ri["base_face"] else set()
                    d = max_abs_diff(rw["bends"], ri["bends"], exclude=excl)
                    if d < SAME_THRESH:
                        n_same += 1
                    else:
                        n_diff += 1
            avg_idl = sum_idl_iters / n_idl_ok if n_idl_ok else 0
            avg_wsh = sum_wsh_iters / n_wsh_ok if n_wsh_ok else 0
            print(f"  {v:>3}  {alpha_deg:>5.1f}  "
                  f"{n_total:>5}  {n_idl_ok:>6}  {n_wsh_ok:>6}  "
                  f"{n_both:>7}  {n_same:>5}  {n_diff:>5}  "
                  f"{avg_idl:>9.2f}  {avg_wsh:>9.2f}",
                  flush=True)


if __name__ == "__main__":
    main()
