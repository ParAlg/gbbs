#!/usr/bin/env python3
"""Sweep HaMIS vs RG-rootset vs RG-specfor over the available graph datasets.

Reports best ### Running Time over 5 rounds (single binary invocation; gbbs
loads the graph once and runs the algorithm N times). Decompression time for
compressed graphs is excluded from the reported figure for HaMIS (printed
separately by the runner).
"""

import os
import re
import subprocess
import sys
import time

GRAPH_DIR = "/mnt/fast_data/graphdata"
ROUNDS = 5

# (display_name, filename, fmt) where fmt in {"adj","bin"}.
# .bytepda graphs are excluded because HaMIS uses indexed get_neighbor(i),
# which decompresses from scratch on every call.
GRAPHS = [
    ("friendster",      "friendster_sym.adj",             "adj"),
    ("Germany",         "Germany_sym.bin",                "bin"),
    ("RoadUSA",         "RoadUSA_sym.bin",                "bin"),
    ("GeoLife_2",       "GeoLifeNoScale_2_sym.bin",       "bin"),
    ("GeoLife_5",       "GeoLifeNoScale_5_sym.bin",       "bin"),
    ("GeoLife_10",      "GeoLifeNoScale_10_sym.bin",      "bin"),
    ("GeoLife_15",      "GeoLifeNoScale_15_sym.bin",      "bin"),
    ("GeoLife_20",      "GeoLifeNoScale_20_sym.bin",      "bin"),
]

BIN_HAMIS = "./bazel-bin/benchmarks/MaximalIndependentSet/HaMIS/HaMIS_main"
BIN_RG    = "./bazel-bin/benchmarks/MaximalIndependentSet/RandomGreedy/MaximalIndependentSet_main"

# Algorithms: (display, binary, extra-args)
ALGS = [
    ("HaMIS",       BIN_HAMIS, []),
    ("RG-rootset",  BIN_RG,    []),
    ("RG-specfor",  BIN_RG,    ["-specfor"]),
]


def format_flags(fmt):
    if fmt == "adj":     return ["-s"]
    if fmt == "bin":     return ["-s", "-b"]
    if fmt == "bytepda": return ["-s", "-c"]
    raise ValueError(fmt)


def run_one(alg_name, binary, extra_args, graph_path, fmt, rounds):
    cmd = [binary] + format_flags(fmt) + extra_args
    cmd += ["-stats", "-rounds", str(rounds)]
    # HaMIS needs to decompress in-memory for bytepda inputs.
    if fmt == "bytepda" and alg_name == "HaMIS":
        cmd += ["-decompress"]
    cmd += [graph_path]
    t0 = time.time()
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        return {"error": "timeout", "wall": time.time() - t0}
    wall = time.time() - t0
    stdout = out.stdout
    stderr = out.stderr
    if out.returncode != 0:
        return {"error": f"rc={out.returncode}", "wall": wall,
                "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:]}
    times = [float(m) for m in re.findall(r"### Running Time:\s+([0-9.e+\-]+)", stdout)]
    sizes = [int(m) for m in re.findall(r"MaximalIndependentSet size:\s+(\d+)", stdout)]
    decomp = re.findall(r"### Decompression Time:\s+([0-9.e+\-]+)", stdout)
    n_match = re.search(r"### n:\s+(\d+)", stdout)
    m_match = re.search(r"### m:\s+(\d+)", stdout)
    return {
        "times": times,
        "best": min(times) if times else None,
        "mis_size": sizes[0] if sizes else None,
        "mis_sizes_consistent": len(set(sizes)) <= 1 if sizes else None,
        "decompress": float(decomp[0]) if decomp else None,
        "n": int(n_match.group(1)) if n_match else None,
        "m": int(m_match.group(1)) if m_match else None,
        "wall": wall,
    }


def main():
    print(f"### Sweep: {len(GRAPHS)} graphs x {len(ALGS)} algorithms, {ROUNDS} rounds each",
          flush=True)
    results = {}  # (graph_name, alg_name) -> result
    for gname, fname, fmt in GRAPHS:
        path = os.path.join(GRAPH_DIR, fname)
        if not os.path.exists(path):
            print(f"!! missing graph: {path}", flush=True)
            continue
        for alg_name, binary, extra in ALGS:
            print(f"\n>>> {gname} :: {alg_name}", flush=True)
            r = run_one(alg_name, binary, extra, path, fmt, ROUNDS)
            results[(gname, alg_name)] = r
            if "error" in r:
                print(f"    ERROR: {r['error']} (wall {r['wall']:.1f}s)", flush=True)
                if "stderr_tail" in r:
                    print(r["stderr_tail"], flush=True)
            else:
                print(f"    n={r['n']} m={r['m']} mis_size={r['mis_size']} "
                      f"best={r['best']:.6f}s times={r['times']} wall={r['wall']:.1f}s",
                      flush=True)
                if r.get("decompress") is not None:
                    print(f"    decompression: {r['decompress']:.4f}s", flush=True)

    # Final table
    print("\n\n=========== SUMMARY TABLE (best ### Running Time over rounds) ===========")
    header = f"{'graph':<16} {'n':>12} {'m':>14}  " + "  ".join(
        f"{a:>22}" for a, _, _ in ALGS)
    print(header)
    for gname, _, _ in GRAPHS:
        any_r = None
        for a, _, _ in ALGS:
            r = results.get((gname, a))
            if r and "error" not in r:
                any_r = r
                break
        if any_r is None:
            continue
        cells = []
        hamis = results.get((gname, "HaMIS"))
        for a, _, _ in ALGS:
            r = results.get((gname, a))
            if not r or "error" in r:
                cells.append("ERR")
                continue
            tag = f"{r['best']*1000:.2f}ms |{r['mis_size']}"
            if a != "HaMIS" and hamis and "error" not in hamis and hamis['best']:
                tag += f" ({r['best']/hamis['best']:.1f}x)"
            cells.append(tag)
        row = f"{gname:<16} {any_r['n']:>12} {any_r['m']:>14}  " + "  ".join(
            f"{c:>22}" for c in cells)
        print(row)


if __name__ == "__main__":
    main()
