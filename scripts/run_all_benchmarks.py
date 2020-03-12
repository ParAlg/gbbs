"""
This script runs all the benchmark binaries on specified input graphs.

Everything listed as a Bazel `cc_binary` under `benchmarks/` must be listed in
either `UNWEIGHTED_GRAPH_BENCHMARKS`, `WEIGHTED_GRAPH_BENCHMARKS`, or
`IGNORED_BINARIES` below. Forcing this explicit labeling stops contributors from
forgetting to update this file when they add a new benchmark.

The script only checks that the benchmarks run and exit without an error. It
does not check that the output of each benchmark is correct.

This script could be extended to further split benchmarks into ones that process
symmetric graphs versus asymmetric graphs, but currently the input graphs must
be symmetric so that all benchmarks can run on them.

This script should be invoked directly via Python 3. Because this script calls
other Bazel commands, invoking it with `bazel run` won't work.
"""
from typing import List, Optional, Set
import argparse
import fnmatch
import os
import sys
import subprocess

# The script will invoke these benchmark on an unweighted graph.
UNWEIGHTED_GRAPH_BENCHMARKS = [
    "//benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12:DensestSubgraph_main",
    "//benchmarks/ApproximateDensestSubgraph/GreedyCharikar:DensestSubgraph_main",
    "//benchmarks/ApproximateSetCover/MANISBPT11:ApproximateSetCover_main",
    "//benchmarks/BFS/NonDeterministicBFS:BFS_main",
    "//benchmarks/Biconnectivity/TarjanVishkin:Biconnectivity_main",
    "//benchmarks/CliqueCounting:Clique_main",
    "//benchmarks/Connectivity/BFSCC:Connectivity_main",
    "//benchmarks/Connectivity/Framework/mains:bfscc",
    "//benchmarks/Connectivity/Framework/mains:gbbscc",
    "//benchmarks/Connectivity/Framework/mains:jayanti_bfs",
    "//benchmarks/Connectivity/Framework/mains:jayanti_kout",
    "//benchmarks/Connectivity/Framework/mains:jayanti_ldd",
    "//benchmarks/Connectivity/Framework/mains:jayanti_nosample",
    "//benchmarks/Connectivity/Framework/mains:label_propagation",
    "//benchmarks/Connectivity/Framework/mains:liutarjan_bfs",
    "//benchmarks/Connectivity/Framework/mains:liutarjan_kout",
    "//benchmarks/Connectivity/Framework/mains:liutarjan_ldd",
    "//benchmarks/Connectivity/Framework/mains:liutarjan_nosample",
    "//benchmarks/Connectivity/Framework/mains:shiloach_vishkin",
    "//benchmarks/Connectivity/Framework/mains:unite_bfs",
    "//benchmarks/Connectivity/Framework/mains:unite_early_bfs",
    "//benchmarks/Connectivity/Framework/mains:unite_early_kout",
    "//benchmarks/Connectivity/Framework/mains:unite_early_ldd",
    "//benchmarks/Connectivity/Framework/mains:unite_early_nosample",
    "//benchmarks/Connectivity/Framework/mains:unite_kout",
    "//benchmarks/Connectivity/Framework/mains:unite_ldd",
    "//benchmarks/Connectivity/Framework/mains:unite_nd_bfs",
    "//benchmarks/Connectivity/Framework/mains:unite_nd_kout",
    "//benchmarks/Connectivity/Framework/mains:unite_nd_ldd",
    "//benchmarks/Connectivity/Framework/mains:unite_nd_nosample",
    "//benchmarks/Connectivity/Framework/mains:unite_nosample",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_cas_bfs",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_cas_kout",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_cas_ldd",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_cas_nosample",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_lock_bfs",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_lock_kout",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_lock_ldd",
    "//benchmarks/Connectivity/Framework/mains:unite_rem_lock_nosample",
    "//benchmarks/Connectivity/Incremental/mains:jayanti_starting",
    "//benchmarks/Connectivity/Incremental/mains:liutarjan_starting",
    "//benchmarks/Connectivity/Incremental/mains:shiloachvishkin_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_early_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_nd_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_rem_cas_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_rem_lock_starting",
    "//benchmarks/Connectivity/Incremental/mains:unite_starting",
    "//benchmarks/Connectivity/LabelPropagation:Connectivity_main",
    "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity_main",
    "//benchmarks/CycleCounting/Kowalik5Cycle:FiveCycle_main",
    "//benchmarks/DegeneracyOrder/GoodrichPszona11:DegeneracyOrder_main",
    "//benchmarks/GraphColoring/Hasenplaugh14:GraphColoring_main",
    "//benchmarks/KCore/JulienneDBS17:KCore_main",
    "//benchmarks/KTruss:KTruss_main",
    "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition_main",
    "//benchmarks/MaximalIndependentSet/RandomGreedy:MaximalIndependentSet_main",
    "//benchmarks/MaximalIndependentSet/Yoshida:MaximalIndependentSet_main",
    "//benchmarks/MaximalMatching/RandomGreedy:MaximalMatching_main",
    "//benchmarks/MaximalMatching/Yoshida:MaximalMatching_main",
    "//benchmarks/PageRank:PageRank_main",
    "//benchmarks/SSBetweenessCentrality/Brandes:SSBetweennessCentrality_main",
    "//benchmarks/Spanner/MPXV15:Spanner_main",
    "//benchmarks/SpanningForest/BFSSF:SpanningForest_main",
    "//benchmarks/SpanningForest/LabelPropagation:SpanningForest_main",
    "//benchmarks/SpanningForest/SDB14:SpanningForest_main",
    "//benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16:StronglyConnectedComponents_main",
    "//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle_main",
]

# The script will invoke these benchmarks on a weighted graph.
WEIGHTED_GRAPH_BENCHMARKS = [
    "//benchmarks/GeneralWeightSSSP/BellmanFord:BellmanFord_main",
    "//benchmarks/IntegralWeightSSSP/JulienneDBS17:wBFS_main",
    "//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest_main",
    "//benchmarks/MinimumSpanningForest/PBBSMST:MinimumSpanningForest_main",
    "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping_main",
    "//benchmarks/SSWidestPath/JulienneDBS17:SSWidestPath_main",
]

# The script will not invoke these binaries. Shell-style globbing is allowed in
# this list.
IGNORED_BINARIES = [
    # These files take input in a format that's different than that of other
    # benchmarks.
    "//benchmarks/Connectivity/Incremental/mains:*_no_starting",
]


def get_all_benchmark_binaries() -> List[str]:
    """Returns a list of all binaries under `benchmarks/`."""
    return subprocess.run(
        ["bazel", "query", "kind(cc_binary, //benchmarks/...)"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.splitlines()


def check_listed_binaries(
    valid_binaries: List[str], ignored_binaries: List[str]
) -> None:
    """Checks listed binaries for consistency.

    Args:
        valid_binaries: Names of binaries that are considered valid.
        ignored_binaries: Names of binaries that are considered invalid
            and should not be run.

    Raises:
        ValueError: `valid_binaries` and `ignored_binaries` overlap.
        ValueError: `valid_binaries` and `ignored_binaries` combined don't equal
            the set of all benchmark-related binaries as determined by
            `get_all_benchmark_binaries()`.
    """
    remaining_binaries = set(get_all_benchmark_binaries())
    for ignored_binaries_pattern in ignored_binaries:
        conflicting_binaries = fnmatch.filter(valid_binaries, ignored_binaries_pattern)
        if conflicting_binaries:
            raise ValueError(
                "Benchmarks listed as both valid and ignored: {}".format(
                    conflicting_binaries
                )
            )
        binaries_to_ignore = fnmatch.filter(
            remaining_binaries, ignored_binaries_pattern
        )
        if not binaries_to_ignore:
            print(
                "Warning: ignore rule {} has no effect".format(ignored_binaries_pattern)
            )
        remaining_binaries -= set(binaries_to_ignore)

    valid_binaries = set(valid_binaries)
    if remaining_binaries != valid_binaries:
        extra_listed_binaries = valid_binaries - remaining_binaries
        if extra_listed_binaries:
            raise ValueError(
                "Listed benchmarks do not exist: {}".format(extra_listed_binaries)
            )
        missing_listed_binaries = remaining_binaries - valid_binaries
        if missing_listed_binaries:
            raise ValueError(
                "Please update {} to include binaries {}".format(
                    __file__, missing_listed_binaries
                )
            )


def run_all_benchmarks(
    unweighted_graph_benchmarks: List[str],
    weighted_graph_benchmarks: List[str],
    unweighted_graph_file: str,
    weighted_graph_file: str,
    timeout: Optional[int],
) -> List[str]:
    """Runs all benchmarks, returning a list of failing benchmarks.

    Args:
        unweighted_graph_benchmarks: List of all benchmarks to run on the
            unweighted graph.
        weighted_graph_benchmarks: List of all benchmarks to run on the
            weighted graph.
        unweighted_graph_file: File path to the unweighted graph.
        weighted_graph_file: File path to the weighted graph.
        timeout: Benchmarks that run longer than this timeout period in seconds
            are considered to have failed. If this is `None` then the benchmarks
            have no time limit.

    Returns:
        A list of names of benchmarks that exit with a non-zero error code.
    """

    def compile_benchmark(benchmark: str):
        return subprocess.run(
            ["bazel", "build", "--compilation_mode", "opt", benchmark]
        )

    failed_benchmarks = []
    for benchmark in unweighted_graph_benchmarks:
        benchmark_compile = compile_benchmark(benchmark)
        if benchmark_compile.returncode:
            failed_benchmarks.append(benchmark)
            continue
        benchmark_run = subprocess.run(
            [
                "bazel",
                "run",
                "--compilation_mode",
                "opt",
                benchmark,
                "--",
                "-s",
                "-rounds",
                "1",
                unweighted_graph_file,
            ],
            timeout=timeout,
        )
        if benchmark_run.returncode:
            failed_benchmarks.append(benchmark)
    for benchmark in weighted_graph_benchmarks:
        benchmark_compile = compile_benchmark(benchmark)
        if benchmark_compile.returncode:
            failed_benchmarks.append(benchmark)
            continue
        subprocess.run(
            [
                "bazel",
                "run",
                "--compilation_mode",
                "opt",
                benchmark,
                "--",
                "-s",
                "-w",
                "-rounds",
                "1",
                weighted_graph_file,
            ],
            timeout=timeout,
        )
        if benchmark_run.returncode:
            failed_benchmarks.append(benchmark)
    return failed_benchmarks


if __name__ == "__main__":
    check_listed_binaries(
        valid_binaries=UNWEIGHTED_GRAPH_BENCHMARKS + WEIGHTED_GRAPH_BENCHMARKS,
        ignored_binaries=IGNORED_BINARIES,
    )

    parser = argparse.ArgumentParser(
        description=(
            "Runs all benchmarks on the specified input graphs to check "
            "whether the benchmarks run and exit without errors."
        )
    )
    parser.add_argument(
        "--unweighted_graph",
        type=str,
        required=True,
        help=(
            "Absolute path to an unweighted graph on which to run all "
            "unweighted graph benchmarks"
        ),
    )
    parser.add_argument(
        "--weighted_graph",
        type=str,
        required=True,
        help=(
            "Absolute path to a weighted graph on which to run all "
            "weighted graph benchmarks"
        ),
    )
    parser.add_argument(
        "--timeout",
        type=float,
        help=(
            "(seconds) - Fail benchmarks that run longer than this timeout " "period"
        ),
    )
    parsed_args = parser.parse_args()

    failed_benchmarks = run_all_benchmarks(
        unweighted_graph_benchmarks=UNWEIGHTED_GRAPH_BENCHMARKS,
        weighted_graph_benchmarks=WEIGHTED_GRAPH_BENCHMARKS,
        unweighted_graph_file=os.path.expanduser(parsed_args.unweighted_graph),
        weighted_graph_file=os.path.expanduser(parsed_args.weighted_graph),
        timeout=parsed_args.timeout,
    )
    if failed_benchmarks:
        print("Benchmarks failed: {}".format(failed_benchmarks))
        sys.exit(1)
    else:
        print("Success! All benchmarks completed without an error.")
        sys.exit(0)
