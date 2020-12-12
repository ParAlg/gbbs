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

This script should be invoked directly via Python >=3.7. Because this script
calls other Bazel commands, invoking it with `bazel run` won't work.
"""
from typing import List, Optional, Set, Tuple
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
    "//benchmarks/CoSimRank:CoSimRank_main",
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
    "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity_main",
    "//benchmarks/Connectivity/WorkEfficientSDB14:Connectivity_main",
    "//benchmarks/CycleCounting/Kowalik5Cycle:FiveCycle_main",
    "//benchmarks/DegeneracyOrder/BarenboimElkin08:DegeneracyOrder_main",
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
    "//benchmarks/SCAN/IndexBased:SCAN_main",
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
    "//benchmarks/MinimumSpanningForest/Kruskal:MinimumSpanningForest_main",
    "//benchmarks/MinimumSpanningForest/PBBSMST:MinimumSpanningForest_main",
    "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping_main",
    "//benchmarks/SSWidestPath/JulienneDBS17:SSWidestPath_main",
]

# The script will not invoke these binaries. Shell-style globbing is allowed in
# this list.
IGNORED_BINARIES = [
    # These incremental connectivity binaries take input in a format that's
    # different than that of other benchmarks.
    "//benchmarks/Connectivity/Incremental/mains:*_no_starting",
    "//benchmarks/SCAN/IndexBased/experiments:*",
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
    unweighted_graph_file: Optional[str],
    weighted_graph_file: Optional[str],
    are_graphs_compressed: bool,
    timeout: Optional[int],
) -> List[Tuple[str, str]]:
    """Runs all benchmarks, returning a list of failing benchmarks.

    Args:
        unweighted_graph_benchmarks: List of all benchmarks to run on the
            unweighted graph.
        weighted_graph_benchmarks: List of all benchmarks to run on the
            weighted graph.
        unweighted_graph_file: File path to the unweighted graph.
        weighted_graph_file: File path to the weighted graph.
        are_compressed_compressed: Whether the graph files hold compressed
            graphs.
        timeout: Benchmarks that run longer than this timeout period in seconds
            are considered to have failed. If this is `None` then the benchmarks
            have no time limit.

    Returns:
        A list of names of benchmarks that fail along with a failure reason.
    """

    BAZEL_FLAGS = ["--compilation_mode", "opt"]
    gbbs_flags = ["-s", "-rounds", "1"]
    if are_graphs_compressed:
        gbbs_flags += ["-c"]

    benchmarks = []
    if unweighted_graph_file:
        benchmarks += unweighted_graph_benchmarks
    if weighted_graph_file:
        benchmarks += weighted_graph_benchmarks
    # Compile all the benchmarks up front --- it's faster than compiling
    # them individually since Bazel can compile several files in parallel.
    subprocess.run(["bazel", "build"] + BAZEL_FLAGS + ["--keep_going"] + benchmarks)

    failed_benchmarks = []

    def test_benchmark(
        benchmark: str, graph_file: str, additional_gbbs_flags: List[str]
    ) -> None:
        try:
            benchmark_run = subprocess.run(
                ["bazel", "run"]
                + BAZEL_FLAGS
                + [benchmark, "--"]
                + gbbs_flags
                + additional_gbbs_flags
                + [graph_file],
                timeout=timeout,
            )
            if benchmark_run.returncode:
                failed_benchmarks.append(
                    (
                        benchmark,
                        "Exited with error code {}".format(benchmark_run.returncode),
                    )
                )
        except subprocess.TimeoutExpired:
            failed_benchmarks.append((benchmark, "Timeout"))

    if unweighted_graph_file:
        for benchmark in unweighted_graph_benchmarks:
            test_benchmark(
                benchmark=benchmark,
                graph_file=unweighted_graph_file,
                additional_gbbs_flags=[],
            )
    if weighted_graph_file:
        for benchmark in weighted_graph_benchmarks:
            test_benchmark(
                benchmark=benchmark,
                graph_file=weighted_graph_file,
                additional_gbbs_flags=["-w"],
            )

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
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--unweighted_graph",
        "-u",
        type=str,
        help=(
            "Absolute path to an unweighted graph on which to run all "
            "unweighted graph benchmarks. If not provided, those benchmarks "
            "will not be run."
        ),
    )
    parser.add_argument(
        "--weighted_graph",
        "-w",
        type=str,
        help=(
            "Absolute path to a weighted graph on which to run all "
            "weighted graph benchmarks. If not provided, those benchmarks will "
            "not be run."
        ),
    )
    parser.add_argument(
        "--compressed",
        "-c",
        action="store_true",
        help="Add this flag if input graphs are compressed graphs.",
    )
    parser.add_argument(
        "--timeout",
        "-t",
        type=float,
        default=60,
        help="(seconds) - Halt benchmarks that run longer than this time.",
    )
    parsed_args = parser.parse_args()
    if not parsed_args.unweighted_graph and not parsed_args.weighted_graph:
        parser.error(
            "At least one of --unweighted_graph and --weighted_graph is required."
        )

    unweighted_graph_file = (
        os.path.abspath(parsed_args.unweighted_graph)
        if parsed_args.unweighted_graph
        else None
    )
    weighted_graph_file = (
        os.path.abspath(parsed_args.weighted_graph)
        if parsed_args.weighted_graph
        else None
    )

    failed_benchmarks = run_all_benchmarks(
        unweighted_graph_benchmarks=UNWEIGHTED_GRAPH_BENCHMARKS,
        weighted_graph_benchmarks=WEIGHTED_GRAPH_BENCHMARKS,
        unweighted_graph_file=unweighted_graph_file,
        weighted_graph_file=weighted_graph_file,
        are_graphs_compressed=parsed_args.compressed,
        timeout=parsed_args.timeout,
    )
    if failed_benchmarks:
        print("Benchmarks failed: {}".format(failed_benchmarks))
        sys.exit(1)
    else:
        print("Success! All benchmarks completed without an error.")
        sys.exit(0)
