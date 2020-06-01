---
id: sage
title: Sage (Semi-Asymmetric Graph Engine)
---

This page contains information about how to build and run the
NVRAM-based code used in our paper [Sage: Parallel Semi-Asymmetric
Graph Algorithms for NVRAMs](https://arxiv.org/abs/1910.12310).

## Download
Our code is currently included as part of the GBBS benchmarks, under the
top-level `sage` directory. You can download the code by cloning the
gbbs repository through ssh (see our
[github](https://github.com/ParAlg/gbbs) for other download options).

```
git clone git@github.com:ParAlg/gbbs.git
```

## Compile
Compilation requires the `bazel` build tool, which you can find
installation instructions for
[here](https://docs.bazel.build/versions/master/install.html).

After installing bazel, all of the Sage implementations can be built
by running:
```
bazel build //sage/...
```

This command will compile binaries and store them in the top-level
`bazel-bin` directory. For example, the BFS benchmark is located at:
```
bazel-bin/sage/benchmarks/BFS/BFS_main
```

## Graph Formats
The code currently requires users to use either the binary
compressed-sparse row format (binary CSR), or the parallel byte-encoded
compressed-sparse row format (parallel-byte CSR). Both formats are supported by GBBS, and
we provide more information about the formats
[here](https://paralg.github.io/gbbs/docs/formats).

The default format is the binary CSR. Specifying that the input graph
is compressed is done through the command-line argument `-c`.

## Running

