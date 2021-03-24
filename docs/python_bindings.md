---
id: python_bindings
title: Python Bindings
---

Next, we illustrate how to use GBBS from Python using the python
bindings we provide (through [pybind11](https://pybind11.readthedocs.io/en/stable/))

Extending the Python bindings after implementing a new benchmark
requires only a few lines of code to add an extra method to the graph
object exported by the library.
Note that the contents of this file have been tested using Python
3.7.6. You may run into issues using an earlier version.

We first build the bindings using Bazel and add the compiled libraries
to the Python path:
``` sh
$ bazel build //pybindings:gbbs_lib.so
$ export PYTHONPATH=$(pwd)/bazel-bin/pybindings/:$PYTHONPATH
```

## Loading SNAP Graphs

Launch the Python REPL, import the library, and import a downloaded
graph from the SNAP dataset.

``` python
>>> import gbbs
>>> G = gbbs.loadSnap("com-youtube.ungraph.txt", undirected=True)
```

This command creates an uncompressed graph in the GBBS format at the
same location as the input (compression can optionally be enabled
using a separate flag).


## Loading Graphs in the Adjacency Array Format

We also provide commands to load graphs that have already been created
in some valid format.

To load a graph in text format:

``` python
>>> import gbbs
>>> G = gbbs.loadGraph("rMatGraph_J_5_100", undirected=True, compressed=False, binary=True)
```

To load a graph in compressed format:

``` python
>>> import gbbs
>>> G = gbbs.loadGraph("rMatGraph_J_5_100.bytepda", undirected=True, compressed=True, binary=False)
```

To load a graph in the binary format:

``` python
>>> import gbbs
>>> G = gbbs.loadGraph("rMatGraph_J_5_100.binary", undirected=True, compressed=False, binary=True)
```

## Running Graph Algorithms

We can then apply the methods defined on graphs as follows:
```
>>> sims = G.PageRank()
```

Other primitives can be applied similarly. For example:
```
>>> components = G.Connectivity()
>>> print(components[10] == components[82])
True
>>> cores = G.KCore() # Computes coreness values
>>> print(cores[10], cores[82])
(41, 50)
```

