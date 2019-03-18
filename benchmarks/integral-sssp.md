---
title: "Integral-Weight SSSP"
layout: benchmark
permalink: /benchmarks/integral-sssp/
comments: false
mathjax: true
sidebar:
  nav: "benchmarks"
---

## Input
\\(G=(V, E, w)\\), a weighted graph in the [adjacency graph
format](/benchmarks/formats/), with integral edge weights, \\(src \in
V\\).

## Output
\\(D\\), a [mapping](/benchmarks/definitions/) where \\(D[v]\\) is the
shortest path distance from \\(src\\) to \\(v\\) in \\(G\\) and \\(\infty\\) if
\\(v\\) is unreachable.

## Compilation
```
$ make wBFS
g++ -I../ -I../src -mcx16 -ldl -std=c++17 -march=native -Wall -DLONG  -O3 -DAMORTIZEDPD -fconcepts -DCONCEPTS -L/usr/local/lib -ljemalloc -DHOMEGROWN -pthread -o wBFS wBFS.C
```

## Running wBFS

### soc-LiveJournal, src = 10012
```
$ numactl -i all ./wBFS -rounds 10 -src 10012 -s -m inputs/soc-LiveJournal1_sym.adj
...
```

### com-Orkut, src = 10012
```
$ numactl -i all ./wBFS -rounds 10 -src 10012 -s -m inputs/com-orkut.adj
...
```

### Twitter, src = 10012
```
$ numactl -i all ./wBFS -rounds 10 -src 10012 -s -m inputs/twitter_SJ.adj
...
```

### ClueWeb, src = 100000
```
$ numactl -i all ./wBFS -rounds 10 -src 100000 -s -c -m inputs/clueweb_sym.bytepda
...
```

### Hyperlink2014, src = 100000
```
$ numactl -i all ./wBFS -rounds 10 -src 100000 -s -c -m inputs/hyperlink2014.bytepda
...
```

### Hyperlink2012, src = 100000
```
$ numactl -i all ./wBFS -rounds 10 -src 100000 -s -c -m inputs/hyperlink2012.bytepda
...
```
