---
title: "Single-Source Betweenness Centrality"
layout: benchmark
permalink: /benchmark/bc
comments: false
mathjax: true
sidebar:
  nav: "benchmark"
---

## Input
$$G=(V, E)$$, an unweighted graph in the [adjacency graph
format](/benchmark/formats/), and a source, $$src \in V$$.

## Output
$$S$$, a [mapping](/benchmark/definitions/) from each vertex $$v$$ to the
centrality contribution from all $$(src, t)$$ shortest paths that pass through
$$v$$.

## Compilation
```
$ make BC
g++ -I../ -I../src -mcx16 -ldl -std=c++17 -march=native -Wall -DLONG  -O3 -DAMORTIZEDPD -fconcepts -DCONCEPTS -L/usr/local/lib -ljemalloc -DHOMEGROWN -pthread -o BC BC.C
```

## Running BC

### soc-LiveJournal, src = 10012
```
$ numactl -i all ./BC -rounds 10 -src 10012 -s -m inputs/soc-LiveJournal1_sym.adj
...
```

### com-Orkut, src = 10012
```
$ numactl -i all ./BC -rounds 10 -src 10012 -s -m inputs/com-orkut.adj
...
```

### Twitter, src = 10012
```
$ numactl -i all ./BC -rounds 10 -src 10012 -s -m inputs/twitter_SJ.adj
...
```

### ClueWeb, src = 100000
```
$ numactl -i all ./BC -rounds 10 -src 100000 -s -c -m inputs/clueweb_sym.bytepda
...
```

### Hyperlink2014, src = 100000
```
$ numactl -i all ./BC -rounds 10 -src 100000 -s -c -m inputs/hyperlink2014.bytepda
...
```

### Hyperlink2012, src = 100000
```
$ numactl -i all ./BC -rounds 10 -src 100000 -s -c -m inputs/hyperlink2012.bytepda
...
```
