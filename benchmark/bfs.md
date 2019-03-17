---
title: "Breadth-First Search"
layout: benchmark
permalink: /benchmark/bfs
comments: false
mathjax: true
sidebar:
  nav: "benchmark"
benchname: "Breadth-First Search"
---

### Input
$$G=(V, E)$$, an unweighted graph in the [adjacency graph
format](/benchmark/formats/), and a source, $$src \in V$$.

### Output
$$D$$, a [mapping](/benchmark/definitions/) where $$D[v]$$ is the
shortest path distance from $$src$$ to $$v \in V$$ and $$\infty$$ if
$$v$$ is unreachable, and $$reach$$, the number of reachable vertices.


#### Todo:
Change to output distance instead of the parent.

### Compilation
```
$ make BFS
g++ -I../ -I../src -mcx16 -ldl -std=c++17 -march=native -Wall -DLONG  -O3 -DAMORTIZEDPD -fconcepts -DCONCEPTS -L/usr/local/lib -ljemalloc -DHOMEGROWN -pthread -o BFS BFS.C
```

### Running BFS

#### soc-LiveJournal, src = 10012
```
$ numactl -i all ./BFS -rounds 10 -src 10012 -s -m inputs/soc-LiveJournal1_sym.adj
...
Reachable: 4843953
PBBS time: Running time: 0.0181
time per iter: 0.01836
```

#### com-Orkut, src = 10012
```
$ numactl -i all ./BFS -rounds 10 -src 10012 -s -m inputs/com-orkut.adj
...
Reachable: 3072441
PBBS time: Running time: 0.0092
time per iter: 0.01049
```

#### Twitter, src = 10012
```
$ numactl -i all ./BFS -rounds 10 -src 10012 -s -m inputs/twitter_SJ.adj
...
Reachable: 41652230
PBBS time: Running time: 0.0805
time per iter: 0.08433
```

#### ClueWeb, src = 100000
```
$ numactl -i all ./BFS -rounds 10 -src 100000 -s -c -m inputs/clueweb_sym.bytepda
...
```

#### Hyperlink2014, src = 100000
```
$ numactl -i all ./BFS -rounds 10 -src 100000 -s -c -m inputs/hyperlink2014.bytepda
...
```

#### Hyperlink2012, src = 100000
```
$ numactl -i all ./BFS -rounds 10 -src 100000 -s -c -m inputs/hyperlink2012.bytepda
...
```
