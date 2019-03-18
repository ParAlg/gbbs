---
title: "Our Results"
layout: table_archive
permalink: /results/
comments: false
mathjax: true
---

{% include mathjax.html %}

This page presents our results for the GBBS benchmark on a diverse set of
real-world graphs.

Our focus is primarily on the Hyperlink2012 graph, which has 224 billion
(symmetric) edges. This graph is currently (as of March, 2019) the largest
publicly available real-world graph.

We also include our results on the Hyperlink2014 graph, which is roughly half
the size of Hyperlink2012 and the Clueweb graph, which is about one third the
size of Hyperlink2012.

# Hyperlink2012
{% include hyperlink2012_ourresults.html %}

<small>
<b>T1</b> is the 1-thread running time in seconds. <b>T72 (h)</b> is the 72-thread
running time in seconds with 2-way hyperthreading enabled.
</small>

# Hyperlink2014
{% include hyperlink2014_ourresults.html %}

<small>
<b>T1</b> is the 1-thread running time in seconds. <b>T72 (h)</b> is the 72-thread
running time in seconds with 2-way hyperthreading enabled.
</small>

# Clueweb
{% include clueweb_ourresults.html %}

<small>
<b>T1</b> is the 1-thread running time in seconds. <b>T72 (h)</b> is the 72-thread
running time in seconds with 2-way hyperthreading enabled.
</small>

# Experimental Setup
We run all of our experiments on a 72-core Dell PowerEdge R930 (with two-way
hyper-threading) with $$4\times 2.4$$GHz Intel 18-core E7-8867 v4 Xeon
processors (with a 4800MHz bus and 45MB L3 cache) and 1TBof main memory.
Our programs use Cilk Plus to express parallelism and are compiled with the
`g++` compiler (version 5.4.1) with the `-O3` flag.

By using Cilk's work-stealing scheduler we are able obtain an expected running
time of $$W/P + O(D)$$ for an algorithm with $$W$$ work and $$D$$ depth on $$P$$
processors [\[Blumofe1999\]](http://supertech.csail.mit.edu/papers/steal.pdf).

For the parallel experiments, we use the command `numactl -i all` to balance the
memory allocations across the sockets.  All of the speedup numbers we report are
the running times of our parallel implementation on 72-cores with
hyper-threading over the running time of the implementation on a single thread.

