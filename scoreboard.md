---
title: "Scoreboard"
layout: table_archive
permalink: /scoreboard/
comments: false
mathjax: true
---

{% include mathjax.html %}

To put our results in context, we compare our 72-core running times to running
times reported by existing work in the literature. Our goal is to compare the
performance achieved state-of-the-art graph processing systems (distributed,
external-memory, shared-memory, etc.) for a set of fundamental problems over the
_same graph_.

We note that the times reported here are those reported by the authors in their
papers.

# Hyperlink2012

### BFS

{:.display}
| Time (s) | Hyper-threads | Memory (Gb) | Flash (Gb) | HDD (Gb) | Type | System |
| -------- | -----         | ------      | -----      | --- | ---- | ------ |
| 16.7     | 144           | 350*        | --         | -- | SM | GBBS |
| 208      | 64            | 22*         | 1500       | -- | EM | FlashGraph |
| 2500     | 32            | 64          | 2500       | -- | EM | BigSparse |


<small>
<b>\*</b> indicates that the authors reported the _peak-memory_ used by the system, which may be smaller than the total memory of the machine(s) used.
</small>


### Betweenness Centrality

{:.display}
| Time (s) | Hyper-threads | Memory (Gb) | Flash (Gb) | HDD (Gb) | Type | System |
| ------------- | ----- | ------ | ----- | --- | ---- | ------ |
| 35.2 | 144 | 350* | -- | -- | Shared Mem. | GBBS |
| 208  | 64  | 22*  | 1500 | -- | External Mem. | BigSparse |
| 2500 | 32  | 64  | 2500 | -- | External Mem. | BigSparse |

