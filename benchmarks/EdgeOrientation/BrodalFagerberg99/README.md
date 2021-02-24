Some notation and data structure description, basically copied from
Henzinger et al. (2020). (https://arxiv.org/pdf/2002.10142.pdf)

notation:
levels : log^2(n) many, notated as l. level of a vertex is l(v)
groups : log(n) many. group of a particular level if g(l)

invariants:
* (1) each v \in V has at most 5 * 2^{g(l(v))} neighbors at its own, or higher
    levels. (not too many neighbors above)
* (2) each v \in V where l(v) > 1 has at least 2^{g(l(v) - 1)} neighbors in
    level >= l(v) - 1. (enough neighbors in the level below and above)

vertices that do not satisfy these invariants are said to be **dirty**.

The Henzinger et al. algorithm works by moving vertices that are dirty to
a level where both invariants are satisfied.

data structures maintained:
For each vertex v, and each level l' < l(v), maintain the following:
- doubly-linked list neighbors(v, l') containing all neighbors of V in
  level l'

Furthermore, we have another doubly-linked list up_neighbors(v) containing all
neighbors v in levels l' >= l(v).

For each edge (v,u), we store a pointer to the position of u in neighbors or
up_neighbors (and vice versa).

In addition, we store the size of each list (neighbors and up_neighbors) which
lets us test invariants (1) and (2) in O(1) work.
