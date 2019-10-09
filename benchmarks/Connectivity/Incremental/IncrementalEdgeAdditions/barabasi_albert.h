#pragma once


/* Sequential barbasi_albert generator.
 * Outputs a sequence of *directed* updates; m per new node. */
template <class Graph>
auto generate_barabasi_albert(
    Graph& base_graph,
    uint64_t nodes_to_add,
    uint64_t edges_per_node,
    bool symmetrize=false) {
  /* Sample according to degree distribution */

  size_t base_n = base_graph.n;
  size_t n = base_n + nodes_to_add;
  /* Initialize initial degree disribution. */
  auto degree_distribution = pbbs::sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) {
    if (i < base_n) {
      degree_distribution[i] = base_graph.get_vertex(i).getOutDegree();
    } else {
      degree_distribution[i] = 0;
    }
  });

  /* Ask : how to generate according to dist? */


  /* Output: directed edge updates, or undirected edge updates? Different
   * algorithms may require diff. things. Use _symmetrize_ flag to indicate
   * whether we should double and remove duplicates. */
}

/* As in Kanat's paper, show a figure showing for each point in the stream the
 * number of connected components in the stream */
