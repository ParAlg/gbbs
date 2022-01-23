#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include "gbbs/gbbs.h"

namespace gbbs {
template <class Graph>
void randomReorder(Graph& GA, std::string& outfile) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;
  auto perm = parlay::random_permutation<uintE>(n);

  auto edges = sequence<uintE>(m);
  auto offs = sequence<uintT>(n);
  parallel_for(
      0, n, [&](size_t i) { offs[perm[i]] = GA.get_vertex(i).out_degree(); });
  size_t tot = parlay::scan_inplace(make_slice(offs));
  std::cout << "m = " << m << " tot = " << tot << std::endl;

  parallel_for(0, n, [&](size_t i) {
    size_t off = offs[perm[i]];
    size_t next_off = (perm[i] == (n - 1)) ? m : offs[perm[i] + 1];
    auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      uintE ngh_perm = perm[ngh];
      edges[off++] = ngh_perm;
    };
    GA.get_vertex(i).out_neighbors().map(map_f, false);
    std::sort(edges.begin() + off, edges.begin() + next_off,
              std::less<uintE>());
  });

  std::ofstream file(outfile, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << outfile << std::endl;
    exit(0);
  }
  file << "AdjacencyGraph" << std::endl;
  file << n << std::endl;
  file << m << std::endl;
  auto off_chars = parlay::sequence_to_string(offs);
  auto edges_chars = parlay::sequence_to_string(edges);
  file.write(off_chars.begin(), off_chars.size());
  file.write(edges_chars.begin(), edges_chars.size());

  file.close();
  std::cout << "Done" << std::endl;
}

template <class Graph>
double Reorderer(Graph& GA, commandLine P) {
  auto outfile =
      P.getOptionValue("-of", "/ssd0/graphs/bench_experiments/out.adj");
  randomReorder(GA, outfile);
  exit(0);
  return 1.0;
}
}  // namespace gbbs

generate_main(gbbs::Reorderer, false);
