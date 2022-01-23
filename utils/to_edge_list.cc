#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include "gbbs/gbbs.h"

namespace gbbs {
/* Format:
 * Emits a lexicographically sorted list of edges. The format takes (u,v)
 * and writes:
 * "u v\n" */
template <class Graph>
void print_edge_list(Graph& GA, std::string& outfile, bool direct_sym,
                     bool multistep_header) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto edges = sequence<std::tuple<uintE, uintE>>(m);

  auto offs = sequence<size_t>(n);
  parallel_for(0, n, [&](size_t i) {
    size_t ctr = 0;
    if (direct_sym) {
      auto f = [&](const uintE& u, const uintE& v, const W& wgh) {
        if (u < v) ctr++;
      };
      GA.get_vertex(i).out_neighbors().map(f, false);
    } else {
      ctr = GA.get_vertex(i).out_degree();
    }
    offs[i] = ctr;
  });
  size_t m_out = parlay::scan_inplace(make_slice(offs));

  parallel_for(0, n, [&](size_t i) {
    size_t off = offs[i];
    size_t ctr = 0;
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (direct_sym) {
        if (u < v) {
          edges[off + ctr] = std::make_tuple(u, v);
          ctr++;
        }
      } else {
        edges[off + ctr] = std::make_tuple(u, v);
        ctr++;
      }
    };
    GA.get_vertex(i).out_neighbors().map(map_f, false);
  });

  std::ofstream file(outfile, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << outfile << std::endl;
    exit(0);
  }
  if (multistep_header) {
    file << n << " " << m_out << std::endl;
  }
  auto edges_chars = parlay::sequence_to_string(edges);
  file.write(edges_chars.begin(), edges_chars.size());

  file.close();
  std::cout << "Done" << std::endl;
}

/* Format:
 * Emits a lexicographically sorted list of edges.
 *  n m              # header (once)
 * "u v 1\n"         # per edge*/
template <class Graph>
void print_edge_list_matrixmarket(Graph& GA, std::string& outfile) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto edges = sequence<std::tuple<uintE, uintE>>(m);

  auto offs = sequence<size_t>(n);
  parallel_for(0, n, [&](size_t i) {
    size_t ctr = 0;
    auto f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) ctr++;
    };
    GA.get_vertex(i).out_neighbors().map(f, false);
    offs[i] = ctr;
  });
  size_t m_out = parlay::scan_inplace(make_slice(offs));

  parallel_for(0, n, [&](size_t i) {
    size_t off = offs[i];
    size_t ctr = 0;
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        edges[off + ctr] = std::make_tuple(u, v);
        ctr++;
      }
    };
    GA.get_vertex(i).out_neighbors().map(map_f, false);
  });

  std::ofstream file(outfile, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << outfile << std::endl;
    exit(0);
  }
  file << n << " " << m_out << std::endl;
  auto edges_chars = parlay::sequence_to_string(edges);
  file.write(edges_chars.begin(), edges_chars.size());

  file.close();
  std::cout << "Done" << std::endl;
}

template <class Graph>
double Reorderer(Graph& GA, commandLine P) {
  auto outfile = P.getOptionValue("-of", "");
  auto direct_sym =
      P.getOptionValue("-direct_sym"); /* only emit { (u,v) | u < v } */
  auto multistep_header =
      P.getOptionValue("-multistep_header"); /* only emit { (u,v) | u < v } */
  auto matrixmarket_format = P.getOptionValue(
      "-matrixmarket_format"); /* only emit { (u,v) | u < v } */
  if (matrixmarket_format) {
    print_edge_list_matrixmarket(GA, outfile);
  } else {
    print_edge_list(GA, outfile, direct_sym, multistep_header);
  }
  exit(0);
  return 1.0;
}
}  // namespace gbbs

generate_main(gbbs::Reorderer, false);
