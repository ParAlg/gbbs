#pragma once

#include <math.h>

#include <limits>

#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"

#define INDUCED_STACK_THR 5000

namespace gbbs {

struct HybridSpace_lw {
  // Number of induced neighbors, per recursive level
  uintE* num_induced = nullptr;
  // Induced neighbors, per recurisve level
  uintE* induced = nullptr;
  // Number of edges in induced subgraph on induced neighbors, from first
  // recursive level
  uintE num_edges = 0;
  // Edges in induced subgraph on induced neighbors, from first recursive level
  uintE* induced_edges = nullptr;
  // Degrees of each induced neighbor in induced subgrpah, from first recursive
  // level
  uintE* induced_degs = nullptr;
  // Label each induced neighbor with latest recursive level it was included in
  char* labels = nullptr;
  // Label each vertex for fast intersect for first recursive level
  uintE* old_labels = nullptr;
  // Number of vertices
  size_t nn = 0;
  bool use_old_labels = true;  // use label intersect for first recursive level
  bool free_relabel = true;    // if true, take ownership of induced_edges,
                               // induced_degs, and relabel
  uintE* relabel = nullptr;    // if counting per vertex, must undo local
                               // relabeling for 0 to alpha indices
  bool use_base = false;       // if true, counting per vertex
  HybridSpace_lw() {}

  void alloc(size_t max_induced, size_t k, size_t n, bool _use_old_labels,
             bool _use_base, bool _free_relabel = true) {
    use_old_labels = _use_old_labels;
    use_base = _use_base;
    free_relabel = _free_relabel;
    if (induced == nullptr)
      induced = (uintE*)malloc(sizeof(uintE) * k * max_induced);
    if (free_relabel && induced_degs == nullptr)
      induced_degs = (uintE*)malloc(sizeof(uintE) * max_induced);
    if (labels == nullptr) labels = (char*)calloc(max_induced, sizeof(char));
    if (free_relabel && induced_edges == nullptr)
      induced_edges = (uintE*)malloc(sizeof(uintE) * max_induced * max_induced);
    if (num_induced == nullptr) num_induced = (uintE*)malloc(sizeof(uintE) * k);
    if (free_relabel && use_old_labels && old_labels == nullptr)
      old_labels = (uintE*)calloc(n, sizeof(uintE));
    if (free_relabel && use_base && relabel == nullptr)
      relabel = (uintE*)malloc(sizeof(uintE) * max_induced);
  }

  void copy(HybridSpace_lw* space) {
    // Shallow copy
    if (use_base) relabel = space->relabel;
    induced_edges = space->induced_edges;
    induced_degs = space->induced_degs;
    nn = space->nn;
  }

  template <class Graph, class F>
  void setup(Graph& DG, size_t k, size_t i, F f) {
    if (use_old_labels)
      setup_labels(DG, DG, k, i, f);
    else
      setup_intersect(DG, DG, k, i, f);
  }

  template <class Graph>
  void setup(Graph& DG, size_t k, size_t i) {
    auto f = [&](const uintE& src, const uintE& u) { return true; };
    if (use_old_labels)
      setup_labels(DG, DG, k, i, f);
    else
      setup_intersect(DG, DG, k, i, f);
  }

  template <class Graph, class Graph2, class F>
  void setup(Graph& G, Graph2& DG, size_t k, size_t i, F f) {
    if (use_old_labels)
      setup_labels(G, DG, k, i, f);
    else
      setup_intersect(G, DG, k, i, f);
  }

  // Perform first level recursion, using space-efficient intersection
  template <class Graph, class Graph2, class F>
  void setup_intersect(Graph& DG, Graph2& DG2, size_t k, size_t i, F f) {
    using W = typename Graph::weight_type;

    // Set up relabeling if counting per vertex
    if (use_base) {
      size_t j = 0;
      auto map_base_f = [&](const uintE& src, const uintE& v, const W& wgh) {
        relabel[j] = v;
        j++;
      };
      DG.get_vertex(i).out_neighbors().map(map_base_f, false);
    }

    // Set up first level induced neighborhood (neighbors of vertex i, relabeled
    // from 0 to degree of i)
    nn = DG.get_vertex(i).out_degree();
    parallel_for(0, nn, [&](size_t j) { induced_degs[j] = 0; });
    num_induced[0] = nn;
    parallel_for(0, nn, [&](size_t j) { induced[j] = j; });

    size_t j = 0;
    auto map_f = [&](const uintE& src, const uintE& v, const W& wgh) {
      // Return if neighbor v is invalid
      if (!f(src, v)) {
        j++;
        return;
      }
      // For a neighbor v of i, intersect N(v) with N(i)
      // These are the edges adjacent to v in the induced neighborhood of i
      // Note that v is relabeled to j under the relabeling from 0 to degree of
      // i
      // Store these edges in induced_edges[j*nn]
      // Store the number of edge (degree) in induced_degs[j]
      size_t v_deg = DG2.get_vertex(v).out_degree();
      auto i_iter = DG.get_vertex(i).out_neighbors().get_iter();
      auto v_iter = DG2.get_vertex(v).out_neighbors().get_iter();
      size_t i_iter_idx = 0;
      size_t v_iter_idx = 0;

      while (i_iter_idx < nn && v_iter_idx < v_deg) {
        // Check if the neighbors of v and i match
        if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
          // Check if the edges (i, i_iter.cur()) and (v, v_iter.cur()) are
          // valid
          if (f(i, std::get<0>(i_iter.cur())) &&
              f(v, std::get<0>(i_iter.cur()))) {
            // Save the neighbor and increment the induced degree on j
            induced_edges[j * nn + induced_degs[j]] = i_iter_idx;
            induced_degs[j]++;
          }
          // Move to next neighbors
          i_iter_idx++;
          v_iter_idx++;
          if (i_iter.has_next()) i_iter.next();
          if (v_iter.has_next()) v_iter.next();
        } else if (std::get<0>(i_iter.cur()) < std::get<0>(v_iter.cur())) {
          i_iter_idx++;
          if (i_iter.has_next()) i_iter.next();
        } else {
          v_iter_idx++;
          if (v_iter.has_next()) v_iter.next();
        }
      }
      j++;
    };
    DG.get_vertex(i).out_neighbors().map(map_f, false);

    // Count total number of edges in induced neighborhood
    auto deg_seq = gbbs::make_slice(induced_degs, nn);
    num_edges = parlay::reduce(deg_seq);
  }

  // Perform first level recursion, using linear space to intersect
  template <class Graph, class Graph2, class F>
  void setup_labels(Graph& DG, Graph2& DG2, size_t k, size_t i, F f) {
    using W = typename Graph::weight_type;

    // Set up first level induced neighborhood (neighbors of vertex i, relabeled
    // from 0 to degree of i)
    nn = DG.get_vertex(i).out_degree();
    parallel_for(0, nn, [&](size_t j) { induced_degs[j] = 0; });
    num_induced[0] = nn;
    parallel_for(0, nn, [&](size_t j) { induced[j] = j; });

    size_t o = 0;
    auto map_label_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      // Return if edge is invalid
      if (!f(src, ngh)) {
        o++;
        return;
      }
      // Set up label for intersection
      old_labels[ngh] = o + 1;
      // Set up relabeling if counting per vertex
      if (use_base) {
        relabel[o] = ngh;
      }
      o++;
    };
    DG.get_vertex(i).out_neighbors().map(map_label_f, false);

    size_t j = 0;
    auto map_f = [&](const uintE& src, const uintE& v, const W& wgh) {
      // Return if neighbor v is invalid
      if (!f(src, v)) {
        j++;
        return;
      }
      // For a neighbor v of i, intersect N(v) with N(i)
      // These are the edges adjacent to v in the induced neighborhood of i
      // Note that v is relabeled to j under the relabeling from 0 to degree of
      // i
      // Store these edges in induced_edges[j*nn]
      // Store the number of edge (degree) in induced_degs[j]
      auto map_nbhrs_f = [&](const uintE& src_v, const uintE& v_nbhr,
                             const W& wgh_v) {
        // Check if the two-hop edge (v, v_nbhr) is valid
        if (!f(src_v, v_nbhr)) return;
        // Check if the neighbor of v is also a valid neighbor of i
        if (old_labels[v_nbhr] > 0) {
          // Save the neighbor and increment the induced degree on j
          induced_edges[j * nn + induced_degs[j]] = old_labels[v_nbhr] - 1;
          induced_degs[j]++;
        }
      };
      DG2.get_vertex(v).out_neighbors().map(map_nbhrs_f, false);
      j++;
    };
    DG.get_vertex(i).out_neighbors().map(map_f, false);

    // Reset the array used for intersecting
    auto map_relabel_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = 0;
    };
    DG.get_vertex(i).out_neighbors().map(map_relabel_f, false);

    // Count total number of edges in induced neighborhood
    auto deg_seq = gbbs::make_slice(induced_degs, nn);
    num_edges = parlay::reduce(deg_seq);
  }

  template <class Graph, class F>
  void setup_edge(Graph& DG, size_t k, size_t i, size_t l, F f) {
    if (use_old_labels)
      setup_labels_edge(DG, k, i, l, f);
    else
      setup_intersect_edge(DG, k, i, l, f);
  }

  template <class Graph>
  void setup_edge(Graph& DG, size_t k, size_t i, size_t l) {
    auto f = [&](const uintE& src, const uintE& u) { return true; };
    if (use_old_labels)
      setup_labels_edge(DG, k, i, l, f);
    else
      setup_intersect_edge(DG, k, i, l, f);
  }

  // Perform first and second level recursion, using linear space to intersect
  template <class Graph, class F>
  void setup_labels_edge(Graph& DG, size_t k, size_t i, size_t l, F f) {
    using W = typename Graph::weight_type;

    // Mark induced neighbors in the first recursive level
    auto map_label_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      // Return if edge is invalid
      if (!f(src, ngh)) return;
      old_labels[ngh] = UINT_E_MAX;
    };
    DG.get_vertex(i).out_neighbors().map(map_label_f, false);

    // Label induced neighbors in the second recursive level (considering
    // neighbors of i and l)
    size_t o = 0;
    auto lmap_label_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      // Return if edge is invalid
      if (!f(src, ngh)) {
        return;
      }
      if (old_labels[ngh] == UINT_E_MAX) {
        // Set up label for intersection
        induced[o] = ngh;
        old_labels[ngh] = o + 1;
        // Set up relabeling if counting per vertex
        if (use_base) {
          relabel[o] = ngh;
        }
        o++;
      }
    };
    DG.get_vertex(l).out_neighbors().map(lmap_label_f, false);

    // Reset induced neighbors from the first recursive level but not the second
    auto remap_label_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      if (old_labels[ngh] == UINT_E_MAX) old_labels[ngh] = 0;
    };
    DG.get_vertex(i).out_neighbors().map(remap_label_f, true);

    // Set up second level induced neighborhood (neighbors of vertex i and
    // vertex l, relabeled)
    nn = o;
    num_induced[0] = o;
    parallel_for(0, nn, [&](size_t p) { induced_degs[p] = 0; });

    for (size_t j = 0; j < nn; j++) {
      auto v = induced[j];
      // For a neighbor v in the second level induced neighborhood, intersect
      // N(v) with the second level induced neighborhood
      // Note that v is relabeled to j under the relabeling
      // Store these edges in induced_edges[j*nn]
      // Store the number of edge (degree) in induced_degs[j]
      auto map_nbhrs_f = [&](const uintE& src_v, const uintE& v_nbhr,
                             const W& wgh_v) {
        // Return if neighbor v_nbhr is invalid
        if (!f(src_v, v_nbhr)) return;
        if (old_labels[v_nbhr] > 0) {
          induced_edges[j * nn + induced_degs[j]] = old_labels[v_nbhr] - 1;
          induced_degs[j]++;
        }
      };
      DG.get_vertex(v).out_neighbors().map(map_nbhrs_f, false);
    }

    // Set up vertices in the second level induced neighborhood
    for (size_t p = 0; p < nn; p++) {
      induced[p] = p;
    }

    // Reset the array used for intersecting
    auto lremap_label_f = [&](const uintE& src, const uintE& ngh,
                              const W& wgh) { old_labels[ngh] = 0; };
    DG.get_vertex(l).out_neighbors().map(lremap_label_f, true);

    // Count total number of edges in induced neighborhood
    auto deg_seq = gbbs::make_slice(induced_degs, nn);
    num_edges = parlay::reduce(deg_seq);
  }

  // Perform first and second level recursion, using space-efficient
  // intersection
  template <class Graph, class F>
  void setup_intersect_edge(Graph& DG, size_t k, size_t i, size_t l, F f) {
    // Retrieve induced neighbors in intersection of N(i) and N(l)
    size_t j = 0;
    size_t v_deg = DG.get_vertex(l).out_degree();
    size_t i_deg = DG.get_vertex(i).out_degree();
    auto i_iter = DG.get_vertex(i).out_neighbors().get_iter();
    auto v_iter = DG.get_vertex(l).out_neighbors().get_iter();
    size_t i_iter_idx = 0;
    size_t v_iter_idx = 0;
    while (i_iter_idx < i_deg && v_iter_idx < v_deg) {
      if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
        if (f(i, std::get<0>(i_iter.cur())) &&
            f(l, std::get<0>(i_iter.cur()))) {
          induced[j] = i_iter_idx;
          j++;
        }
        i_iter_idx++;
        v_iter_idx++;
        if (i_iter.has_next()) i_iter.next();
        if (v_iter.has_next()) v_iter.next();
      } else if (std::get<0>(i_iter.cur()) < std::get<0>(v_iter.cur())) {
        i_iter_idx++;
        if (i_iter.has_next()) i_iter.next();
      } else {
        v_iter_idx++;
        if (v_iter.has_next()) v_iter.next();
      }
    }

    // Set up second level induced neighborhood
    nn = j;
    num_induced[0] = j;
    // Set up relabeling if counting per vertex
    if (use_base) {
      parallel_for(0, nn, [&](size_t p) { relabel[p] = induced[p]; });
    }
    parallel_for(0, nn, [&](size_t p) { induced_degs[p] = 0; });

    for (size_t p = 0; p < nn; p++) {
      // For a neighbor v = induced[p] in the second level induced neighborhood,
      // intersect N(v) with the second level induced neighborhood
      // Note that v is relabeled to p under the relabeling
      // Store these edges in induced_edges[p*nn]
      // Store the number of edge (degree) in induced_degs[p]
      size_t u_deg = DG.get_vertex(induced[p]).out_degree();
      auto u_iter = DG.get_vertex(induced[p]).out_neighbors().get_iter();
      size_t u_iter_idx = 0;
      i_iter_idx = 0;
      while (i_iter_idx < nn && u_iter_idx < u_deg) {
        if (induced[i_iter_idx] == std::get<0>(u_iter.cur())) {
          if (f(induced[p], induced[i_iter_idx])) {
            induced_edges[p * nn + induced_degs[p]] = i_iter_idx;
            induced_degs[p]++;
          }
          i_iter_idx++;
          u_iter_idx++;
          if (u_iter.has_next()) u_iter.next();
        } else if (induced[i_iter_idx] < std::get<0>(u_iter.cur())) {
          i_iter_idx++;
        } else {
          u_iter_idx++;
          if (u_iter.has_next()) u_iter.next();
        }
      }
    }

    // Set up vertices in the second level induced neighborhood
    for (size_t p = 0; p < nn; j++) {
      induced[p] = p;
    }

    // Count total number of edges in induced neighborhood
    auto deg_seq = gbbs::make_slice(induced_degs, nn);
    num_edges = parlay::reduce(deg_seq);
  }

  static void init() {}
  static void finish() {}

  ~HybridSpace_lw() {
    if (labels) {
      free(labels);
      labels = nullptr;
    }
    if (induced) {
      free(induced);
      induced = nullptr;
    }
    if (free_relabel && induced_edges) {
      free(induced_edges);
      induced_edges = nullptr;
    }
    if (free_relabel && induced_degs) {
      free(induced_degs);
      induced_degs = nullptr;
    }
    if (num_induced) {
      free(num_induced);
      num_induced = nullptr;
    }
    if (use_old_labels && old_labels) {
      free(old_labels);
      old_labels = nullptr;
    }
    if (use_base && relabel && free_relabel) {
      free(relabel);
      relabel = nullptr;
    }
    if (!free_relabel) {
      relabel = nullptr;
      induced_edges = nullptr;
      induced_degs = nullptr;
    }
  }
};

// The other space types are deprecated
/****************************************************************************/
/****************************************************************************/

struct InducedSpace_lw {
  // Array storing size of induced neighbor list, for each level of recursion.
  uintE* num_induced = nullptr;
  // Array storing induced neighbor list, for all levels of recursion. Levels go
  // from i=0 ... k-1. For each i, the induced neighbor list at level i is
  // induced + num_induced[0]*i (num_induced[0] is an upper-bound on the maximum
  // number of induced neighbors at any level).
  uintE* induced = nullptr;

  int* intersect = nullptr;
  InducedSpace_lw() {}

  void alloc(size_t max_deg, size_t k, size_t n) {
    if (!induced) induced = (uintE*)malloc(k * max_deg * sizeof(uintE));
    if (!num_induced) num_induced = (uintE*)malloc(k * sizeof(uintE));
    if (!intersect) {
      intersect = (int*)malloc(n * sizeof(int));
      for (size_t i = 0; i < n; i++) {
        intersect[i] = 0;
      }
    }
  }

  ~InducedSpace_lw() {
    if (induced) {
      free(induced);
      induced = nullptr;
    }
    if (num_induced) {
      free(num_induced);
      num_induced = nullptr;
    }
    if (intersect) {
      free(intersect);
      intersect = nullptr;
    }
  }
};

struct FullSpace_orig_lw {
  uintE* num_induced = nullptr;
  uintE* induced = nullptr;
  uintE* num_edges = 0;
  uintE* induced_edges = nullptr;
  uintE* induced_degs = nullptr;
  uintE* labels = nullptr;
  uintE* old_labels = nullptr;
  size_t nn = 0;

  FullSpace_orig_lw() {}

  void alloc(size_t max_induced, size_t k, size_t n) {
    induced = (uintE*)malloc(sizeof(uintE) * k * max_induced);
    induced_degs = (uintE*)malloc(sizeof(uintE) * k * max_induced);
    labels = (uintE*)malloc(sizeof(uintE) * max_induced);
    induced_edges = (uintE*)malloc(sizeof(uintE) * max_induced * max_induced);
    num_induced = (uintE*)malloc(sizeof(uintE) * k);
    num_edges = (uintE*)malloc(sizeof(uintE) * k);
    if (old_labels == nullptr) old_labels = (uintE*)calloc(n, sizeof(uintE));
  }

  template <class Graph>
  void setup(Graph& DG, size_t k, size_t i) {
    using W = typename Graph::weight_type;
    num_induced[0] = DG.get_vertex(i).out_degree();
    nn = num_induced[0];

    for (size_t j = 0; j < nn; j++) {
      induced[j] = j;
    }
    for (size_t j = 0; j < nn; j++) {
      induced_degs[j] = 0;
    }
    for (size_t j = 0; j < nn; j++) {
      labels[j] = 0;
    }

    size_t o = 0;
    auto map_label_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = o + 1;
      o++;
    };
    DG.get_vertex(i).out_neighbors().map(map_label_f, false);

    size_t j = 0;
    auto map_f = [&](const uintE& src, const uintE& v, const W& wgh) {

      auto map_nbhrs_f = [&](const uintE& src_v, const uintE& v_nbhr,
                             const W& wgh_v) {
        if (old_labels[v_nbhr] > 0) {
          induced_edges[j * nn + induced_degs[j]] = old_labels[v_nbhr] - 1;
          induced_degs[j]++;
        }
      };
      DG.get_vertex(v).out_neighbors().map(map_nbhrs_f, false);

      j++;
    };
    DG.get_vertex(i).out_neighbors().map(map_f, false);

    auto map_relabel_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = 0;
    };
    DG.get_vertex(i).out_neighbors().map(map_relabel_f, false);

    auto deg_seq = gbbs::make_slice(induced_degs, nn);
    num_edges[0] = parlay::reduce(deg_seq);
  }

  static void init() {}
  static void finish() {}

  ~FullSpace_orig_lw() {
    if (labels) {
      free(labels);
      labels = nullptr;
    }
    if (induced) {
      free(induced);
      induced = nullptr;
    }
    if (induced_edges) {
      free(induced_edges);
      induced_edges = nullptr;
    }
    if (induced_degs) {
      free(induced_degs);
      induced_degs = nullptr;
    }
    if (num_edges) {
      free(num_edges);
      num_edges = nullptr;
    }
    if (num_induced) {
      free(num_induced);
      num_induced = nullptr;
    }
    if (old_labels) {
      free(old_labels);
      old_labels = nullptr;
    }
  }
};

}  // namespace gbbs
