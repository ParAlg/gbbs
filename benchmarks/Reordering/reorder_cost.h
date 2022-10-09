#include "gbbs/gbbs.h"

namespace gbbs {

template <class Graph>
parlay::sequence<uintE> DegreeOrdering(const Graph& G) {
  auto by_degree = parlay::sequence<std::pair<uintE, uintE>>::from_function(
      G.num_vertices(), [&](size_t i) {
        return std::make_pair(G.get_vertex(i).out_degree(), (uintE)i);
      });
  parlay::sort_inplace(parlay::make_slice(by_degree));

  auto scattered = parlay::sequence<uintE>(G.num_vertices());
  parallel_for(0, G.num_vertices(), [&](size_t i) {
    auto vtx_id = by_degree[i].second;
    scattered[vtx_id] = i;
  });
  return scattered;
}

template <class Graph, class Perm>
double LogCost(const Graph& G, const Perm& p) {
  size_t n = G.num_vertices();
  auto logs = sequence<double>(n, 0.0);
  using W = typename Graph::weight_type;

  parallel_for(0, n, [&](size_t i) {
    auto v = G.get_vertex(i);
    double cost = 0.0;
    long id_u = p[i];
    auto map_f = [&](uintE u, uintE v, const W& wgh) {
      long id_v = p[v];
      cost += log((double)abs(id_v - id_u) + 1);
    };
    v.out_neighbors().map(map_f, false);
    logs[i] = cost;
  });
  double log_cost = parlay::reduce(logs) / (G.num_edges() * log(2.0));
  return log_cost;
}

template <class Graph, class Perm>
double LogGapCost(const Graph& G, const Perm& p) {
  size_t n = G.num_vertices();
  auto logs = sequence<double>(n, 0.0);
  using W = typename Graph::weight_type;

  parallel_for(0, n, [&](size_t i) {
    auto v = G.get_vertex(i);
    double cost = 0.0;
    size_t d = v.out_degree();
    auto nghs = parlay::sequence<uintE>(d);
    size_t k = 0;
    auto map_f = [&](uintE u, uintE v, const W& wgh) {
      long id_v = p[v];
      nghs[k] = id_v;
      ++k;
    };
    v.out_neighbors().map(map_f, false);
    parlay::sort_inplace(parlay::make_slice(nghs));

    uintE prev = 0;
    double our_id = p[i];
    for (size_t j = 0; j < nghs.size(); ++j) {
      double ngh_id = nghs[j];
      if (i == 0) {
        cost += log(abs(our_id - ngh_id) + 1);
        prev = ngh_id;
      } else {
        double gap = ngh_id - prev;
        cost += log(gap + 1);
        prev = ngh_id;
      }
    }
    logs[i] = cost;
  });
  double log_cost = parlay::reduce(logs) / (G.num_edges() * log(2.0));
  return log_cost;
}

template <class Graph, class Perm>
void PrintOrderingCosts(const Graph& G, const Perm& p) {
  auto random_order = parlay::random_permutation(G.num_vertices());
  auto orig_order =
      parlay::delayed_tabulate<>(G.num_vertices(), [&](size_t i) { return i; });
  auto degree_order = DegreeOrdering(G);
  std::cout << "LogCost(p) = " << LogCost(G, p) << std::endl;
  std::cout << "LogCost(orig_order) = " << LogCost(G, orig_order) << std::endl;
  std::cout << "LogCost(degree_order) = " << LogCost(G, degree_order)
            << std::endl;
  std::cout << "LogCost(random_perm) = " << LogCost(G, random_order)
            << std::endl;

  std::cout << "LogGapCost(p) = " << LogGapCost(G, p) << std::endl;
  std::cout << "LogGapCost(orig_order) = " << LogGapCost(G, orig_order)
            << std::endl;
  std::cout << "LogGapCost(degree_order) = " << LogGapCost(G, degree_order)
            << std::endl;
  std::cout << "LogGapCost(random_perm) = " << LogGapCost(G, random_order)
            << std::endl;
}

}  // namespace gbbs
