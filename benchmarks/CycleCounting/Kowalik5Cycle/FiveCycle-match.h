#include <algorithm>
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "gbbs/gbbs.h"
#include "gbbs/graph.h"
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include <unordered_map>


struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return pbbs::hash_combine(h1, h2);  
    }
};
using Unordered_map = std::unordered_map<std::pair<uintE, uintE>, int, pair_hash>;

template <class Graph>
inline ulong is_edge(uintE u, uintE v, Graph& G){
  auto deg = G.get_vertex(u).getOutDegree();
  if (deg == 0) return 0;
  uintE* nghs = (uintE*) G.get_vertex(u).getOutNeighbors();
  auto custom_less_v = [&] (uintE arg) { return v < arg;};
  auto index_leq_v = _binary_search(nghs, custom_less_v, deg);
  if (nghs[index_leq_v] == v && index_leq_v < deg) return 1;
  return 0;
}

template <class Graph>
inline ulong Count5CycleVertex_match(Graph& G, size_t i,
                                     Unordered_map& wedges){
    ulong count = 0;
    long s[5]= {i, -1, -1, -1, -1};
    auto v0 = G.get_vertex(i);
    uintE* nghs = (uintE*) v0.getOutNeighbors();
    for (uintE ei = 0; ei < v0.getOutDegree(); ei++) {
      s[1] = nghs[ei];
      auto v1 = G.get_vertex(s[1]);
      uintE* nghs_v1 = (uintE*) v1.getOutNeighbors();
      for (uintE ej = 0; ej < v0.getOutDegree(); ej++){
        if (ej == ei) continue;
        s[2] = nghs[ej];
        auto v2 = G.get_vertex(s[2]);
        uintE* nghs_v2 = (uintE*) v2.getOutNeighbors();
        for (uintE e1i = 0; e1i < v1.getOutDegree(); e1i++){
          if (nghs_v1[e1i] == s[2]) continue; 
          s[3] = nghs_v1[e1i];
          auto v3 = G.get_vertex(s[3]);
          uintE* nghs_v3 = (uintE*) v3.getOutNeighbors();
          count += wedges[std::make_pair(s[2] < s[3]? s[2] : s[3], s[2] < s[3] ? s[3] : s[2] )];
          count -= is_edge(i, s[2], G) & is_edge(i, s[3], G);
          count -= is_edge(s[1], s[2], G); 
          for (uintE e2i = 0; e2i < v2.getOutDegree(); e2i++){
            if (nghs_v2[e2i] == s[1] || nghs_v2[e2i] == s[3]) continue;
            count += is_edge(s[3], nghs_v2[e2i], G);
          }
          for (uintE e3i =0; e3i < v3.getOutDegree(); e3i++){
            if (nghs_v3[e3i] == s[2]) continue;
            count += is_edge(nghs_v3[e3i], s[2], G);
          }
        }
      }
    }
    return count;
}

template <class Graph>
inline ulong Count5Cycle_match(Graph& GA, long order_type = 0, size_t block_size = 500000, double epsilon = 0.1) {
  timer tp; tp.start();
  using W = typename Graph::weight_type;
  sequence<uintE> rank;
  sequence<uintT> order_to_vertex;

  if (order_type == 0) {
    rank = goodrichpszona_degen::DegeneracyOrder_intsort(GA, epsilon);
    order_to_vertex = sequence<uintT>(GA.n, [&](size_t i){return 0;} );
    par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t v)
                { order_to_vertex[rank[v]] = v; });
  } else if (order_type == 1) {
    rank = barenboimelkin_degen::DegeneracyOrder(GA, epsilon);
    order_to_vertex = sequence<uintT>(GA.n, [&](size_t i){return 0;} );
    par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t v)
               { order_to_vertex[rank[v]] = v; });
  } else if (order_type == 2) {
    order_to_vertex = orderNodesByDegree(GA, GA.n);
    rank = sequence<uintE>(GA.n, [&](size_t i){return 0;} );
    par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                { rank[order_to_vertex[i]] = i; });
  }
  debug(cout << "Rank abd Order done\n"; fflush(stdout););
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  debug(cout << "Relabel done\n"; fflush(stdout););

  auto out_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u < v;
  };
  auto OUTG = GDO.filterGraph(GDO, out_direction); // only keeps track of out edges
  debug(cout << "Filter done\n"; fflush(stdout););
  double ttp = tp.stop();
  debug(std::cout << "##### Preprocessing: " << ttp << std::endl;);

  timer ts; ts.start();
  auto parallel_work = sequence<size_t>(GA.n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return GDO.get_vertex(v).getOutDegree();
    };
    par_for(0, GA.n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = GDO.get_vertex(i).template reduceOutNgh<size_t>(i, map_f, monoid); // summing the degrees of the neighbors for each vertex?
    });
  }

  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());
  size_t n_blocks = total_work/block_size + 1;
  debug(std::cout << "##### number of blocks: " << n_blocks << std::endl;);
  size_t work_per_block = total_work / n_blocks;
  n_blocks = (total_work/work_per_block) + 1;
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  double tts = ts.stop();
  debug(std::cout << "##### Scheduling: " << tts << std::endl;);

  timer tc; tc.start();
  sequence<ulong> cycleCounts = sequence<ulong>(n_blocks * eltsPerCacheLine, [&](size_t s) { return 0; });
  Unordered_map inout_wedges;
  std::pair<uintE, uintE> wedge;
  for (uintT i = 0; i < OUTG.n; i++){
    auto vi = OUTG.get_vertex(i);
    uintE degree = vi.getOutDegree();
    uintE* nghs = (uintE*) vi.getOutNeighbors();
    for (uintE ei=0; ei < degree; ei++) {
      for (uintE ej=ei+1; ej < degree; ej++){
        wedge = std::make_pair(nghs[ej],nghs[ei]); // nghs[ej] < nghs[ei]
        inout_wedges[wedge]++;
      }
    }
  }

  auto run_block = [&](size_t start_ind, size_t end_ind, size_t block_index){
    for (size_t i = start_ind; i < end_ind; i++) { 
      ulong temp = Count5CycleVertex_match(OUTG, i, inout_wedges); 
      cycleCounts[block_index * eltsPerCacheLine] += temp;  
    }
  };

  par_for(0, n_blocks, [&] (size_t i) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = pbbslib::binary_search(parallel_work, start, less_fn);
    size_t end_ind = pbbslib::binary_search(parallel_work, end, less_fn);
    run_block(start_ind, end_ind, i);
  });

  ulong total = pbbslib::reduce_add(cycleCounts);
  double ttc = tc.stop();
  debug(std::cout << "##### Actual counting: " << ttc << std::endl;);
  std::cout << "##### Result: " << total << "," << ttp << "," << tts << "," << ttc << std::endl;

  GDO.del();
  OUTG.del();
  return total;
}