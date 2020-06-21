// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <algorithm>
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "gbbs/gbbs.h"
#include "gbbs/graph.h"
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"


namespace gbbs {
struct U_FastReset {
  uintE* U = nullptr;
  uintE* distinct = nullptr;
  ulong num_distinct = 0; 

  void alloc(size_t n) {
    if (!U) U = (uintE*) calloc(n, sizeof(uintE)); //malloc(n * sizeof(uintE));

    if (!distinct) distinct = (uintE*) calloc(n, sizeof(uintE));
    // auto V = sequence<uintE>(n, [&](size_t s) { return 0; });
    // U = (uintE*) &V;
    // auto distinct_s = sequence<uintE>(n, [&](size_t s) { return 0; });
    // distinct = (uintE*) &distinct_s;
  }

  void reset(){
    for (ulong i = 0; i < num_distinct; i++){
      U[distinct[i]] = 0;
      distinct[i] = 0;
    }
    num_distinct = 0;
  }

  void del() {
    if (U) { free(U); U = nullptr; }
    if (distinct) { free(distinct); distinct = nullptr; }
  }

  ~U_FastReset() { del(); }
};

constexpr const size_t binary_search_base = 16;

template <typename T, typename F>
inline size_t _linear_search(T* I, const F& less, size_t n) {
  for (size_t i = 0; i < n; i++)
    if (!less(I[i])) return i;
  return n;
}

// return index to first key where less is false
template <typename T, typename F>
inline size_t _binary_search(T* I, const F& less, size_t n) {
  size_t start = 0;
  size_t end = n; 
  while (end - start > binary_search_base) {
    size_t mid = start + (end - start) / 2;
    if (!less(I[mid]))
      end = mid;
    else
      start = mid + 1;
  }
  return start + _linear_search(I + start, less, end - start);
}


// arg order_to_vertex[i] = vertex in original graph G that will be vertex i in the new graph.
// new_order[i] = the newe label for vertex i
/* Filters a symmetric graph, G, with a predicate function pred.  Note
 * that the predicate does not have to be symmetric, i.e. f(u,v) is
 * not necesssarily equal to f(v,u), but we only represent the out-edges of this
 * (possibly) directed graph. For convenience in cases where the graph needed is
 * symmetric, we coerce this to a symmetric_graph. */
template <template <class W> class vertex, class W,
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
inline symmetric_graph<symmetric_vertex, W> relabel_graph(symmetric_graph<vertex, W>& G, sequence<uintT>& order_to_vertex) {
  using w_vertex = vertex<W>;
  size_t n = G.n;

  // vertex_to_order[i] = the newe label (order) for vertex i
  sequence<uintE> vertex_to_order = sequence<uintE>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
    vertex_to_order[order_to_vertex[i]] = i;
  });

  auto outOffsets = sequence<uintT>(n + 1);
  outOffsets[n] = 0;
  for (size_t i = 0; i < n; i++) { // for vertex of rank/order i.
    outOffsets[i] = G.get_vertex(order_to_vertex[i]).getOutDegree();
  }
  //
  uintE outEdgeCount = pbbslib::scan_add_inplace(outOffsets);

  using edge = std::tuple<uintE, W>;
  auto cmp_by_dest_order = [](const edge& e1, const edge& e2) {
      return std::get<0>(e1) > std::get<0>(e2);
  };

  auto out_edges = sequence<edge>(outEdgeCount);
  parallel_for(0, n, [&] (size_t i) {
    w_vertex u = G.get_vertex(order_to_vertex[i]);
    size_t out_offset = outOffsets[i];
    uintE d = u.getOutDegree();
    if (d > 0) {
      edge* nghs = u.getOutNeighbors();
      edge* new_nghs = out_edges.begin() + out_offset;
      for (uintE j= 0; j < d; j++) {
        new_nghs[j] = std::make_tuple(vertex_to_order[std::get<0>(nghs[j])], std::get<1>(nghs[j]));
      }

      // neighbor with largest index first.
      pbbslib::sample_sort_inplace(pbbslib::make_sequence(new_nghs, d), cmp_by_dest_order);
    }
  }, 1);

  auto out_vdata = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    out_vdata[i].offset = outOffsets[i];
    out_vdata[i].degree = outOffsets[i+1]-outOffsets[i];
  });
  outOffsets.clear();

  auto out_edge_arr = out_edges.to_array();
  return symmetric_graph<symmetric_vertex, W>(
      out_vdata, G.n, outEdgeCount,
      [=] { pbbslib::free_arrays(out_vdata, out_edge_arr); },
      out_edge_arr);
}

// TODO: add relabel_graph with different signatures
template <template <class W> class vertex, class W,
          typename std::enable_if<
              std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value,
              int>::type = 0>
inline auto relabel_graph(symmetric_graph<vertex, W>& G, sequence<uintT>& order_to_vertex) -> decltype(G) {
  std::cout << "Relabel graph not implemented for byte representation" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W,
    typename std::enable_if<std::is_same<vertex<W>, asymmetric_vertex<W>>::value,
                            int>::type = 0>
inline auto relabel_graph(asymmetric_graph<vertex, W>& G, sequence<uintT>& order_to_vertex) -> decltype(G) {
  std::cout << "Relabel graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline auto relabel_graph(asymmetric_graph<vertex, W>& G, sequence<uintT>& order_to_vertex) -> decltype(G) {
  std::cout << "Relabel graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

// highest degree go first
template <class Graph>
inline sequence<uintT> orderNodesByDegree(Graph& G, size_t n) {
  sequence<uintT> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });

  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).getOutDegree() > G.get_vertex(v).getOutDegree();
  });
  //par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
  //                { r[o[i]] = i; });
  t.stop();
  debug(t.reportTotal("Rank time"););
  return o;
}


// ****************************************************************************************
// Relabels graph by degree and return it
// ****************************************************************************************
template <class Graph>
inline auto Preprocess(Graph& GA, sequence<uintE>& rank, long order_type = 0, double epsilon = 0.1) {
  //using W = typename Graph::weight_type;
  // sequence<uintE> rank;
  // relabel the graph first. then do degeneracyorder
  auto order_to_vertex = orderNodesByDegree(GA, GA.n);
  cout << "Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering
  cout << "Relabel done\n"; fflush(stdout);

  if (order_type == 0) rank = goodrichpszona_degen::DegeneracyOrder_intsort(GDO, epsilon);
  else if (order_type == 1) rank = barenboimelkin_degen::DegeneracyOrder(GDO, epsilon);
  else if (order_type == 2) rank = sequence<uintE>(GA.n, [&](size_t s) { return s; });
  cout << "Rank done\n"; fflush(stdout);
  return GDO;
}

// ****************************************************************************************
// Count the number of 5-cycles for a vertex in teh Kowalik algorithm
// ****************************************************************************************
template <class Graph, class Graph2>
inline ulong Count5CycleVertex(Graph& GDO, Graph2& DGDO, U_FastReset* V, size_t i){
  ulong temp = 0;
  uintE* U = V->U;
  auto vi = GDO.get_vertex(i);
  uintE degree = vi.getOutDegree();
  uintE* nghs = (uintE*) vi.getOutNeighbors(); 
  
  if (degree == 0) return 0; 

  // uintE viOutDegree = DGDO.get_vertex(i).getOutDegree();
  // uintE* outnghs_vi = (uintE*) DGDO.get_vertex(i).getOutNeighbors();
  // auto outnghs_vi_seq =  sequence<uintE>(outnghs_vi, viOutDegree);
  // auto custom_less_i = [&](uintE arg) { return i < arg; }; // for binary searching.

  uintE u, uDegree, w, wOutDegree, x;


  for (uintE j = 0; (j < degree) && ((u = nghs[j]) > i); j++) {
    auto vu = GDO.get_vertex(u);
    uintE* nghs_u = (uintE*) vu.getOutNeighbors();
    uDegree = vu.getOutDegree();
    for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
      if (U[w] == 0) { V->distinct[V->num_distinct] = w; V->num_distinct++; }
      U[w] += 1;
    }
  } // end of line 7. 
  
  for (uintE j = 0; (j < degree) && ((u = nghs[j]) > i); j++) {
    auto vu = GDO.get_vertex(u);
    uintE* nghs_u = (uintE*) vu.getOutNeighbors();
    uDegree = vu.getOutDegree();
    for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
      U[w] -= 1;
    }

    for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
      auto vw = DGDO.get_vertex(w);
      uintE* outnghs_w = (uintE*) vw.getOutNeighbors();
      wOutDegree = vw.getOutDegree();
      int w_vi_neighbors = 0;

      // This part check if w->vi or vi->w is in DGDO, hence w_vi_neighbors
      // TODO: binary search!

      // stuff for binary search
      // auto outnghs_w_seq =  sequence<uintE>(outnghs_w, wOutDegree);
      auto custom_less_w = [&](uintE arg) { return w < arg; };
      //uintE index_leq_i, index_leq_w; 
      // // auto less_fn = std::less<size_t>();
      // if (wOutDegree > 0
      //       && (index_leq_i = _binary_search(outnghs_w, custom_less_i, wOutDegree)) < wOutDegree 
      //        && outnghs_w[index_leq_i] == i) {
      //   w_vi_neighbors = 1;
      // } else if ((index_leq_w = _binary_search(outnghs_vi, custom_less_w, viOutDegree)) < viOutDegree   
      //        && outnghs_vi[index_leq_w] == w) {
      //   w_vi_neighbors = 1;
      // }
      uintE index_leq_w;
      if ((index_leq_w = _binary_search(nghs, custom_less_w, degree)) < degree 
             && nghs[index_leq_w] == w) {
        w_vi_neighbors = 1;
      }



      // for (uintE l = 0; (l < wOutDegree) && ((x = outnghs_w[l]) >= i); l++)  {
      //   if (x == i) {w_vi_neighbors = 1; break;}
      // }
      // for (uintE l = 0; (l < viOutDegree) && ((x = outnghs_vi[l]) >= w); l++)  {
      //   if (x == w) w_vi_neighbors = 1; 
      // }
      
      for (uintE l = 0; (l < wOutDegree) && ((x = outnghs_w[l]) > i); l++) {
        if (x != u) {
            temp += U[x] - w_vi_neighbors;
        }
      }

      // auto wnghs = outnghs_w_seq.to_array();
    }

    for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
      U[w] += 1;
    }

  }

  // auto inghs = outnghs_vi_seq.to_array();

  V->reset();
  return temp;
}





// ****************************************************************************************
// This implements the 5-cycle counting algorithm from Kowalik
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle(Graph& GA, long order_type = 0, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;
  std::cout << "INFO: testing resetting U instead of reallocating (parallel with scheduling)." << std::endl;
  sequence<uintE> rank;
  auto GDO = Preprocess(GA, rank, order_type, epsilon);

  auto direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DGDO = GDO.filterGraph(GDO, direction); // only keeps track of out  uintEs
  cout << "Filter done\n"; fflush(stdout);


  // Start of work scheduling code.
  timer t; t.start();
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

  // auto parallel_work = sequence<size_t>(GA.n);
  // {
  //   auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
  //     return pre_parallel_work[v];
  //   };
  //   par_for(0, GA.n, [&] (size_t i) {
  //     auto monoid = pbbslib::addm<size_t>();
  //     parallel_work[i] = GDO.get_vertex(i).template reduceOutNgh<size_t>(i, map_f, monoid); // summing the degrees of the neighbors for each vertex?
  //   });
  // }

  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

  size_t block_size = 5000000;
  size_t n_blocks = total_work/block_size + 1;
  std::cout << "##### number of blocks: " << n_blocks << std::endl;
  size_t work_per_block = total_work / n_blocks;
  n_blocks = (total_work/work_per_block) + 1;
  double tt = t.stop();
  std::cout << "##### Scheduling: " << tt << std::endl;

  timer t2; t2.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  //sequence<ulong> cycleCounts = sequence<ulong>(GA.n * eltsPerCacheLine, [&](size_t s) { return 0; });
  sequence<ulong> cycleCounts = sequence<ulong>(n_blocks * eltsPerCacheLine, [&](size_t s) { return 0; });

  auto run_intersection = [&](size_t start_ind, size_t end_ind, size_t block_index, U_FastReset* V){ //sequence<uintE>* V) {

    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      ulong temp = Count5CycleVertex(GDO, DGDO, V, i); 
      cycleCounts[block_index * eltsPerCacheLine] += temp;  

      // for (uintE j = 0; (j < degree) && ((u = nghs[j]) > i); j++) {
      //   auto vu = GDO.get_vertex(u);
      //   uintE* nghs_u = (uintE*) vu.getOutNeighbors();
      //   uDegree = vu.getOutDegree();
      //   for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
      //     U[w] = 0; 
      //   }
      // } 
    }
  };

  // auto init_U = [&](sequence<uintE>* U){ *U = sequence<uintE>(GA.n, [&](size_t s) { return 0; }); };
  // auto finish_U = [&](sequence<uintE>* U){ return; };
  auto init_V = [&](U_FastReset* V){ V->alloc(GA.n); };
  auto finish_V = [&](U_FastReset* V){ return; };
  parallel_for_alloc<U_FastReset>(init_V, finish_V, 0, n_blocks,   // Testing space reuse
                                      [&] (size_t i, U_FastReset* V) {
  // par_for(0, n_blocks, 1, [&] (size_t i) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = pbbslib::binary_search(parallel_work, start, less_fn);
    size_t end_ind = pbbslib::binary_search(parallel_work, end, less_fn);
    run_intersection(start_ind, end_ind, i, V);
  });

  ulong total = pbbslib::reduce_add(cycleCounts);
  // ulong total = 0;
  // ulong temp;
  // std::cout << "### Number of neighbors of vertex, " << std::endl;

  // for (size_t i = 0; i < GA.n; i++) {
  //   std::cout << GDO.get_vertex(i).getOutDegree() << ", ";
  //   // temp = total + cycleCounts[i];
  //   // assert(ULONG_MAX - cycleCounts[i] >= total);
  //   // total = total + cycleCounts[i];
  // }
  // std::cout << std::endl;
  double tt2 = t2.stop();
  std::cout << "##### Actual counting: " << tt2 << std::endl;
  GDO.del();
  DGDO.del();
  return total;
}



// ****************************************************************************************
// This is a PURE SERIAL version of the five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_serial(Graph& GA, long order_type = 0, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;

  std::cout << "INFO: testing resetting U instead of reallocating." << std::endl;

  sequence<uintE> rank;
  //auto rank = sequence<uintE>(GA.n, [&](size_t s) { return s; });
  // relabel the graph first. then do degeneracyorder

  auto order_to_vertex = orderNodesByDegree(GA, GA.n);
  cout << "Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering
  //auto GDO = GA;
  cout << "Relabel done\n"; fflush(stdout);

  if (order_type == 0) rank = goodrichpszona_degen::DegeneracyOrder_intsort(GDO, epsilon);
  else if (order_type == 1) rank = barenboimelkin_degen::DegeneracyOrder(GDO, epsilon);
  else if (order_type == 2) rank = sequence<uintE>(GA.n, [&](size_t s) { return s; });
  cout << "Rank done\n"; fflush(stdout);

  auto direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DGDO = GDO.filterGraph(GDO, direction); // only keeps track of out  uintEs
  cout << "Filter done\n"; fflush(stdout);

  timer t; t.start();
  ulong cycleCount = 0;

  U_FastReset V1; 
  U_FastReset* V = &V1; 
  V->alloc(GA.n);
  for (size_t i = 0; i < GA.n; i++) { 
    auto temp = Count5CycleVertex(GDO, DGDO, V, i);
    cycleCount += temp;
  }

  double tt = t.stop();
  std::cout << "##### Actual counting: " << tt << std::endl;
  GDO.del();
  DGDO.del();
  return cycleCount;
}


// ****************************************************************************************
// This is a PARALLEL WITHOUT WORK SCHEDULING version of the five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_no_scheduling(Graph& GA, long order_type = 0, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;
  std::cout << "INFO: testing resetting U instead of reallocating." << std::endl;

  sequence<uintE> rank;
  auto GDO = Preprocess(GA, rank, order_type, epsilon);

  auto direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DGDO = GDO.filterGraph(GDO, direction); // only keeps track of out  uintEs
  cout << "Filter done\n"; fflush(stdout);

  timer t2; t2.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  sequence<ulong> cycleCounts = sequence<ulong>(GA.n * eltsPerCacheLine, [&](size_t s) { return 0; });


  auto init_V = [&](U_FastReset* V){ V->alloc(GA.n); };
  auto finish_V = [&](U_FastReset* V){  if (V != nullptr) { delete V; } };

  parallel_for_alloc<U_FastReset>(init_V, finish_V, 0, GA.n,   // Testing space reuse
                                      [&] (size_t i, U_FastReset* V) {
 
      auto temp = Count5CycleVertex(GDO, DGDO, V, i);
      cycleCounts[i * eltsPerCacheLine] += temp;

  });

  ulong total = pbbslib::reduce_add(cycleCounts);

  double tt2 = t2.stop();
  std::cout << "##### Actual counting: " << tt2 << std::endl;
  GDO.del();
  DGDO.del();
  return total;
}


// ****************************************************************************************
// This is a EXPERIMENTAL version of the five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_experiment(Graph& GA, long order_type = 0, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  sequence<uintE> rank ;
  // relabel the graph first. then do degeneracyorder
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
  cout << "Rank abd Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  //auto GDO = GA;
  cout << "Relabel done\n"; fflush(stdout);

  auto direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    //return rank[u] < rank[v];
    return u < v;
  };
  auto DGDO = GDO.filterGraph(GDO, direction); // only keeps track of out  uintEs

  cout << "Filter done\n"; fflush(stdout);

  timer t; t.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  sequence<ulong> cycleCounts = sequence<ulong>(72 * eltsPerCacheLine, [&](size_t s) { return 0; });


  par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
    auto vi = GDO.get_vertex(i);
    uintE degree = vi.getOutDegree();
    // uintE* nghs = (uintE*) vi.getOutNeighbors();

    if (degree == 0) return;

    uintE viOutDegree  = DGDO.get_vertex(i).getOutDegree();
    uintE* outnghs_vi = (uintE*) DGDO.get_vertex(i).getOutNeighbors();

    uintE u, uDegree, w, wOutDegree, x;

    //for (uintE j = 0; (j < degree) /*&& ((u = nghs[j]) > i)*/; j++) {
    for (uintE j = 0; j < viOutDegree; j++) {
      u = outnghs_vi[j]; // nghs[j];
      auto vu = GDO.get_vertex(u);
      uintE* nghs_u = (uintE*) vu.getOutNeighbors();
      uDegree = vu.getOutDegree();
      for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
        U[w] += 1;
      }
    }
    //for (uintE j = 0; (j < degree) /*&& ((u = nghs[j]) > i) */; j++) {
    for (uintE j = 0; j < viOutDegree; j++) {
      u = outnghs_vi[j];  //nghs[j];
      auto vu = GDO.get_vertex(u);
      uintE* nghs_u = (uintE*) vu.getOutNeighbors();
      uDegree = vu.getOutDegree();
      for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
        U[w] -= 1;
      }

      // TODO: maybe optimize this chunk further
      for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
        w = nghs_u[k];
        auto vw = DGDO.get_vertex(w);
        uintE* outnghs_w = (uintE*) vw.getOutNeighbors();
        wOutDegree = vw.getOutDegree();
        int w_vi_neighbors = 0;

        // for (uintE l = 0; (l < wOutDegree) /* && ((x = outnghs_w[l]) >= i) */ ; l++)  {
        //   if (x == i) {w_vi_neighbors = 1; break;}
        // }
        for (uintE l = 0; (l < viOutDegree) && ((x = outnghs_vi[l]) >= w); l++)  {
          if (x == w) w_vi_neighbors = 1;
        }

        for (uintE l = 0; (l < wOutDegree) /*&& ((x = outnghs_w[l]) > i)*/; l++) {
          // x has to be greater than i.
          if ((x = outnghs_w[l]) != u) {
              cycleCounts[worker_id()] += U[x] - w_vi_neighbors;
          }
        }
      }

      for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
        w = nghs_u[k];
        U[w] += 1;
      }
    }

  });

  ulong cycleCount = pbbslib::reduce_add(cycleCounts);

  double tt = t.stop();
  std::cout << "##### Actual counting: " << tt << std::endl;
  GDO.del();
  DGDO.del();
  return cycleCount;
}



// ****************************************************************************************
// This is the ESCAPE five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_ESCAPE(Graph& GA, long order_type = 0, double epsilon = 0.1) {

  //std::cout << "JUST TO MAKE SURE I AM REALLY USING ESCAPE" << std::endl;
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;

  //auto GA2 = relabel_graph(GA, sequence<uintT>(GA.n, [&](size_t i){return i;}));
  sequence<uintE> rank ;
  //auto rank = sequence<uintE>(GA.n, [&](size_t s) { return s; });
  // relabel the graph first. then do degeneracyorder
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
  cout << "Rank abd Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  //auto GDO = GA;
  cout << "Relabel done\n"; fflush(stdout);

  auto out_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    //return rank[u] < rank[v];
    return u < v;
  };
  auto OUTG = GDO.filterGraph(GDO, out_direction); // only keeps track of out edges

  auto in_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    //return rank[u] < rank[v];
    return u > v;
  };
  auto ING = GDO.filterGraph(GDO, in_direction); // only keeps track of in edges

  cout << "Filter done\n"; fflush(stdout);

  timer t; t.start();
  ulong cycleCount = 0;
  auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
  //par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
  for (uintE i = 0; i < GA.n; i++) {
    ulong tmp = 0;
    //auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
    auto vi = GDO.get_vertex(i);
    uintE degree = vi.getOutDegree();
    uintE* nghs = (uintE*) vi.getOutNeighbors();


    if (degree == 0) continue; //return;
    auto nghs_seq = sequence<uintE>(nghs, degree);


    uintE viOutDegree  = OUTG.get_vertex(i).getOutDegree();
    uintE* outnghs_vi = (uintE*) OUTG.get_vertex(i).getOutNeighbors();
    uintE viInDegree  = ING.get_vertex(i).getOutDegree();
    uintE* innghs_vi = (uintE*) ING.get_vertex(i).getOutNeighbors();

    uintE vj, vk, vl;
    // i -> j -> k
    for (uintE j = 0; j < viOutDegree; j++){
      vj = outnghs_vi[j];
      uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
      uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

      for (uintE k = 0; k < vjOutDegree;  k++) {
        vk = outnghs_vj[k];
        U[vk] += 1;
      }

    }

    // i -> j <- k or i -> j -> k
    for (uintE j = 0; j < viInDegree; j++){
      vj = innghs_vi[j];
      uintE* nghs_vj = (uintE*) GDO.get_vertex(vj).getOutNeighbors();
      uintE vj_degree = GDO.get_vertex(vj).getOutDegree();
      for (uintE k = 0; k < vj_degree; k++){
        U[nghs_vj[k]] += 1;
      }
    }

    U[i] = 0;


    for (uintE j = 0; j < viInDegree; j++) {
      vj = innghs_vi[j];
      uintE* innghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
      uintE vjInDegree = ING.get_vertex(vj).getOutDegree();
      uintE* nghs_vj = (uintE*) GDO.get_vertex(vj).getOutNeighbors();
      uintE vj_degree = GDO.get_vertex(vj).getOutDegree();

      if (vj_degree == 0) continue;
      auto nghs_vj_seq = sequence<uintE>(nghs_vj, vj_degree);


      for (uintE k = 0; k < vjInDegree; k++) {
        vk = innghs_vj[k];
        uintE* outnghs_vk = (uintE*) OUTG.get_vertex(vk).getOutNeighbors();
        uintE vkOutDegree = OUTG.get_vertex(vk).getOutDegree();
        auto custom_less_k = [&](uintE arg) { return vk < arg; };
        for (uintE ell = 0; ell < vkOutDegree; ell++) {
          vl = outnghs_vk[ell];
          if (vl == i || vl == vj) continue;
          tmp += U[vl];

          auto custom_less_l = [&](uintE arg) { return vl < arg; };

          uintE index_leq_vk, index_leq_vl;
          if ( ((index_leq_vk = pbbslib::binary_search(nghs_seq, custom_less_k)) < degree )
              && (nghs_seq[index_leq_vk] == vk))  tmp--;

          if (((index_leq_vl = pbbslib::binary_search(nghs_vj_seq, custom_less_l)) < vj_degree)
            && (nghs_vj_seq[index_leq_vl] == vl)) tmp--;
        }

      }

      nghs_vj = nghs_vj_seq.to_array();

    }

    // cycleCounts[worker_id()] += tmp;
    cycleCount += tmp;


    // U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
    for (uintE j = 0; j < viOutDegree; j++){
      vj = outnghs_vi[j];
      uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
      uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

      for (uintE k = 0; k < vjOutDegree;  k++) {
        vk = outnghs_vj[k];
        U[vk] = 0;
      }

    }

    // i -> j <- k or i -> j -> k
    for (uintE j = 0; j < viInDegree; j++){
      vj = innghs_vi[j];
      uintE* nghs_vj = (uintE*) GDO.get_vertex(vj).getOutNeighbors();
      uintE vj_degree = GDO.get_vertex(vj).getOutDegree();
      for (uintE k = 0; k < vj_degree; k++){
        U[nghs_vj[k]] = 0;
      }
    }

    nghs = nghs_seq.to_array();
  }

  //ulong cycleCount = pbbslib::reduce_add(cycleCounts);

  double tt = t.stop();
  std::cout << "##### Actual counting: " << tt << std::endl;
  GDO.del();
  ING.del();
  OUTG.del();
  return cycleCount;
}



// TODO: Make this parallel with work-scheduling
// ****************************************************************************************
// This is the ESCAPE five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_ESCAPE_par(Graph& GA, long order_type = 0, double epsilon = 0.1) {

  //std::cout << "JUST TO MAKE SURE I AM REALLY USING ESCAPE" << std::endl;
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;

  //auto GA2 = relabel_graph(GA, sequence<uintT>(GA.n, [&](size_t i){return i;}));
  sequence<uintE> rank ;
  //auto rank = sequence<uintE>(GA.n, [&](size_t s) { return s; });
  // relabel the graph first. then do degeneracyorder
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
  cout << "Rank abd Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  //auto GDO = GA;
  cout << "Relabel done\n"; fflush(stdout);

  auto out_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    //return rank[u] < rank[v];
    return u < v;
  };
  auto OUTG = GDO.filterGraph(GDO, out_direction); // only keeps track of out edges

  auto in_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    //return rank[u] < rank[v];
    return u > v;
  };
  auto ING = GDO.filterGraph(GDO, in_direction); // only keeps track of in edges

  cout << "Filter done\n"; fflush(stdout);


  timer t; t.start();
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

  size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
  size_t work_per_block = total_work / n_blocks;
  n_blocks = (total_work/work_per_block) + 1;
  double tt = t.stop();
  std::cout << "##### Scheduling: " << tt << std::endl;


  timer t2; t2.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  //sequence<ulong> cycleCounts = sequence<ulong>(72 * eltsPerCacheLine, [&](size_t s) { return 0; });
  sequence<ulong> cycleCounts = sequence<ulong>(n_blocks * eltsPerCacheLine, [&](size_t s) { return 0; });

  auto run_intersection = [&](size_t start_ind, size_t end_ind, size_t block_index) {
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
      ulong tmp = 0;
      //auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
      auto vi = GDO.get_vertex(i);
      uintE degree = vi.getOutDegree();
      uintE* nghs = (uintE*) vi.getOutNeighbors();


      if (degree == 0) continue; //return;
      auto nghs_seq = sequence<uintE>(nghs, degree);


      uintE viOutDegree  = OUTG.get_vertex(i).getOutDegree();
      uintE* outnghs_vi = (uintE*) OUTG.get_vertex(i).getOutNeighbors();
      uintE viInDegree  = ING.get_vertex(i).getOutDegree();
      uintE* innghs_vi = (uintE*) ING.get_vertex(i).getOutNeighbors();

      uintE vj, vk, vl;
      // i -> j -> k
      for (uintE j = 0; j < viOutDegree; j++){
        vj = outnghs_vi[j];
        uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

        for (uintE k = 0; k < vjOutDegree;  k++) {
          vk = outnghs_vj[k];
          U[vk] += 1;
        }
      }
      // i -> j <- k or i -> j -> k
      for (uintE j = 0; j < viInDegree; j++){
        vj = innghs_vi[j];
        uintE* nghs_vj = (uintE*) GDO.get_vertex(vj).getOutNeighbors();
        uintE vj_degree = GDO.get_vertex(vj).getOutDegree();
        for (uintE k = 0; k < vj_degree; k++){
          U[nghs_vj[k]] += 1;
        }
      }
      U[i] = 0;

      for (uintE j = 0; j < viInDegree; j++) {
        vj = innghs_vi[j];
        uintE* innghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
        uintE vjInDegree = ING.get_vertex(vj).getOutDegree();
        uintE* nghs_vj = (uintE*) GDO.get_vertex(vj).getOutNeighbors();
        uintE vj_degree = GDO.get_vertex(vj).getOutDegree();

        if (vj_degree == 0) continue;
        auto nghs_vj_seq = sequence<uintE>(nghs_vj, vj_degree);

        for (uintE k = 0; k < vjInDegree; k++) {
          vk = innghs_vj[k];
          uintE* outnghs_vk = (uintE*) OUTG.get_vertex(vk).getOutNeighbors();
          uintE vkOutDegree = OUTG.get_vertex(vk).getOutDegree();
          auto custom_less_k = [&](uintE arg) { return vk < arg; };
          for (uintE ell = 0; ell < vkOutDegree; ell++) {
            vl = outnghs_vk[ell];
            if (vl == i || vl == vj) continue;
            tmp += U[vl];

            auto custom_less_l = [&](uintE arg) { return vl < arg; };

            uintE index_leq_vk, index_leq_vl;
            if ( ((index_leq_vk = pbbslib::binary_search(nghs_seq, custom_less_k)) < degree )
                && (nghs_seq[index_leq_vk] == vk))  tmp--;

            if (((index_leq_vl = pbbslib::binary_search(nghs_vj_seq, custom_less_l)) < vj_degree)
              && (nghs_vj_seq[index_leq_vl] == vl)) tmp--;
          }
        }
        nghs_vj = nghs_vj_seq.to_array();
      }
      // cycleCounts[worker_id()] += tmp;
      cycleCounts[block_index * eltsPerCacheLine] += tmp;


      // U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
      for (uintE j = 0; j < viOutDegree; j++){
        vj = outnghs_vi[j];
        uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

        for (uintE k = 0; k < vjOutDegree;  k++) {
          vk = outnghs_vj[k];
          U[vk] = 0;
        }

      }

      // i -> j <- k or i -> j -> k
      for (uintE j = 0; j < viInDegree; j++){
        vj = innghs_vi[j];
        uintE* nghs_vj = (uintE*) GDO.get_vertex(vj).getOutNeighbors();
        uintE vj_degree = GDO.get_vertex(vj).getOutDegree();
        for (uintE k = 0; k < vj_degree; k++){
          U[nghs_vj[k]] = 0;
        }
      }

      nghs = nghs_seq.to_array();
    }
  };


  par_for(0, n_blocks, 1, [&] (size_t i) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = pbbslib::binary_search(parallel_work, start, less_fn);
    size_t end_ind = pbbslib::binary_search(parallel_work, end, less_fn);
    run_intersection(start_ind, end_ind, i);
  });

  ulong cycleCount = pbbslib::reduce_add(cycleCounts);

  double tt2 = t2.stop();
  std::cout << "##### Actual counting: " << tt2 << std::endl;
  GDO.del();
  ING.del();
  OUTG.del();
  return cycleCount;
}
}  // namespace gbbs