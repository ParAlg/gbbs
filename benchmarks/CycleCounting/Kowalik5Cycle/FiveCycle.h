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
#include "ligra/ligra.h"
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"


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
  // also keep INNEIGHBORS.

  // TODO: turn this back to parallel_for
  parallel_for(0, n, [&] (size_t i) {
  //for (size_t i  =0; i < n ; i++){
    // for vertex i in new graph, fill in its neighbor array.
    // then sort the neighbor array.
    w_vertex u = G.get_vertex(order_to_vertex[i]);
    size_t out_offset = outOffsets[i];
    uintE d = u.getOutDegree();
    if (d > 0) {
      edge* nghs = u.getOutNeighbors(); 
      edge* new_nghs = out_edges.begin() + out_offset;
      // parallelize???
      for (uintE j= 0; j < d; j++) {
        new_nghs[j] = std::make_tuple(vertex_to_order[std::get<0>(nghs[j])], std::get<1>(nghs[j]));
      }

      // neighbor with largest index first.
      pbbslib::sample_sort_inplace(pbbslib::make_sequence(new_nghs, d), cmp_by_dest_order);

      // auto pred_c = [&](const edge& e) {
      //   return pred(i, std::get<0>(e), std::get<1>(e));
      // };
      // auto n_im_f = [&](size_t i) { return nghs[i]; };
      // auto n_im = pbbslib::make_sequence<edge>(d, n_im_f);
      // pbbslib::filter_out(n_im, pbbslib::make_sequence(dir_nghs, d), pred_c, pbbslib::no_flag);
    }
  //}
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
      get_deletion_fn(out_vdata, out_edge_arr),
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
  // uintE* r = pbbslib::new_array_no_init<uintE>(n);
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

// constexpr const size_t _binary_search_base = 16;
// template <typename T>
// inline size_t reverse_binary_search(T* I, T v, size_t n){
//   //cout << "binary searching on: ";
//   // for (int j = 0; j < n ; j++) {
//   //   cout << I[j] << ", ";
//   // }
//   // cout<< endl;
//   size_t start = 0;
//   size_t end = n;
//   while (end-start > _binary_search_base) {
//     size_t mid = start + (end-start)/2;
//     if (!(I[mid]>v)) end = mid;
//     else start = mid + 1;
//   }
//   size_t i;
//   for (i = start; i < end; i++) {
//     if (!(I[i]>v)) {
//       //cout << v << " is found at index " << i << endl; 
//       return i;
//     }
//   }
//   //cout << v << " is found at index " << i << endl; 
//   if (i != n){
//     return i;
//   } else {
//     return n-1;
//   }
// }


// This file implements the 5-cycle counting algorithm from
//
// 
template <class Graph>
inline ulong Count5Cycle(Graph& GA, long order_type = 0, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;

  sequence<uintE> rank;
  //auto rank = sequence<uintE>(GA.n, [&](size_t s) { return s; });
  uintE n = GA.n; 
  // relabel the graph first. then do degeneracyorder

  auto order_to_vertex = orderNodesByDegree(GA, GA.n);
  cout << "Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering
  //auto GDO = GA;
  cout << "Relabel done\n"; fflush(stdout);

  if (order_type == 0) rank = goodrichpszona_degen::DegeneracyOrder_intsort(GDO, epsilon);
  else if (order_type == 1) rank = barenboimelkin_degen::DegeneracyOrder(GDO, epsilon);
  cout << "Rank done\n"; fflush(stdout);

  auto direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DGDO = filter_graph(GDO, direction); // only keeps track of out  uintEs
  cout << "Filter done\n"; fflush(stdout);

  timer t; t.start();
  auto parallel_work = sequence<size_t>(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return GDO.get_vertex(v).getOutDegree();
    };
    par_for(0, n, [&] (size_t i) {
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

  // sequence<timer> part_one_timers = sequence<timer>(72, [&](size_t s){return timer();});
  // sequence<timer> part_two_timers = sequence<timer>(72, [&](size_t s){return timer();});
  // sequence<timer> part_three_timers = sequence<timer>(72, [&](size_t s){return timer();});

  timer t2; t2.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  //sequence<ulong> cycleCounts = sequence<ulong>(GA.n * eltsPerCacheLine, [&](size_t s) { return 0; });
  sequence<ulong> cycleCounts = sequence<ulong>(n_blocks * eltsPerCacheLine, [&](size_t s) { return 0; });

  timer tp1, tp2, tp3, tp4; 

  auto run_intersection = [&](size_t start_ind, size_t end_ind, size_t block_index) {

    // timer tp1 = part_one_timers[worker_id()];
    // timer tp2 = part_two_timers[worker_id()];
    // timer tp3 = part_three_timers[worker_id()];
    ulong temp;
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
      temp = 0;
      auto vi = GDO.get_vertex(i);
      uintE degree = vi.getOutDegree();
      uintE* nghs = (uintE*) vi.getOutNeighbors(); 
      
      if (degree == 0) continue; 

      uintE viOutDegree = DGDO.get_vertex(i).getOutDegree();
      uintE* outnghs_vi = (uintE*) DGDO.get_vertex(i).getOutNeighbors();
      //auto outnghs_vi_seq =  sequence<uintE>(outnghs_vi, viOutDegree);
      //auto custom_less_i = [&](uintE arg) { return i < arg; }; // for binary searching.

      uintE u, uDegree, w, wOutDegree, x;
      
      if (i < 10) tp1.start(); // DEBUG

      for (uintE j = 0; (j < degree) && ((u = nghs[j]) > i); j++) {
        auto vu = GDO.get_vertex(u);
        uintE* nghs_u = (uintE*) vu.getOutNeighbors();
        uDegree = vu.getOutDegree();
        for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
          U[w] += 1;
        }
      } // end of line 7. 
      
      if (i < 10) tp1.stop(); // DEBUG

      //timer tp2 = part_two_timers[worker_id()];
      if (i < 10) tp2.start(); // DEBUG
     
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
          // auto custom_less_w = [&](uintE arg) { return w < arg; };
          // uintE index_leq_i, index_leq_w; 
          // // auto less_fn = std::less<size_t>();
          // if (wOutDegree > 0
          //       && (index_leq_i = pbbslib::binary_search(outnghs_w_seq, custom_less_i)) < wOutDegree 
          //        && outnghs_w[index_leq_i] == i) {
          //   w_vi_neighbors = 1;
          // } else if ((index_leq_w = pbbslib::binary_search(outnghs_vi_seq, custom_less_w)) < wOutDegree   
          //        && outnghs_vi[index_leq_w] == w) {
          //   w_vi_neighbors = 1;
          // }
          
          if (i < 10) tp4.start();

          for (uintE l = 0; (l < wOutDegree) && ((x = outnghs_w[l]) >= i); l++)  {
            if (x == i) {w_vi_neighbors = 1; break;}
          }
          for (uintE l = 0; (l < viOutDegree) && ((x = outnghs_vi[l]) >= w); l++)  {
            if (x == w) w_vi_neighbors = 1; 
          }
          if (i < 10) tp4.stop();

          if (i < 10) tp3.start(); // DEBUG
          for (uintE l = 0; (l < wOutDegree) && ((x = outnghs_w[l]) > i); l++) {
            if (x != u) {
               temp += U[x] - w_vi_neighbors;
            }
          }
          
          if (i < 10) tp3.stop(); // DEBUG
          
          //outnghs_w_seq.to_array();
        }

        for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
          U[w] += 1;
        }

      }
      cycleCounts[block_index * eltsPerCacheLine] += temp;  
      if (i < 10){
        tp2.stop();
        std::cout << "===== vertex " << i << " =====" << std::endl; 
        tp1.reportTotal("Part 1: ");
        tp2.reportTotal("Part 2: ");
        tp3.reportTotal("Part 3: ");
        tp4.reportTotal("Part 4: ");
        tp1.reset(); 
        tp2.reset();
        tp3.reset();
        tp4.reset();
      }
      //outnghs_vi_seq.to_array();

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


  // timer stuff....
  // double part_one_time = 0;
  // double part_two_time = 0;
  // double part_three_time = 0;
  // std::cout << "#### each timer for part 1: ";
  // for (int i = 0; i < 72; i++) {
  //   std::cout << part_one_timers[i].get_total() << ", ";
  //   part_one_time += part_one_timers[i].get_total();
  //   part_two_time += part_two_timers[i].get_total();
  //   part_three_time += part_three_timers[i].get_total();
  // }
  // std::cout << std::endl; 
  // std::cout << "#### Time spent in part 1: " << part_one_time << std::endl;
  // std::cout << "#### Time spent in part 2: " << part_two_time << std::endl;
  // std::cout << "#### Time spent in part 3: " << part_three_time << std::endl;


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
