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

      // // DEBUG
      //   cout << "vertex i has " << d << " neighbors, before reordering: ";
      //   for (uintE k = 0; k < d; k++){
      //     cout << std::get<0>(new_nghs[k]) << " ";
      //   }
      //   cout << endl;
      // // DEBUG

      // neighbor with largest index first.
      pbbslib::sample_sort_inplace(pbbslib::make_sequence(new_nghs, d), cmp_by_dest_order);


        // // DEBUG
        // cout << "vertex i has " << d << " neighbors, after reordering: ";
        // for (uintE k = 0; k < d; k++){
        //   cout << std::get<0>(new_nghs[k]) << " ";
        // }
        // cout << endl;
        // // DEBUG

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


// This file implements the 5-cycle counting algorithm from
//
// 
template <class Graph>
inline ulong Count5Cycle(Graph& GA, long order_type = 0, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  // using edge = typename Graph::edge_type; //std::tuple<uintE, W>;

  sequence<uintE> rank;
  uintE n = GA.n; 
  // relabel the graph first. then do degeneracyorder

  auto order_to_vertex = orderNodesByDegree(GA, GA.n);
  cout << "Order done\n"; fflush(stdout);
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  // for (size_t i = 0; i < 5; i++) {
  //   cout << "node " << i << " has " << GDO.get_vertex(i).getOutDegree() << " neighbors: ";
  //   for (size_t j = 0; j < GDO.get_vertex(i).getOutDegree(); j++) {
  //     cout << *((uintE*) GDO.get_vertex(i).getOutNeighbors() + j) << " ";
  //   }
  //   cout << endl;
  // }

  cout << "Relabel done\n"; fflush(stdout);

  if (order_type == 0) rank = goodrichpszona_degen::DegeneracyOrder_intsort(GDO, epsilon);
  else if (order_type == 1) rank = barenboimelkin_degen::DegeneracyOrder(GDO, epsilon);

  cout << "Rank done\n"; fflush(stdout);
  //double tt_rank = t_rank.stop();

  // rank[v] = degeneracy order of vertex v
  auto direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DGDO = filter_graph(GDO, direction); // only keeps track of out  uintEs
  cout << "Filter done\n"; fflush(stdout);

  // for (size_t i = 0; i < 5; i++) {
  //   cout << "filtered node " << i << " has " << DGDO.get_vertex(i).getOutDegree() << " neighbors: ";
  //   for (size_t j = 0; j < DGDO.get_vertex(i).getOutDegree(); j++) {
  //     cout << *((uintE*) DGDO.get_vertex(i).getOutNeighbors() + j) << " ";
  //   }
  //   cout << endl;
  // }


  // degree_order[i] = vertex label that has the i-th highest degree
  // sequence<uintE> degree_order = orderNodesByDegree(GA, GA.n);

  // degree_rank[i] = rank of vertex i by degree (high to low)
  // uintE* degree_rank = pbbslib::new_array_no_init<uintE>(n);
  // par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
  //                 { degree_rank[degree_order[i]] = i; });

  // // specify the point at which to 
  // using  uintE = std::tuple<uintE, W>;
  // sequence<uintE> valid_neighbor_start_index(n);
  // par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
  //   auto u = GDO.get_vertex(i);
  //   uintE outDegree_i = u.getOutDegree();
  //    uintE* nghs = u.getOutNeighbors(); // QUESTION: is this of type  uintE* or uintE????
  //    uintE* end = nghs + outDegree_i;
  //   uintE valid_index = 0;
  //   for ( uintE* x = nghs; valid_index < outDegree_i; x++, valid_index++){
  //     if (i < std::get<0>(x))
  //       break;
  //   }
  //   valid_neighbor_start_index[i] = valid_index;
  // });
  // => ALWAYS ITERATE FROM THE BACK AND BRERAK WHEN THE NEIGHBOR IS < i.

  sequence<ulong> cycleCounts = sequence<ulong>(GA.n, [&](size_t s) { return 0; });

  par_for(0, GA.n, pbbslib::kSequentialForThreshold,
   [&] (size_t i) { 

      auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });

      auto vi = GDO.get_vertex(i);
      uintE degree = vi.getOutDegree();
      uintE* nghs = (uintE*) vi.getOutNeighbors(); // QUESTION: is this of type  uintE* or uintE????

      uintE viOutDegree = DGDO.get_vertex(i).getOutDegree();
      uintE* outnghs_vi = (uintE*) DGDO.get_vertex(i).getOutNeighbors();

      uintE u, uDegree, w, wOutDegree, x;

      for (uintE j = 0; (j < degree) && ((u = nghs[j]) > i); j++) {
        auto vu = GDO.get_vertex(u);
        uintE* nghs_u = (uintE*) vu.getOutNeighbors();
        uDegree = vu.getOutDegree();
        for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
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

          // check if w->vi or vi->w is in DGDO
          for (uintE l = 0; (l < wOutDegree) && ((x = outnghs_w[l]) >= i); l++)  {
            if (x == i) w_vi_neighbors = 1; 
          }
          for (uintE l = 0; (l < viOutDegree) && ((x = outnghs_vi[l]) >= w); l++)  {
            if (x == w) w_vi_neighbors = 1; 
          }

          for (uintE l = 0; (l < wOutDegree) && ((x = outnghs_w[l]) > i); l++) {
            if (x != u) {
              cycleCounts[i] += U[x] - w_vi_neighbors;
            }
          }
        }

        for (uintE k = 0; (k < uDegree) && ((w = nghs_u[k]) > i); k++) {
          U[w] += 1;
        }

      }
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
  return total;
}
