// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
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

#pragma once

#include <math.h>

#include "bucket.h"
#include "edge_map_reduce.h"
#include "ligra.h"
#include "pbbslib/dyn_arr.h"
#include "pbbslib/list_allocator.h"
#include "intersect.h"
#include "radix_wrapper.h"
#include "benchmarks/ApproximateDensestSubgraph/GreedyCharikar/DensestSubgraph.h"
#include "benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/DensestSubgraph.h"

template <class vertex>
inline uintE* rankNodes(vertex* V, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n);
  sequence<uintE> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return V[u].getOutDegree() < V[v].getOutDegree();
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  debug(t.reportTotal("Rank time"););
  return r;
}

// TODO densest subgraph paper to approximate alpha, 
template<class Graph>
inline sequence<uintE> DensestAppDegenOrder(Graph& GA, double epsilon=0.1, bool approx=false) {
  double alpha = approx ? CharikarAppxDensestSubgraph(GA) : WorkEfficientDensestSubgraph(GA, epsilon);
  const size_t n = GA.n;
  const size_t deg_cutoff = std::max((size_t) (ceil(alpha * epsilon)), (size_t) 1);
  auto sortD = sequence<uintE>(n, [&](size_t i) {
    return i;
  });
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).getOutDegree(); });
  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto get_deg =
      [&](uintE& p) -> uintE { return D[p] < deg_cutoff; };
  size_t start = 0;
  while (start < n) {
    // move all vert with deg < deg_cutoff in the front
    //integer_sort_inplace(sortD.slice(start, n), get_deg);
    radix::parallelIntegerSort(sortD.begin() + start, n - start, get_deg);
    auto BS = pbbs::delayed_seq<size_t>(n - start, [&] (size_t i) -> size_t {
      return D[sortD[i + start]] < deg_cutoff ? i + start : 0;});
    size_t end = pbbs::reduce(BS, pbbs::maxm<size_t>());
    if (end == start) end++; //TODO step?
    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
        D[v] -= edgesRemoved;
      return Maybe<std::tuple<uintE, uintE> >();
    };
    auto active =
        vertexSubset(n, end - start, sortD.begin() + start);
    auto moved = em.template edgeMapCount_sparse<uintE>(active, apply_f);

    moved.del();
    start = end;
  }
  auto ret = sequence<uintE>::no_init(n);
  parallel_for (0,n,[&] (size_t j) { ret[sortD[j]] = j; });
  return ret;
}

// Goodrich (2+epsilon) approx for degeneracy ordering where epsilon > 0
// Returns vertice sorted in degeneracy order
template<class Graph>
inline sequence<uintE> AppKCore(Graph& GA, double epsilon=0.001) {
  const size_t n = GA.n;
  const size_t ns = std::max((size_t) (ceil((n*epsilon) / (2+epsilon))), (size_t) 1);

  auto sortD = sequence<uintE>(n, [&](size_t i) {
    return i;
  });
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).getOutDegree(); });
  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto get_deg =
      [&](uintE& p) -> uintE { return D[p]; };
  for (size_t start = 0; start < n; start += ns) {
    // sort vertices in GA by degree, from start to n
    //integer_sort_inplace(sortD.slice(start, n), get_deg);
    radix::parallelIntegerSort(sortD.begin() + start, n - start, get_deg);
    uintE deg_max = D[sortD[std::min(ns + start, n)]];
    
    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      //if (D[v] >= deg_max) {
        D[v] -= edgesRemoved;
      //  return wrap(v, D[v]);
      //}
      return Maybe<std::tuple<uintE, uintE> >();
    };
    auto active =
        vertexSubset(n, std::min(ns + start, n) - start, sortD.begin() + start);
    auto moved = em.template edgeMapCount_sparse<uintE>(active, apply_f);

    moved.del();
  }
  auto ret = sequence<uintE>::no_init(n);
  parallel_for (0,n,[&] (size_t j) { ret[sortD[j]] = j; });
  return ret;
}

// keep track of induced subgraph as you go up -- store edge lists
// this would be P alpha k space; P k if we edidn't have induced subgraph, but longer to do k way intersect instead of 2 way intersect (k factor in work)

// TODO Pnk space without ordering; induced subgraphs have to be stored in hash tables

// can preinitialize k arrays of size n for each processor, and reuse when you do mem allocations -- check
// which processor is doing allocation and get the space assoc w/that processor
template <class Graph>
inline size_t KCliqueDir_rec(Graph& DG, size_t k_idx, size_t k, sequence<uintE> base) {
  // intersect outneighbors of verts in base
  auto lst_intersect = kintersect(DG, base, k_idx); // TODO hash table?, vectors, induced subgraph
  size_t num_intersect = lst_intersect.size();

  if (k_idx == k) {
    return num_intersect;
  }
  //auto counts = sequence<size_t>(num_intersect);
  size_t total_ct = 0;
  // then, for each v in the intersection
  for (size_t i=0; i < num_intersect; ++i) {
    base[k_idx] = lst_intersect[i]; //if par here, must duplicate base
    total_ct += KCliqueDir_rec(DG, k_idx+1, k, base);
  }
  // TODO leave this for now, unroll loop + reuse base in that -- write a program to generate loop unrolling
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}

template <class Graph>
inline size_t KCliqueDir(Graph& DG, size_t k) {
  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
  //for (size_t i = 0; i < DG.n; ++i) {
    auto base_idxs = sequence<uintE>::no_init(k);
    base_idxs[0] = i;
    tots[i] = KCliqueDir_rec(DG, 1, k, base_idxs);
  });
  return pbbslib::reduce_add(tots);
}
// STACK allocatee when < 1000 or something
// threshold alpha
// pbbslib list allocator


// base must have space for k if count_only = false
template <class Graph, class S, class F, class G>
inline size_t KCliqueIndDir_rec(Graph& DG, size_t k_idx, size_t k, S induced, F lstintersect,
  sequence<uintE> base, G g_f, bool count_only = true) {
  size_t num_intersect = induced.size();

  if (k_idx == k) {
    g_f(base);
    return num_intersect;
  }
  //auto counts = sequence<size_t>(num_intersect);
  // then, for each v in the intersection

  // optimization if counting and not listing
  if (k_idx + 1 == k && count_only) {
    auto counts = sequence<size_t>::no_init(num_intersect);
    parallel_for (0, num_intersect, [&] (size_t i) {
      counts[i] = std::get<1>(lstintersect(DG, induced[i], induced, false));
    });
    return pbbslib::reduce_add(counts);
  }

  size_t total_ct = 0;
  for (size_t i=0; i < num_intersect; ++i) {
    if (!count_only) base[k_idx] = induced[i];
    auto new_induced = std::get<0>(lstintersect(DG, induced[i], induced, true));
    total_ct += KCliqueIndDir_rec(DG, k_idx+1, k, new_induced, lstintersect, base, g_f, count_only);
  }
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}

template <class Graph, class F, class G>
inline size_t KCliqueIndDir(Graph& DG, size_t k, F lstintersect, G g_f, bool count_only = true) {
  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
    sequence<uintE> base = sequence<uintE>();
    if (!count_only) {
      base = sequence<uintE>::no_init(k);
      base[0] = i;
    }
    auto induced = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(i).getOutNeighbors()), DG.get_vertex(i).getOutDegree());
    tots[i] = KCliqueIndDir_rec(DG, 1, k, induced, lstintersect, base, g_f, count_only);
  });
  return pbbslib::reduce_add(tots);
}


// induced
// generated
// -i 0 (simple gbbs intersect), -i 1 (set intersect), -i 2 (simd intersect)

// todo approx work and do some kind of break in gen if too much
template <class Graph>
inline size_t KClique(Graph& GA, size_t k, double epsilon=0.001,
  bool induced = true, bool gen = true, long inter = 0, long order = 0) {
  using W = typename Graph::weight_type;
  assert (k >= 1);
  if (k == 1) return GA.n;
  else if (k == 2) return GA.m;

  sequence<uintE> rank;
  timer t_rank; t_rank.start();
  if (order == 0) rank = AppKCore(GA, epsilon);
  else if (order == 1) rank = DensestAppDegenOrder(GA, epsilon, false);
  else rank = DensestAppDegenOrder(GA, epsilon, true);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filter_graph(GA, pack_predicate);
  auto nop_f = [&] (sequence<uintE> b) {return;};

  timer t; t.start();
  size_t count = 0;
  if (!induced && !gen) count = KCliqueDir(DG, k-1);
  else if (induced && !gen) {
    if (inter == 0) count = KCliqueIndDir(DG, k-1, lstintersect_par_struct{}, nop_f, true);
    //else if (inter == 1) count = KCliqueIndDir(DG, k-1, lstintersect_set_struct{}, nop_f, true);
    else if (inter == 2) count = KCliqueIndDir(DG, k-1, lstintersect_vec_struct{}, nop_f, true);
  }
  else if (induced && gen) {
    if (inter == 0) count = KCliqueIndGenDir(DG, k-1, lstintersect_par_struct{}, nop_f, true);
    //else if (inter == 1) count = KCliqueIndGenDir(DG, k-1, lstintersect_set_struct{}, nop_f, true);
    else if (inter == 2) count = KCliqueIndGenDir(DG, k-1, lstintersect_vec_struct{}, nop_f, true);
  }
  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";
  return count;
}

//**********************************************************************GENERATED

// TODO keep array of size order alpha per processor???

template <class Graph, class F, class G>
inline size_t KCliqueIndGenDir(Graph& DG, size_t k, F lstintersect, G g_f, bool count_only = true) {
using induced_alloc = list_allocator<uintE[INDUCED_STACK_THR]>; 
 induced_alloc::init();
 switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(a).getOutNeighbors()), DG.get_vertex(a).getOutDegree());
 auto sizea = induceda.size();
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 if (count_only) {
 auto tupleb = lstintersect(DG, induceda[b], induceda, false);
 sizeb = std::get<1>(tupleb);
 } else {
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(DG, induceda[b], induceda, true, *ptr_inducedb);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto ptr_base = induced_alloc::alloc();
 auto base = sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 for (size_t xx = 0; xx < sizeb; xx++) {
 base[2] = inducedb[xx];
 g_f(base);
 }
 base.to_array(); induced_alloc::free(ptr_base);
 inducedb.to_array(); induced_alloc::free(ptr_inducedb);
 }
 storeb[b] = sizeb;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 storeb.to_array(); induced_alloc::free(ptr_storeb);
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(a).getOutNeighbors()), DG.get_vertex(a).getOutDegree());
 auto sizea = induceda.size();
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(DG, induceda[b], induceda, true, *ptr_inducedb);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto ptr_storec = induced_alloc::alloc();
 auto storec = sequence<uintE>(*ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 if (count_only) {
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, false);
 sizec = std::get<1>(tuplec);
 } else {
 auto ptr_inducedc = induced_alloc::alloc();
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true, *ptr_inducedc);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto ptr_base = induced_alloc::alloc();
 auto base = sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 for (size_t xx = 0; xx < sizec; xx++) {
 base[3] = inducedc[xx];
 g_f(base);
 }
 base.to_array(); induced_alloc::free(ptr_base);
 inducedc.to_array(); induced_alloc::free(ptr_inducedc);
 }
 storec[c] = sizec;
 });
 inducedb.to_array(); induced_alloc::free(ptr_inducedb);
 storeb[b] = pbbslib::reduce_add(storec);
 storec.to_array(); induced_alloc::free(ptr_storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 storeb.to_array(); induced_alloc::free(ptr_storeb);
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea);
 break; }
 default:
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(a).getOutNeighbors()), DG.get_vertex(a).getOutDegree());
 auto sizea = induceda.size();
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(DG, induceda[b], induceda, true, *ptr_inducedb);
 auto inducedb = std::get<0>(tupleb);
 size_t sizeb = std::get<1>(tupleb);
 auto ptr_storec = induced_alloc::alloc();
 auto storec = sequence<uintE>(*ptr_storec, sizeb);
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto ptr_inducedc = induced_alloc::alloc();
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true, *ptr_inducedc);
 auto inducedc = std::get<0>(tuplec);
 auto base = sequence<uintE>();
 if (!count_only) {
 auto ptr_base = induced_alloc::alloc();
 base = sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 storec[c] = KCliqueIndDir_rec(DG, 3, k, inducedc, lstintersect, base, g_f, count_only);
 base.to_array(); induced_alloc::free(ptr_base); }
 else storec[c] = KCliqueIndDir_rec(DG, 3, k, inducedc, lstintersect, base, g_f, count_only);
 });
 inducedb.to_array(); induced_alloc::free(ptr_inducedb);
 storeb[b] = pbbslib::reduce_add(storec);
 storec.to_array(); induced_alloc::free(ptr_storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 storeb.to_array(); induced_alloc::free(ptr_storeb);
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea); 
 }
}

//size_t temp(graph<vertex<W>>& DG, size_t k) {
//  auto storea = pbbslib::no_init<size_t>(DG.n);
//  parallel_for (0, DG.n, [&] (size_t a) {
//    auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
//    auto storeb = pbbslib::no_init<size_t>(induceda.size());
//    parallel_for (0, induceda.size(), [&] (size_t b) {
//      auto inducedb = lstintersect(DG, induceda[b], induceda);
      // for k = 3, store size of induced here; to be returned
//      storeb[b] = inducedb.size();
      // if we were recursing further, we would invoke the recursive version here
//    });
    // reduce storeb; store that in storea
//    storea[a] = pbbslib::reduce_add(storeb);
//  });
//  return pbbslib::reduce_add(storea);
//}

// Ideas:
// hash table for outvert of G -- would make it technically work-efficient, to do intersection
// where should we parallelize? first level only? through? be careful of space usage
// would approx kcore be faster if we used buckets instead? instead of the sort?