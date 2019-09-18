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
#include "simdinter/include/intersection.h"

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

// Goodrich (2+epsilon) approx for degeneracy ordering where epsilon > 0
// Returns vertice sorted in degeneracy order
template<template <typename W> class vertex, class W>
inline sequence<uintE> AppKCore(graph<vertex<W>>& GA, double epsilon=0.1) {
  const size_t n = GA.n;
  const size_t ns = std::max((size_t) (ceil((n*epsilon) / (2+epsilon))), (size_t) 1);

  auto sortD = sequence<uintE>(n, [&](size_t i) {
    return i;
  });
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });
  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto get_deg =
      [&](uintE& p) -> uintE { return D[p]; };
  for (size_t start = 0; start < n; start += ns) {
    // sort vertices in GA by degree, from start to n
    integer_sort_inplace(sortD.slice(start, n), get_deg);
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


template <template <typename W> class vertex, class W>
inline sequence<uintE> kintersect(graph<vertex<W>>& DG, sequence<uintE> base, size_t num) {
  if (num == 1) {
    uintT deg = DG.V[base[0]].getOutDegree();
    uintE* ngh = (uintE*)(DG.V[base[0]].getOutNeighbors());
    return pbbslib::make_sequence<uintE>(ngh, deg);
  }

  // find base index with min outdeg
  auto base_idxs = sequence<size_t>::no_init(num);
  parallel_for (0,num,[&] (size_t i) { base_idxs[i] = i; });
  auto base_deg_f = [&](size_t i, size_t j) -> size_t {
    return DG.V[base[i]].getOutDegree() < DG.V[base[j]].getOutDegree() ? i : j;
  };
  size_t min_base = pbbslib::reduce(base_idxs, pbbslib::make_monoid(base_deg_f, 0));
  size_t min_deg = DG.V[base[min_base]].getOutDegree();
  auto min_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.V[base[min_base]].getOutNeighbors()), min_deg);

  // set up array to mark where min seq intersects with other outneighbors in base
  auto marks = sequence<bool>(min_deg * (num - 1), false);
  // mark intersections
  parallel_for (0,num,[&] (size_t i) {
    if (i != min_base) {
    size_t idx = i > min_base ? i - 1 : i;
    auto merge_f = [&] (uintE ngh) {
      size_t j = pbbslib::binary_search(min_seq, ngh, std::less<size_t>());
      marks[min_deg * idx + j] = true;
    };
    auto a_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.V[base[i]].getOutNeighbors()), DG.V[base[i]].getOutDegree());
    intersection::merge(min_seq, a_seq, merge_f);
    }
  });

  // merge marks to determine which intersects all (can make more efficient here)
  auto merge_seq = sequence<bool>(min_deg, true);
  parallel_for (0,min_deg,[&] (size_t i) {
    for (size_t j=0; j < num - 1; ++j) {
      merge_seq[i] = merge_seq[i] && marks[min_deg * j + i];
    }
  });
  // pack out intersected vertices
  auto out = pbbslib::pack(min_seq, merge_seq);
  return out;
  // want to intersect V[base[i]] outneighbors
  // figure out which one represents the smallest
  // intersect each against the smallest
  // any element in the smallest that has num-1 marks is in our merged list

}

// keep track of induced subgraph as you go up -- store edge lists -- write  another func
// this would be P alpha k space; P k if we edidn't have induced subgraph, but longer to do k way intersect instead of 2 way intersect (k factor in work)

// TODO Pnk space without ordering; induced subgraphs have to be stored in hash tables

// can preinitialize k arrays of size n for each processor, and reuse when you do mem allocations -- check
// which processor is doing allocation and get the space assoc w/that processor
template <template <class W> class vertex, class W>
inline size_t KCliqueDir_rec(graph<vertex<W>>& DG, size_t k_idx, size_t k, sequence<uintE> base) {
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

template <template <class W> class vertex, class W>
inline size_t KCliqueDir(graph<vertex<W>>& DG, size_t k) {
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

struct lstintersect_par_struct {
  template <template <typename W> class vertex, class W, class S>
  sequence<uintE> operator()(graph<vertex<W>>& DG, uintE vtx, S induced) const {
    return lstintersect_par(DG, vtx, induced);
  }
};

template <template <typename W> class vertex, class W, class S>
inline sequence<uintE> lstintersect_par(graph<vertex<W>>& DG, uintE vtx, S induced) {
  auto vtx_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.V[vtx].getOutNeighbors()), DG.V[vtx].getOutDegree());
  size_t index = 0;
  auto out = sequence<uintE>::no_init(std::min(induced.size(), vtx_seq.size()));
  auto merge_f = [&] (uintE ngh) {
    out[index] = ngh;
    index++;
  };
  intersection::merge(vtx_seq, induced, merge_f);
  out.shrink(index);
  return out;
  // want to intersect V[base[i]] outneighbors
  // figure out which one represents the smallest
  // intersect each against the smallest
  // any element in the smallest that has num-1 marks is in our merged list
}

struct lstintersect_set_struct {
  template <template <typename W> class vertex, class W, class S>
  std::vector<uintE> operator()(graph<vertex<W>>& DG, uintE vtx, S induced) const {
    return lstintersect_set(DG, vtx, induced);
  }
};

// make sure set intersection is stable
template <template <typename W> class vertex, class W, class S>
inline std::vector<uintE> lstintersect_set(graph<vertex<W>>& DG, uintE vtx, S induced) {
  auto vtx_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.V[vtx].getOutNeighbors()), DG.V[vtx].getOutDegree());
  std::vector<uintE> out;
  std::set_intersection(induced.begin(), induced.end(), vtx_seq.begin(), vtx_seq.end(), std::back_inserter(out));
  return out;
}

struct lstintersect_vec_struct {
  template <template <typename W> class vertex, class W, class S>
  sequence<uintE> operator()(graph<vertex<W>>& DG, uintE vtx, S induced) const {
    return lstintersect_vec(DG, vtx, induced);
  }
};

// TODO radix sort in place
template <template <typename W> class vertex, class W, class S>
inline sequence<uintE> lstintersect_vec(graph<vertex<W>>& DG, uintE vtx, S induced) {
  auto vtx_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.V[vtx].getOutNeighbors()), DG.V[vtx].getOutDegree());
  auto out = sequence<uintE>::no_init(std::min(induced.size(), vtx_seq.size()));
  SIMDCompressionLib::intersectionfunction inter = SIMDCompressionLib::IntersectionFactory::getFromName("simd");
  size_t out_size = inter(induced.begin(), induced.size(), vtx_seq.begin(), vtx_seq.size(), out.begin());
  out.shrink(out_size);
  return out;
}

template <template <class W> class vertex, class W, class S, class F>
inline size_t KCliqueIndDir_rec(graph<vertex<W>>& DG, size_t k_idx, size_t k, S induced, F lstintersect) {
  size_t num_intersect = induced.size();
  // intersect outneighbors of verts in base
  //auto lst_intersect = lstintersect(DG, base, k_idx); // TODO hash table?, vectors, induced subgraph
  //size_t num_intersect = lst_intersect.size();
  if (k_idx == k) {
    return num_intersect;
  }
  //auto counts = sequence<size_t>(num_intersect);
  size_t total_ct = 0;
  // then, for each v in the intersection
  for (size_t i=0; i < num_intersect; ++i) {
    auto new_induced = lstintersect(DG, induced[i], induced);
    total_ct += KCliqueIndDir_rec(DG, k_idx+1, k, new_induced, lstintersect);
  }
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}

template <template <class W> class vertex, class W, class F>
inline size_t KCliqueIndDir(graph<vertex<W>>& DG, size_t k, F lstintersect) {
  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
  //for (size_t i = 0; i < DG.n; ++i) {
    auto induced = pbbslib::make_sequence<uintE>((uintE*)(DG.V[i].getOutNeighbors()), DG.V[i].getOutDegree());
    tots[i] = KCliqueIndDir_rec(DG, 1, k, induced, lstintersect);
  });
  return pbbslib::reduce_add(tots);
}

template <template <class W> class vertex, class W>
inline size_t KClique(graph<vertex<W>>& GA, size_t k, double epsilon=0.1, bool induced = true, bool gen = true, long inter = 0) {
  assert (k >= 1);
  if (k == 1) return GA.n;
  else if (k == 2) return GA.m;

  timer t_rank; t_rank.start();
  auto rank = AppKCore(GA, epsilon);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filter_graph<vertex, W>(GA, pack_predicate);

  timer t; t.start();
  size_t count = 0;
  if (!induced && !gen) count = KCliqueDir(DG, k-1);
  else if (induced && !gen) {
    if (inter == 0) count = KCliqueIndDir(DG, k-1, lstintersect_par_struct{});
    else if (inter == 1) count = KCliqueIndDir(DG, k-1, lstintersect_set_struct{});
    else if (inter == 2) count = KCliqueIndDir(DG, k-1, lstintersect_vec_struct{});
  }
  else if (induced && gen) {
    if (inter == 0) count = KCliqueIndGenDir(DG, k-1, lstintersect_par_struct{});
    else if (inter == 1) count = KCliqueIndGenDir(DG, k-1, lstintersect_set_struct{});
    else if (inter == 2) count = KCliqueIndGenDir(DG, k-1, lstintersect_vec_struct{});
  }
  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";
  return count;
}

// TODO For listing vs counting? have a func passed in?

//**********************************************************************GENERATED

template <template <class W> class vertex, class W, class F>
inline size_t KCliqueIndGenDir(graph<vertex<W>>& DG, size_t k, F lstintersect) {
switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 storeb[b] = inducedb.size();
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 storec[c] = inducedc.size();
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 4:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 stored[d] = inducedd.size();
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 5:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 auto storee = sequence<size_t>::no_init(inducedd.size());
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 auto inducede = lstintersect(DG, inducedd[e], inducedd);
 storee[e] = inducede.size();
 });
 stored[d] = pbbslib::reduce_add(storee);
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 6:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 auto storee = sequence<size_t>::no_init(inducedd.size());
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 auto inducede = lstintersect(DG, inducedd[e], inducedd);
 auto storef = sequence<size_t>::no_init(inducede.size());
 parallel_for (0, inducede.size(), [&] (size_t f) {
 auto inducedf = lstintersect(DG, inducede[f], inducede);
 storef[f] = inducedf.size();
 });
 storee[e] = pbbslib::reduce_add(storef);
 });
 stored[d] = pbbslib::reduce_add(storee);
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 7:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 auto storee = sequence<size_t>::no_init(inducedd.size());
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 auto inducede = lstintersect(DG, inducedd[e], inducedd);
 auto storef = sequence<size_t>::no_init(inducede.size());
 parallel_for (0, inducede.size(), [&] (size_t f) {
 auto inducedf = lstintersect(DG, inducede[f], inducede);
 auto storeg = sequence<size_t>::no_init(inducedf.size());
 parallel_for (0, inducedf.size(), [&] (size_t g) {
 auto inducedg = lstintersect(DG, inducedf[g], inducedf);
 storeg[g] = inducedg.size();
 });
 storef[f] = pbbslib::reduce_add(storeg);
 });
 storee[e] = pbbslib::reduce_add(storef);
 });
 stored[d] = pbbslib::reduce_add(storee);
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 8:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 auto storee = sequence<size_t>::no_init(inducedd.size());
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 auto inducede = lstintersect(DG, inducedd[e], inducedd);
 auto storef = sequence<size_t>::no_init(inducede.size());
 parallel_for (0, inducede.size(), [&] (size_t f) {
 auto inducedf = lstintersect(DG, inducede[f], inducede);
 auto storeg = sequence<size_t>::no_init(inducedf.size());
 parallel_for (0, inducedf.size(), [&] (size_t g) {
 auto inducedg = lstintersect(DG, inducedf[g], inducedf);
 auto storeh = sequence<size_t>::no_init(inducedg.size());
 parallel_for (0, inducedg.size(), [&] (size_t h) {
 auto inducedh = lstintersect(DG, inducedg[h], inducedg);
 storeh[h] = inducedh.size();
 });
 storeg[g] = pbbslib::reduce_add(storeh);
 });
 storef[f] = pbbslib::reduce_add(storeg);
 });
 storee[e] = pbbslib::reduce_add(storef);
 });
 stored[d] = pbbslib::reduce_add(storee);
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 9:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 auto storee = sequence<size_t>::no_init(inducedd.size());
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 auto inducede = lstintersect(DG, inducedd[e], inducedd);
 auto storef = sequence<size_t>::no_init(inducede.size());
 parallel_for (0, inducede.size(), [&] (size_t f) {
 auto inducedf = lstintersect(DG, inducede[f], inducede);
 auto storeg = sequence<size_t>::no_init(inducedf.size());
 parallel_for (0, inducedf.size(), [&] (size_t g) {
 auto inducedg = lstintersect(DG, inducedf[g], inducedf);
 auto storeh = sequence<size_t>::no_init(inducedg.size());
 parallel_for (0, inducedg.size(), [&] (size_t h) {
 auto inducedh = lstintersect(DG, inducedg[h], inducedg);
 auto storei = sequence<size_t>::no_init(inducedh.size());
 parallel_for (0, inducedh.size(), [&] (size_t i) {
 auto inducedi = lstintersect(DG, inducedh[i], inducedh);
 storei[i] = inducedi.size();
 });
 storeh[h] = pbbslib::reduce_add(storei);
 });
 storeg[g] = pbbslib::reduce_add(storeh);
 });
 storef[f] = pbbslib::reduce_add(storeg);
 });
 storee[e] = pbbslib::reduce_add(storef);
 });
 stored[d] = pbbslib::reduce_add(storee);
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 default:
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto storeb = sequence<size_t>::no_init(induceda.size());
 parallel_for (0, induceda.size(), [&] (size_t b) {
 auto inducedb = lstintersect(DG, induceda[b], induceda);
 auto storec = sequence<size_t>::no_init(inducedb.size());
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 auto inducedc = lstintersect(DG, inducedb[c], inducedb);
 auto stored = sequence<size_t>::no_init(inducedc.size());
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 auto inducedd = lstintersect(DG, inducedc[d], inducedc);
 auto storee = sequence<size_t>::no_init(inducedd.size());
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 auto inducede = lstintersect(DG, inducedd[e], inducedd);
 auto storef = sequence<size_t>::no_init(inducede.size());
 parallel_for (0, inducede.size(), [&] (size_t f) {
 auto inducedf = lstintersect(DG, inducede[f], inducede);
 auto storeg = sequence<size_t>::no_init(inducedf.size());
 parallel_for (0, inducedf.size(), [&] (size_t g) {
 auto inducedg = lstintersect(DG, inducedf[g], inducedf);
 auto storeh = sequence<size_t>::no_init(inducedg.size());
 parallel_for (0, inducedg.size(), [&] (size_t h) {
 auto inducedh = lstintersect(DG, inducedg[h], inducedg);
 auto storei = sequence<size_t>::no_init(inducedh.size());
 parallel_for (0, inducedh.size(), [&] (size_t i) {
 auto inducedi = lstintersect(DG, inducedh[i], inducedh);
 storei[i] = KCliqueIndDir_rec(DG, 9, k, inducedi, lstintersect);
 });
 storeh[h] = pbbslib::reduce_add(storei);
 });
 storeg[g] = pbbslib::reduce_add(storeh);
 });
 storef[f] = pbbslib::reduce_add(storeg);
 });
 storee[e] = pbbslib::reduce_add(storef);
 });
 stored[d] = pbbslib::reduce_add(storee);
 });
 storec[c] = pbbslib::reduce_add(stored);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
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

// Note: gbbs main runner is broken