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
#include "intersect.h"

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

// keep track of induced subgraph as you go up -- store edge lists
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

// base must have space for k if count_only = false
template <template <class W> class vertex, class W, class S, class F, class G>
inline size_t KCliqueIndDir_rec(graph<vertex<W>>& DG, size_t k_idx, size_t k, S induced, F lstintersect,
  sequence<uintE> base, G g, bool count_only = true) {
  size_t num_intersect = induced.size();

  if (k_idx == k) {
  //  g(base);
    return num_intersect;
  }
  //auto counts = sequence<size_t>(num_intersect);
  size_t total_ct = 0;
  // then, for each v in the intersection

  // optimization if counting and not listing
  if (k_idx + 1 == k && count_only) {
    for (size_t i=0; i < num_intersect; ++i) {
      total_ct += std::get<1>(lstintersect(DG, induced[i], induced, false));
    }
    return total_ct;
  }

  for (size_t i=0; i < num_intersect; ++i) {
    if (!count_only) base[k_idx] = induced[i];
    auto new_induced = std::get<0>(lstintersect(DG, induced[i], induced, true));
    total_ct += KCliqueIndDir_rec(DG, k_idx+1, k, new_induced, lstintersect, base, g, count_only);
  }
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}

template <template <class W> class vertex, class W, class F, class G>
inline size_t KCliqueIndDir(graph<vertex<W>>& DG, size_t k, F lstintersect, G g, bool count_only = true) {
  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
    sequence<uintE> base = sequence<uintE>();
    if (!count_only) {
      base = sequence<uintE>::no_init(k);
      base[0] = i;
    }
    auto induced = pbbslib::make_sequence<uintE>((uintE*)(DG.V[i].getOutNeighbors()), DG.V[i].getOutDegree());
    tots[i] = KCliqueIndDir_rec(DG, 1, k, induced, lstintersect, base, g, count_only);
  });
  return pbbslib::reduce_add(tots);
}

// todo approx work and do some kind of break in gen if too much
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
    if (inter == 0) count = KCliqueIndDir(DG, k-1, lstintersect_par_struct{}, nop_f{}, true);
    else if (inter == 1) count = KCliqueIndDir(DG, k-1, lstintersect_set_struct{}, nop_f{}, true);
    else if (inter == 2) count = KCliqueIndDir(DG, k-1, lstintersect_vec_struct{}, nop_f{}, true);
  }
  else if (induced && gen) {
    if (inter == 0) count = KCliqueIndGenDir(DG, k-1, lstintersect_par_struct{}, nop_f{}, true);
    else if (inter == 1) count = KCliqueIndGenDir(DG, k-1, lstintersect_set_struct{}, nop_f{}, true);
    else if (inter == 2) count = KCliqueIndGenDir(DG, k-1, lstintersect_vec_struct{}, nop_f{}, true);
  }
  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";
  return count;
}

//**********************************************************************GENERATED

// TODO keep array of size order alpha per processor???

template <template <class W> class vertex, class W, class F, class G>
inline size_t KCliqueIndGenDir(graph<vertex<W>>& DG, size_t k, F lstintersect, G g, bool count_only = true) {
switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 if (count_only) {
 auto tupleb = lstintersect(DG, induceda[b], induceda, false);
 sizeb = std::get<1>(tupleb);
 } else {
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 for (size_t xx = 0; xx < sizeb; xx++) {
 base[2] = inducedb[xx];
 }
 }
 storeb[b] = sizeb;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = pbbslib::make_sequence<uintE>((uintE*)(DG.V[a].getOutNeighbors()), DG.V[a].getOutDegree());
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 if (count_only) {
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, false);
 sizec = std::get<1>(tuplec);
 } else {
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 for (size_t xx = 0; xx < sizec; xx++) {
 base[3] = inducedc[xx];
 }
 }
 storec[c] = sizec;
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
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 if (count_only) {
 auto tupled = lstintersect(DG, inducedc[d], inducedc, false);
 sized = std::get<1>(tupled);
 } else {
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 for (size_t xx = 0; xx < sized; xx++) {
 base[4] = inducedd[xx];
 }
 }
 stored[d] = sized;
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
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto storee = sequence<size_t>::no_init(sized);
 parallel_for (0, sized, [&] (size_t e) {
 size_t sizee = 0;
 if (count_only) {
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, false);
 sizee = std::get<1>(tuplee);
 } else {
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, true);
 auto inducede = std::get<0>(tuplee);
 sizee = std::get<1>(tuplee);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 base[4] = inducedd[e];
 for (size_t xx = 0; xx < sizee; xx++) {
 base[5] = inducede[xx];
 }
 }
 storee[e] = sizee;
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
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto storee = sequence<size_t>::no_init(sized);
 parallel_for (0, sized, [&] (size_t e) {
 size_t sizee = 0;
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, true);
 auto inducede = std::get<0>(tuplee);
 sizee = std::get<1>(tuplee);
 auto storef = sequence<size_t>::no_init(sizee);
 parallel_for (0, sizee, [&] (size_t f) {
 size_t sizef = 0;
 if (count_only) {
 auto tuplef = lstintersect(DG, inducede[f], inducede, false);
 sizef = std::get<1>(tuplef);
 } else {
 auto tuplef = lstintersect(DG, inducede[f], inducede, true);
 auto inducedf = std::get<0>(tuplef);
 sizef = std::get<1>(tuplef);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 base[4] = inducedd[e];
 base[5] = inducede[f];
 for (size_t xx = 0; xx < sizef; xx++) {
 base[6] = inducedf[xx];
 }
 }
 storef[f] = sizef;
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
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto storee = sequence<size_t>::no_init(sized);
 parallel_for (0, sized, [&] (size_t e) {
 size_t sizee = 0;
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, true);
 auto inducede = std::get<0>(tuplee);
 sizee = std::get<1>(tuplee);
 auto storef = sequence<size_t>::no_init(sizee);
 parallel_for (0, sizee, [&] (size_t f) {
 size_t sizef = 0;
 auto tuplef = lstintersect(DG, inducede[f], inducede, true);
 auto inducedf = std::get<0>(tuplef);
 sizef = std::get<1>(tuplef);
 auto storeg = sequence<size_t>::no_init(sizef);
 parallel_for (0, sizef, [&] (size_t g) {
 size_t sizeg = 0;
 if (count_only) {
 auto tupleg = lstintersect(DG, inducedf[g], inducedf, false);
 sizeg = std::get<1>(tupleg);
 } else {
 auto tupleg = lstintersect(DG, inducedf[g], inducedf, true);
 auto inducedg = std::get<0>(tupleg);
 sizeg = std::get<1>(tupleg);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 base[4] = inducedd[e];
 base[5] = inducede[f];
 base[6] = inducedf[g];
 for (size_t xx = 0; xx < sizeg; xx++) {
 base[7] = inducedg[xx];
 }
 }
 storeg[g] = sizeg;
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
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto storee = sequence<size_t>::no_init(sized);
 parallel_for (0, sized, [&] (size_t e) {
 size_t sizee = 0;
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, true);
 auto inducede = std::get<0>(tuplee);
 sizee = std::get<1>(tuplee);
 auto storef = sequence<size_t>::no_init(sizee);
 parallel_for (0, sizee, [&] (size_t f) {
 size_t sizef = 0;
 auto tuplef = lstintersect(DG, inducede[f], inducede, true);
 auto inducedf = std::get<0>(tuplef);
 sizef = std::get<1>(tuplef);
 auto storeg = sequence<size_t>::no_init(sizef);
 parallel_for (0, sizef, [&] (size_t g) {
 size_t sizeg = 0;
 auto tupleg = lstintersect(DG, inducedf[g], inducedf, true);
 auto inducedg = std::get<0>(tupleg);
 sizeg = std::get<1>(tupleg);
 auto storeh = sequence<size_t>::no_init(sizeg);
 parallel_for (0, sizeg, [&] (size_t h) {
 size_t sizeh = 0;
 if (count_only) {
 auto tupleh = lstintersect(DG, inducedg[h], inducedg, false);
 sizeh = std::get<1>(tupleh);
 } else {
 auto tupleh = lstintersect(DG, inducedg[h], inducedg, true);
 auto inducedh = std::get<0>(tupleh);
 sizeh = std::get<1>(tupleh);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 base[4] = inducedd[e];
 base[5] = inducede[f];
 base[6] = inducedf[g];
 base[7] = inducedg[h];
 for (size_t xx = 0; xx < sizeh; xx++) {
 base[8] = inducedh[xx];
 }
 }
 storeh[h] = sizeh;
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
 auto sizea = induceda.size();
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto storee = sequence<size_t>::no_init(sized);
 parallel_for (0, sized, [&] (size_t e) {
 size_t sizee = 0;
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, true);
 auto inducede = std::get<0>(tuplee);
 sizee = std::get<1>(tuplee);
 auto storef = sequence<size_t>::no_init(sizee);
 parallel_for (0, sizee, [&] (size_t f) {
 size_t sizef = 0;
 auto tuplef = lstintersect(DG, inducede[f], inducede, true);
 auto inducedf = std::get<0>(tuplef);
 sizef = std::get<1>(tuplef);
 auto storeg = sequence<size_t>::no_init(sizef);
 parallel_for (0, sizef, [&] (size_t g) {
 size_t sizeg = 0;
 auto tupleg = lstintersect(DG, inducedf[g], inducedf, true);
 auto inducedg = std::get<0>(tupleg);
 sizeg = std::get<1>(tupleg);
 auto storeh = sequence<size_t>::no_init(sizeg);
 parallel_for (0, sizeg, [&] (size_t h) {
 size_t sizeh = 0;
 auto tupleh = lstintersect(DG, inducedg[h], inducedg, true);
 auto inducedh = std::get<0>(tupleh);
 sizeh = std::get<1>(tupleh);
 auto storei = sequence<size_t>::no_init(sizeh);
 parallel_for (0, sizeh, [&] (size_t i) {
 size_t sizei = 0;
 if (count_only) {
 auto tuplei = lstintersect(DG, inducedh[i], inducedh, false);
 sizei = std::get<1>(tuplei);
 } else {
 auto tuplei = lstintersect(DG, inducedh[i], inducedh, true);
 auto inducedi = std::get<0>(tuplei);
 sizei = std::get<1>(tuplei);
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 base[4] = inducedd[e];
 base[5] = inducede[f];
 base[6] = inducedf[g];
 base[7] = inducedg[h];
 base[8] = inducedh[i];
 for (size_t xx = 0; xx < sizei; xx++) {
 base[9] = inducedi[xx];
 }
 }
 storei[i] = sizei;
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
 size_t sizeb = 0;
 auto tupleb = lstintersect(DG, induceda[b], induceda, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, inducedb.size(), [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(DG, inducedb[c], inducedb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, inducedc.size(), [&] (size_t d) {
 size_t sized = 0;
 auto tupled = lstintersect(DG, inducedc[d], inducedc, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 auto storee = sequence<size_t>::no_init(sized);
 parallel_for (0, inducedd.size(), [&] (size_t e) {
 size_t sizee = 0;
 auto tuplee = lstintersect(DG, inducedd[e], inducedd, true);
 auto inducede = std::get<0>(tuplee);
 sizee = std::get<1>(tuplee);
 auto storef = sequence<size_t>::no_init(sizee);
 parallel_for (0, inducede.size(), [&] (size_t f) {
 size_t sizef = 0;
 auto tuplef = lstintersect(DG, inducede[f], inducede, true);
 auto inducedf = std::get<0>(tuplef);
 sizef = std::get<1>(tuplef);
 auto storeg = sequence<size_t>::no_init(sizef);
 parallel_for (0, inducedf.size(), [&] (size_t g) {
 size_t sizeg = 0;
 auto tupleg = lstintersect(DG, inducedf[g], inducedf, true);
 auto inducedg = std::get<0>(tupleg);
 sizeg = std::get<1>(tupleg);
 auto storeh = sequence<size_t>::no_init(sizeg);
 parallel_for (0, inducedg.size(), [&] (size_t h) {
 size_t sizeh = 0;
 auto tupleh = lstintersect(DG, inducedg[h], inducedg, true);
 auto inducedh = std::get<0>(tupleh);
 sizeh = std::get<1>(tupleh);
 auto storei = sequence<size_t>::no_init(sizeh);
 parallel_for (0, inducedh.size(), [&] (size_t i) {
 size_t sizei = 0;
 auto tuplei = lstintersect(DG, inducedh[i], inducedh, true);
 auto inducedi = std::get<0>(tuplei);
 sizei = std::get<1>(tuplei);
 auto base = sequence<uintE>();
 if (!count_only) {
 base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 base[4] = inducedd[e];
 base[5] = inducede[f];
 base[6] = inducedf[g];
 base[7] = inducedg[h];
 base[8] = inducedh[i];
 }
 storei[i] = KCliqueIndDir_rec(DG, 9, k, inducedi, lstintersect, base, g, count_only);
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