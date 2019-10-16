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

#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"
#include "intersect.h"
//#include "radix_wrapper.h"

#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"

#define SIMD_STATE 4

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


// base must have space for k if count_only = false
template <class Graph, class F, class G>
inline size_t KCliqueIndDir_rec(Graph& DG, size_t k_idx, size_t k, uintE* induced, size_t num_intersect, F lstintersect_sub,
  sequence<uintE> base, G g_f, bool count_only = true) {
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
      auto tup = lstintersect(lstintersect_sub, DG, induced[i], induced, num_intersect, false);
      counts[i] = std::get<1>(tup);
    });
    return pbbslib::reduce_add(counts);
  }

  size_t total_ct = 0;
  for (size_t i=0; i < num_intersect; ++i) {
    if (!count_only) base[k_idx] = induced[i];
    auto new_induced_tup = lstintersect(lstintersect_sub, DG, induced[i], induced, num_intersect, true);
    auto new_induced = std::get<0>(new_induced_tup);
    auto new_induced_size = std::get<1>(new_induced_tup);
    if (new_induced_size > 0) {
      total_ct += KCliqueIndDir_rec(DG, k_idx+1, k, new_induced, new_induced_size, lstintersect_sub, base, g_f, count_only);
      pbbs::delete_array<uintE>(new_induced, new_induced_size);
    }
  }
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}
//TODO del array in ind gen dyn

template <class Graph, class F, class G>
inline size_t KCliqueIndDir(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
    if (DG.get_vertex(i).getOutDegree() == 0) {
      tots[i] = 0;
    } else {
    sequence<uintE> base = sequence<uintE>();
    if (!count_only) {
      base = sequence<uintE>::no_init(k);
      base[0] = i;
    }
    tots[i] = KCliqueIndDir_rec(DG, 1, k, (uintE*)(DG.get_vertex(i).getOutNeighbors()), DG.get_vertex(i).getOutDegree(), lstintersect_sub, base, g_f, count_only);
    }
  });
  return pbbslib::reduce_add(tots);
}

template <class Graph, class F, class G>
inline size_t KCliqueIndDir_stack_rec(Graph& DG, size_t k_idx, size_t k, uintE* induced, size_t induced_size, F lstintersect_sub,
  sequence<uintE> base, G g_f, bool count_only = true) {
  if (k_idx == k) {
    g_f(base);
    return induced_size;
  }
  // then, for each v in the intersection
  uintE induced_ptr[INDUCED_STACK_THR];

  // optimization if counting and not listing
  if (k_idx + 1 == k && count_only) {
    auto counts = sequence<size_t>::no_init(induced_size);
    parallel_for (0, induced_size, [&] (size_t i) {
      counts[i] = std::get<1>(lstintersect(lstintersect_sub, DG, induced[i], induced, induced_size, false, induced_ptr));
    });
    return pbbslib::reduce_add(counts);
  }

  size_t total_ct = 0;
  for (size_t i=0; i < induced_size; ++i) {
    if (!count_only) base[k_idx] = induced[i];
    auto new_induced_size = std::get<1>(lstintersect(lstintersect_sub, DG, induced[i], induced, induced_size, true, induced_ptr));
    //if (new_induced_size > 0) assert((induced_ptr + (k_idx * INDUCED_STACK_THR))[0] < DG.n);
    if (new_induced_size > 0) total_ct += KCliqueIndDir_stack_rec(DG, k_idx+1, k, induced_ptr, new_induced_size, lstintersect_sub, base, g_f, count_only);
  }
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}

template <class Graph, class F, class G>
inline size_t KCliqueIndDir_stack(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
  auto idxs = sequence<size_t>::no_init(DG.n);
  parallel_for (0,DG.n,[&] (size_t i) { idxs[i] = DG.get_vertex(i).getOutDegree(); });
  auto base_deg_f = [&](size_t i, size_t j) -> size_t {
    return idxs[i] > idxs[j] ? idxs[i] : idxs[j];
  };
  size_t max_deg = pbbslib::reduce(idxs, pbbslib::make_monoid(base_deg_f, 0));
  assert (max_deg <= INDUCED_STACK_THR);

  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
    if (DG.get_vertex(i).getOutDegree() == 0) {
      tots[i] = 0;
    } else{
    sequence<uintE> base = sequence<uintE>();
    if (!count_only) {
      base = sequence<uintE>::no_init(k);
      base[0] = i;
    }
    uintE induced_ptr[INDUCED_STACK_THR];
    for (size_t j=0; j < DG.get_vertex(i).getOutDegree(); ++j) {
      induced_ptr[j] = DG.get_vertex(i).getOutNeighbor(j);
    }
    //auto induced = pbbslib::make_sequence<uintE>(induced_ptr, DG.get_vertex(i).getOutDegree());
    tots[i] = KCliqueIndDir_stack_rec(DG, 1, k, induced_ptr, DG.get_vertex(i).getOutDegree(), lstintersect_sub, base, g_f, count_only);
    }
  });
  return pbbslib::reduce_add(tots);
}

template <class Graph, class F, class G>
inline size_t KCliqueIndDir_alloc_rec(Graph& DG, size_t k_idx, size_t k, uintE* induced, size_t granularity, size_t induced_size, F lstintersect_sub,
  sequence<uintE> base, G g_f, bool count_only = true) {
  if (k_idx == k) {
    g_f(base);
    return induced_size;
  }
  // then, for each v in the intersection
  auto new_induced = induced + granularity;

  // optimization if counting and not listing
  if (k_idx + 1 == k && count_only) {
    auto counts = sequence<size_t>::no_init(induced_size);
    parallel_for (0, induced_size, [&] (size_t i) {
      counts[i] = std::get<1>(lstintersect(lstintersect_sub, DG, induced[i], induced, induced_size, false, new_induced));
    });
    return pbbslib::reduce_add(counts);
  }

  size_t total_ct = 0;
  for (size_t i=0; i < induced_size; ++i) {
    if (!count_only) base[k_idx] = induced[i];
    auto new_induced_size = std::get<1>(lstintersect(lstintersect_sub, DG, induced[i], induced, induced_size, true, new_induced));
    //if (new_induced_size > 0) assert((induced_ptr + (k_idx * INDUCED_STACK_THR))[0] < DG.n);
    if (new_induced_size >= k - k_idx) total_ct += KCliqueIndDir_alloc_rec(DG, k_idx+1, k, new_induced, granularity, new_induced_size, lstintersect_sub, base, g_f, count_only);
  }
  //auto count_seq = pbbslib::make_sequence<size_t>(counts, active_size);
  //size_t count = pbbslib::reduce_add(count_seq);

  return total_ct;
}

// TODO prune all vert with deg < k
template <class Graph, class F, class G>
inline size_t KCliqueIndDir_alloc(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
  auto idxs = sequence<size_t>::no_init(DG.n);
  parallel_for (0,DG.n,[&] (size_t i) { idxs[i] = DG.get_vertex(i).getOutDegree(); });
  auto base_deg_f = [&](size_t i, size_t j) -> size_t {
    return idxs[i] > idxs[j] ? idxs[i] : idxs[j];
  };
  size_t max_deg = pbbslib::reduce(idxs, pbbslib::make_monoid(base_deg_f, 0));
  assert (k * max_deg <= INDUCED_STACK_THR);

  using induced_alloc = list_allocator<uintE[INDUCED_STACK_THR]>; 
  induced_alloc::init();

  // TODO divide work -- statically or by estimating prefix sum stuff
  auto tots = sequence<size_t>::no_init(DG.n);
  parallel_for (0, DG.n,[&] (size_t i) {
    if (DG.get_vertex(i).getOutDegree() < k) {
      tots[i] = 0;
    } else{
      auto induced_ptr = induced_alloc::alloc();
      sequence<uintE> base = sequence<uintE>();
      if (!count_only) {
        base = sequence<uintE>::no_init(k);
        base[0] = i;
      }
      for (size_t j=0; j < DG.get_vertex(i).getOutDegree(); ++j) {
        (*induced_ptr)[j] = DG.get_vertex(i).getOutNeighbor(j);
      }
      auto induced_deg = DG.get_vertex(i).getOutDegree();
      tots[i] = KCliqueIndDir_alloc_rec(DG, 1, k, *induced_ptr, induced_deg, induced_deg, lstintersect_sub, base, g_f, count_only);
      induced_alloc::free(induced_ptr);
    }
  });

  induced_alloc::finish();
  return pbbslib::reduce_add(tots);
}


// induced
// generated
// -i 0 (simple gbbs intersect), -i 2 (simd intersect), -i 1 (graph set inter)
// -o 0 (goodrich), 1 (barnboimelkin approx), 2 (barenboimelkin exact)

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
  if (order == 0) rank = goodrichpszona_degen::DegeneracyOrder(GA, epsilon);
  else if (order == 1) rank = barenboimelkin_degen::DegeneracyOrder(GA, epsilon, false);
  else rank = barenboimelkin_degen::DegeneracyOrder(GA, epsilon, true);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && !(GA.get_vertex(u).getOutDegree() < k || GA.get_vertex(v).getOutDegree()  < k);
  };
  auto DG = filter_graph(GA, pack_predicate);
  auto nop_f = [&] (sequence<uintE> b) {return;};

  /*auto k_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return !(DG_prev.get_vertex(u).getOutDegree() < k || DG_prev.get_vertex(v).getOutDegree()  < k);
  };
  auto DG = filter_graph(DG_prev, k_predicate);*/

  timer t; t.start();
  size_t count = 0;
  if (!induced && !gen) count = KCliqueIndDir_alloc(DG, k-1, lstintersect_par_struct{}, nop_f, true); //count = KCliqueDir(DG, k-1);
  else if (induced && !gen) {
    if (inter == 0) count = KCliqueIndDir(DG, k-1, lstintersect_par_struct{}, nop_f, true);
    else if (inter == 1) {
      assert (DG.n < INT_MAX);
      count = KCliqueIndDir(DG, k-1, lstintersect_set_struct{}, nop_f, true);
    }
    else if (inter == 2) count = KCliqueIndDir(DG, k-1, lstintersect_vec_struct{}, nop_f, true);
  }
  else if (induced && gen) {
    if (inter == 0) count = KCliqueIndGenDir(DG, k-1, lstintersect_par_struct{}, nop_f, true);
    else if (inter == 1) {
      assert (DG.n < INT_MAX);
      count = KCliqueIndGenDir(DG, k-1, lstintersect_set_struct{}, nop_f, true);
    }
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
inline size_t KCliqueIndGenDir(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
  auto base_idxs = sequence<size_t>::no_init(DG.n);
  parallel_for (0,DG.n,[&] (size_t i) { base_idxs[i] = DG.get_vertex(i).getOutDegree(); });
  auto base_deg_f = [&](size_t i, size_t j) -> size_t {
    return base_idxs[i] > base_idxs[j] ? base_idxs[i] : base_idxs[j];
  };
  size_t max_deg = pbbslib::reduce(base_idxs, pbbslib::make_monoid(base_deg_f, 0));
  if (max_deg <= INDUCED_STACK_THR) return KCliqueIndGenDir_stack(DG, k, lstintersect_sub, g_f, count_only);
  return KCliqueIndGenDir_dyn(DG, k, lstintersect_sub, g_f, count_only);
}

// using stack
template <class Graph, class F, class G>
inline size_t KCliqueIndGenDir_stack(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 uintE ptr_storeb[INDUCED_STACK_THR];
 auto storeb = pbbslib::make_sequence<uintE>(ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 if (count_only) {
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, false);
 sizeb = std::get<1>(tupleb);
 } else {
 uintE ptr_inducedb[INDUCED_STACK_THR];
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, ptr_inducedb);
 sizeb = std::get<1>(tupleb);
 auto inducedb = std::get<0>(tupleb);
 if (sizeb >= k - 1) {
 uintE ptr_base[INDUCED_STACK_THR];
 auto base = pbbslib::make_sequence<uintE>(ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 for (size_t xx = 0; xx < sizeb; xx++) {
 base[2] = inducedb[xx];
 g_f(base);
 }
 }
 }
 storeb[b] = sizeb;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 uintE ptr_storeb[INDUCED_STACK_THR];
 auto storeb = pbbslib::make_sequence<uintE>(ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 uintE ptr_inducedb[INDUCED_STACK_THR];
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, ptr_inducedb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto inducedb = std::get<0>(tupleb);
 uintE ptr_storec[INDUCED_STACK_THR];
 auto storec = pbbslib::make_sequence<uintE>(ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 if (count_only) {
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, false);
 sizec = std::get<1>(tuplec);
 } else {
 uintE ptr_inducedc[INDUCED_STACK_THR];
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true, ptr_inducedc);
 sizec = std::get<1>(tuplec);
 auto inducedc = std::get<0>(tuplec);
 if (sizec >= k - 2) {
 uintE ptr_base[INDUCED_STACK_THR];
 auto base = pbbslib::make_sequence<uintE>(ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 for (size_t xx = 0; xx < sizec; xx++) {
 base[3] = inducedc[xx];
 g_f(base);
 }
 }
 }
 storec[c] = sizec;
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 4:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 uintE ptr_storeb[INDUCED_STACK_THR];
 auto storeb = pbbslib::make_sequence<uintE>(ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 uintE ptr_inducedb[INDUCED_STACK_THR];
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, ptr_inducedb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto inducedb = std::get<0>(tupleb);
 uintE ptr_storec[INDUCED_STACK_THR];
 auto storec = pbbslib::make_sequence<uintE>(ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 uintE ptr_inducedc[INDUCED_STACK_THR];
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true, ptr_inducedc);
 sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto inducedc = std::get<0>(tuplec);
 uintE ptr_stored[INDUCED_STACK_THR];
 auto stored = pbbslib::make_sequence<uintE>(ptr_stored, sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 if (count_only) {
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, false);
 sized = std::get<1>(tupled);
 } else {
 uintE ptr_inducedd[INDUCED_STACK_THR];
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, true, ptr_inducedd);
 sized = std::get<1>(tupled);
 auto inducedd = std::get<0>(tupled);
 if (sized >= k - 3) {
 uintE ptr_base[INDUCED_STACK_THR];
 auto base = pbbslib::make_sequence<uintE>(ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 for (size_t xx = 0; xx < sized; xx++) {
 base[4] = inducedd[xx];
 g_f(base);
 }
 }
 }
 stored[d] = sized;
 });
 storec[c] = pbbslib::reduce_add(stored);} else storec[c] = 0;
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea);
 break; }
 default:
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 uintE ptr_storeb[INDUCED_STACK_THR];
 auto storeb = pbbslib::make_sequence<uintE>(ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 uintE ptr_inducedb[INDUCED_STACK_THR];
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, ptr_inducedb);
 size_t sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto inducedb = std::get<0>(tupleb);
 uintE ptr_storec[INDUCED_STACK_THR];
 auto storec = pbbslib::make_sequence<uintE>(ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 uintE ptr_inducedc[INDUCED_STACK_THR];
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true, ptr_inducedc);
 size_t sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto inducedc = std::get<0>(tuplec);
 uintE ptr_stored[INDUCED_STACK_THR];
 auto stored = pbbslib::make_sequence<uintE>(ptr_stored, sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 uintE ptr_inducedd[INDUCED_STACK_THR];
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, true, ptr_inducedd);
 auto inducedd = std::get<0>(tupled);
 size_t sized = std::get<1>(tupled);
 if (sized >= k - 3) {
 auto base = sequence<uintE>();
 if (!count_only) {
 uintE ptr_base[INDUCED_STACK_THR];
 base = sequence<uintE>(ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 stored[d] = KCliqueIndDir_rec(DG, 4, k, inducedd, sized, lstintersect_sub, base, g_f, count_only);
 } else stored[d] = KCliqueIndDir_rec(DG, 4, k, inducedd, sized, lstintersect_sub, base, g_f, count_only);} else stored[d] = 0;
 });
 storec[c] = pbbslib::reduce_add(stored);} else storec[c] = 0;
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea); 
 }
}

// using list allocator
template <class Graph, class F, class G>
inline size_t KCliqueIndGenDir_alloc(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
using induced_alloc = list_allocator<uintE[INDUCED_STACK_THR]>; 
 induced_alloc::init();
 switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = pbbslib::make_sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 if (count_only) {
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, false);
 sizeb = std::get<1>(tupleb);
 } else {
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, *ptr_inducedb);
 sizeb = std::get<1>(tupleb);
 auto inducedb = std::get<0>(tupleb);
 if (sizeb >= k - 1) {
 auto ptr_base = induced_alloc::alloc();
 auto base = pbbslib::make_sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 for (size_t xx = 0; xx < sizeb; xx++) {
 base[2] = inducedb[xx];
 g_f(base);
 }
 induced_alloc::free(ptr_base);
 }
 induced_alloc::free(ptr_inducedb);
 }
 storeb[b] = sizeb;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 induced_alloc::free(ptr_storeb);
 } else storea[a] = 0;
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = pbbslib::make_sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, *ptr_inducedb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto inducedb = std::get<0>(tupleb);
 auto ptr_storec = induced_alloc::alloc();
 auto storec = pbbslib::make_sequence<uintE>(*ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 if (count_only) {
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, false);
 sizec = std::get<1>(tuplec);
 } else {
 auto ptr_inducedc = induced_alloc::alloc();
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true, *ptr_inducedc);
 sizec = std::get<1>(tuplec);
 auto inducedc = std::get<0>(tuplec);
 if (sizec >= k - 2) {
 auto ptr_base = induced_alloc::alloc();
 auto base = pbbslib::make_sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 for (size_t xx = 0; xx < sizec; xx++) {
 base[3] = inducedc[xx];
 g_f(base);
 }
 induced_alloc::free(ptr_base);
 }
 induced_alloc::free(ptr_inducedc);
 }
 storec[c] = sizec;
 });
 storeb[b] = pbbslib::reduce_add(storec);
 induced_alloc::free(ptr_storec);} else storeb[b] = 0;
 induced_alloc::free(ptr_inducedb);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 induced_alloc::free(ptr_storeb);
 } else storea[a] = 0;
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea);
 break; }
 case 4:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = pbbslib::make_sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, *ptr_inducedb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto inducedb = std::get<0>(tupleb);
 auto ptr_storec = induced_alloc::alloc();
 auto storec = pbbslib::make_sequence<uintE>(*ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto ptr_inducedc = induced_alloc::alloc();
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true, *ptr_inducedc);
 sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto inducedc = std::get<0>(tuplec);
 auto ptr_stored = induced_alloc::alloc();
 auto stored = pbbslib::make_sequence<uintE>(*ptr_stored, sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 if (count_only) {
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, false);
 sized = std::get<1>(tupled);
 } else {
 auto ptr_inducedd = induced_alloc::alloc();
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, true, *ptr_inducedd);
 sized = std::get<1>(tupled);
 auto inducedd = std::get<0>(tupled);
 if (sized >= k - 3) {
 auto ptr_base = induced_alloc::alloc();
 auto base = pbbslib::make_sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 for (size_t xx = 0; xx < sized; xx++) {
 base[4] = inducedd[xx];
 g_f(base);
 }
 induced_alloc::free(ptr_base);
 }
 induced_alloc::free(ptr_inducedd);
 }
 stored[d] = sized;
 });
 storec[c] = pbbslib::reduce_add(stored);
 induced_alloc::free(ptr_stored);} else storec[c] = 0;
 induced_alloc::free(ptr_inducedc);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 induced_alloc::free(ptr_storec);} else storeb[b] = 0;
 induced_alloc::free(ptr_inducedb);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 induced_alloc::free(ptr_storeb);
 } else storea[a] = 0;
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea);
 break; }
 default:
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto ptr_storeb = induced_alloc::alloc();
 auto storeb = pbbslib::make_sequence<uintE>(*ptr_storeb, sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 auto ptr_inducedb = induced_alloc::alloc();
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true, *ptr_inducedb);
 size_t sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto inducedb = std::get<0>(tupleb);
 auto ptr_storec = induced_alloc::alloc();
 auto storec = pbbslib::make_sequence<uintE>(*ptr_storec, sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 auto ptr_inducedc = induced_alloc::alloc();
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true, *ptr_inducedc);
 size_t sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto inducedc = std::get<0>(tuplec);
 auto ptr_stored = induced_alloc::alloc();
 auto stored = pbbslib::make_sequence<uintE>(*ptr_stored, sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 auto ptr_inducedd = induced_alloc::alloc();
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, true, *ptr_inducedd);
 auto inducedd = std::get<0>(tupled);
 size_t sized = std::get<1>(tupled);
 if (sized >= k - 3) {
 auto base = sequence<uintE>();
 if (!count_only) {
 auto ptr_base = induced_alloc::alloc();
 base = sequence<uintE>(*ptr_base, k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 stored[d] = KCliqueIndDir_rec(DG, 4, k, inducedd, sized, lstintersect_sub, base, g_f, count_only);
 induced_alloc::free(ptr_base); }
 else stored[d] = KCliqueIndDir_rec(DG, 4, k, inducedd, sized, lstintersect_sub, base, g_f, count_only);} else stored[d] = 0;
 });
 storec[c] = pbbslib::reduce_add(stored);
 induced_alloc::free(ptr_stored);} else storec[c] = 0;
 induced_alloc::free(ptr_inducedc);
 });
 storeb[b] = pbbslib::reduce_add(storec);
 induced_alloc::free(ptr_storec);} else storeb[b] = 0;
 induced_alloc::free(ptr_inducedb);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 induced_alloc::free(ptr_storeb);
 } else storea[a] = 0;
 });
 induced_alloc::finish();
 return pbbslib::reduce_add(storea); 
 }
}

template <class Graph, class F, class G>
inline size_t KCliqueIndGenDir_dyn(Graph& DG, size_t k, F lstintersect_sub, G g_f, bool count_only = true) {
switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 if (count_only) {
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, false);
 sizeb = std::get<1>(tupleb);
 } else {
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 for (size_t xx = 0; xx < sizeb; xx++) {
 base[2] = inducedb[xx];
 g_f(base);
 }
 }
 if (inducedb) pbbs::delete_array<uintE>(inducedb, sizeb);
 }
 storeb[b] = sizeb;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 if (count_only) {
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, false);
 sizec = std::get<1>(tuplec);
 } else {
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 for (size_t xx = 0; xx < sizec; xx++) {
 base[3] = inducedc[xx];
 g_f(base);
 }
 }
 if (inducedc) pbbs::delete_array<uintE>(inducedc, sizec);
 }
 storec[c] = sizec;
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 if (inducedb) pbbs::delete_array<uintE>(inducedb, sizeb);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea);
 break; }
 case 4:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true);
 auto inducedb = std::get<0>(tupleb);
 sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true);
 auto inducedc = std::get<0>(tuplec);
 sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 if (count_only) {
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, false);
 sized = std::get<1>(tupled);
 } else {
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, true);
 auto inducedd = std::get<0>(tupled);
 sized = std::get<1>(tupled);
 if (sized >= k - 3) {
 auto base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 for (size_t xx = 0; xx < sized; xx++) {
 base[4] = inducedd[xx];
 g_f(base);
 }
 }
 if (inducedd) pbbs::delete_array<uintE>(inducedd, sized);
 }
 stored[d] = sized;
 });
 storec[c] = pbbslib::reduce_add(stored);} else storec[c] = 0;
 if (inducedc) pbbs::delete_array<uintE>(inducedc, sizec);
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 if (inducedb) pbbs::delete_array<uintE>(inducedb, sizeb);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 return pbbslib::reduce_add(storea);
 break; }
 default:
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = (uintE*)(DG.get_vertex(a).getOutNeighbors());
 auto sizea = DG.get_vertex(a).getOutDegree();
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 auto tupleb = lstintersect(lstintersect_sub, DG, induceda[b], induceda, sizea, true);
 auto inducedb = std::get<0>(tupleb);
 size_t sizeb = std::get<1>(tupleb);
 if (sizeb >= k - 1) {
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 auto tuplec = lstintersect(lstintersect_sub, DG, inducedb[c], inducedb, sizeb, true);
 auto inducedc = std::get<0>(tuplec);
 size_t sizec = std::get<1>(tuplec);
 if (sizec >= k - 2) {
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 auto tupled = lstintersect(lstintersect_sub, DG, inducedc[d], inducedc, sizec, true);
 auto inducedd = std::get<0>(tupled);
 size_t sized = std::get<1>(tupled);
 if (sized >= k - 3) {
 auto base = sequence<uintE>();
 if (!count_only) {
 base = sequence<uintE>::no_init(k);
 base[0] = a;
 base[1] = induceda[b];
 base[2] = inducedb[c];
 base[3] = inducedc[d];
 }
 stored[d] = KCliqueIndDir_rec(DG, 4, k, inducedd, sized, lstintersect_sub, base, g_f, count_only);} else stored[d] = 0;
 });
 storec[c] = pbbslib::reduce_add(stored);} else storec[c] = 0;
 if (inducedc) pbbs::delete_array<uintE>(inducedc, sizec);
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 if (inducedb) pbbs::delete_array<uintE>(inducedb, sizeb);
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
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