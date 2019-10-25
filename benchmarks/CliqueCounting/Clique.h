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
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

#define SIMD_STATE 4

template <class Graph>
inline uintE* rankNodes(Graph& G, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n);
  sequence<uintE> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).getOutDegree() < G.get_vertex(v).getOutDegree();
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  debug(t.reportTotal("Rank time"););
  return r;
}

// induced_space must have: num_induced, .del()
template <class IN, class Graph, class I, class F, class G, class H>
inline size_t KCliqueDir_rec(Graph& DG, size_t k_idx, size_t k, I& induced_space,
  F intersect_op, G intersect_op_type, sequence<uintE>& base, H base_op, bool count_only = true) {
  size_t num_induced = induced_space.num_induced;
  if (num_induced == 0) return 0;

  if (k_idx == k) {
    base_op(base);
    return num_induced;
  }

  // optimization if counting and not listing
  if (k_idx + 1 == k && count_only) {
    if (induced_space.full_flag) return induced_space.num_edges;
    auto counts = sequence<size_t>::no_init(num_induced);
    parallel_for (0, num_induced, [&] (size_t i) {
      auto new_induced_space = IN();
      counts[i] = intersect_op(DG, k_idx, k, i, induced_space, intersect_op_type, base, count_only, false, new_induced_space);
    });
    return pbbslib::reduce_add(counts);
  }

  size_t total_ct = 0;
  for (size_t i=0; i < num_induced; ++i) {
    auto new_induced_space = IN();
    size_t new_num_induced = intersect_op(DG, k_idx, k, i, induced_space, intersect_op_type, base, count_only, true, new_induced_space);
  
    if (new_num_induced > 0) { // >= k - k_idx
      total_ct += KCliqueDir_rec<IN>(DG, k_idx + 1, k, new_induced_space, intersect_op, intersect_op_type, base, base_op, count_only);
      new_induced_space.del();
    }
  }

  return total_ct;
}

template <class IN, class I, class Graph, class F, class G, class H>
inline size_t KCliqueDir(Graph& DG, size_t k, F intersect_op, G intersect_op_type, H base_op, bool count_only = true) {
  IN::init();

  auto tots = sequence<size_t>::no_init(DG.n);
  //size_t stepSize=15;
  //for (size_t i=0; i < DG.n; ++i) {
  //for(size_t step = 0; step < (DG.n+stepSize-1)/stepSize; step++) {
   //parallel_for(step*stepSize, std::min((size_t) (step+1)*stepSize,(size_t) DG.n), [&] (size_t i) {
  parallel_for (0, DG.n,[&] (size_t i) {
    if (DG.get_vertex(i).getOutDegree() == 0) {
      tots[i] = 0;
    } else {
      sequence<uintE> base = sequence<uintE>();
      if (!count_only) {
        base = sequence<uintE>::no_init(k);
        base[0] = i;
      }
      I induced_space = I(DG, k, i);
      if (induced_space.num_induced == 0) tots[i] = 0;
      else tots[i] = KCliqueDir_rec<IN>(DG, 1, k, induced_space, intersect_op, intersect_op_type, base, base_op, count_only);
      induced_space.del();
    }
  });
  //}

  IN::finish();
  return pbbslib::reduce_add(tots);
}

//template <class IN, class I, class Graph, class F, class H>
//size_t KCliqueDirGen(Graph& DG, size_t k, F intersect_op, H base_op, bool count_only);

template <class Graph>
void assert_induced_stack_thr(Graph& DG, size_t k = 1) {
  auto idxs = sequence<size_t>::no_init(DG.n);
  parallel_for (0,DG.n,[&] (size_t i) { idxs[i] = DG.get_vertex(i).getOutDegree(); });
  auto base_deg_f = [&](size_t i, size_t j) -> size_t {
    return idxs[i] > idxs[j] ? idxs[i] : idxs[j];
  };
  size_t max_deg = pbbslib::reduce(idxs, pbbslib::make_monoid(base_deg_f, 0));
  assert (max_deg*k <= INDUCED_STACK_THR);
}

template <class Graph, class F>
size_t assemble_induced_KCliqueDir(Graph& DG, size_t k, F inter_use, long subspace_type, bool count_only) {
  auto nop_f = [] (sequence<uintE> b) {return;};
  //auto lstintersect = [&](auto& DGA, size_t k_idx, size_t i, auto& induced_space, sequence<uintE>& base, bool to_save, auto& new_induced_space) {return lstintersect_induced(DGA, k_idx, k-1, i, induced_space, inter_use, base, count_only, to_save, new_induced_space);};
  auto lstintersect = lstintersect_induced_struct2{};
  if (subspace_type == 0) return KCliqueDir<InducedSpace_dyn, InducedSpace_dyn>(DG, k-1, lstintersect, inter_use, nop_f, count_only);
  else if (subspace_type == 1) {
    assert_induced_stack_thr(DG);
    return KCliqueDir<InducedSpace_alloc, InducedSpace_dyn>(DG, k-1, lstintersect, inter_use, nop_f, count_only);
  }
  else if (subspace_type == 2) {
    assert_induced_stack_thr(DG);
    return KCliqueDir<InducedSpace_stack, InducedSpace_dyn>(DG, k-1, lstintersect, inter_use, nop_f, count_only);
  }
  else if (subspace_type == 3) {
    assert_induced_stack_thr(DG, k);
    return KCliqueDir<InducedSpace_rec, InducedSpace_stack_setup>(DG, k-1, lstintersect, inter_use, nop_f, count_only);
  }
}
/*
template <class Graph, class F>
size_t assemble_induced_KCliqueDirGen(Graph& DG, size_t k, F inter_use, long subspace_type, bool count_only) {
  auto nop_f = [] (sequence<uintE> b) {return;};
  auto lstintersect = lstintersect_induced_struct{k, inter_use, count_only};
  if (subspace_type == 0) return KCliqueDirGen<InducedSpace_dyn, InducedSpace_dyn>(DG, k-1, lstintersect, nop_f, count_only);
    else if (subspace_type == 1) {
      assert_induced_stack_thr(DG);
      return KCliqueDirGen<InducedSpace_alloc, InducedSpace_dyn>(DG, k-1, lstintersect, nop_f, count_only);
    }
    else {
      assert_induced_stack_thr(DG);
      return KCliqueDirGen<InducedSpace_stack, InducedSpace_dyn>(DG, k-1, lstintersect, nop_f, count_only);
    }
}*/

// induced
// generated
// -i 0 (simple gbbs intersect), -i 2 (simd intersect), -i 1 (graph set inter)
// -o 0 (goodrich), 1 (barnboimelkin approx), 2 (barenboimelkin exact)

// todo approx work and do some kind of break in gen if too much
// TODO get rid of duplicates in edge lists????
template <class Graph>
inline size_t KClique(Graph& GA, size_t k, long order_type = 0, double epsilon = 0.1, 
bool gen_type = true, long space_type = 0, long subspace_type = 0, long inter_type = 0) {
  using W = typename Graph::weight_type;
  assert (k >= 1);
  if (k == 1) return GA.n;
  else if (k == 2) return GA.m;

  sequence<uintE> rank;
  timer t_rank; t_rank.start();
  if (order_type == 0) rank = goodrichpszona_degen::DegeneracyOrder_intsort(GA, epsilon);
  else if (order_type == 1) rank = barenboimelkin_degen::DegeneracyOrder(GA, epsilon);
  else if (order_type == 2) {
    rank = sequence<uintE>(GA.n, [&](size_t i) { return i; });
    auto kcore = KCore(GA);
    auto get_core = [&](uintE& p) -> uintE { return kcore[p]; };
    integer_sort_inplace(rank.slice(), get_core);
  }
  else if (order_type == 3) rank = pbbslib::make_sequence(rankNodes(GA, GA.n), GA.n);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= k-1 && GA.get_vertex(v).getOutDegree() >= k-1;
  };
  auto DG = filter_graph(GA, pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;

  // Done preprocessing
  timer t; t.start();
  size_t count = 0;
  bool count_only = true;

  if (!gen_type && space_type == 0) {
    if (inter_type == 0){
      count = assemble_induced_KCliqueDir(DG, k, lstintersect_par_struct{}, subspace_type, count_only);
    }
    else if (inter_type == 1) {
      assert (DG.n < INT_MAX);
      count = assemble_induced_KCliqueDir(DG, k, lstintersect_set_struct{}, subspace_type, count_only);
    }
    else if (inter_type == 2){
      count = assemble_induced_KCliqueDir(DG, k, lstintersect_vec_struct{}, subspace_type, count_only);
    }
    else {
      count = assemble_induced_KCliqueDir(DG, k, lstintersect_simple_struct{}, subspace_type, count_only);
    }
  }
  /*else if (!gen_type && space_type == 1) {
    auto nop_f = [] (sequence<uintE> b) {return;};
    auto inter_use = lstintersect_vec_struct{};
    auto lstintersect = [&](auto& DGA, size_t k_idx, size_t i, auto& induced_space, sequence<uintE>& base, bool to_save, auto& new_induced_space) {return lstintersect_full(DGA, k_idx, k-1, i, induced_space, inter_use, base, count_only, to_save, new_induced_space);};
    if (subspace_type == 0) count = KCliqueDir<FullSpace_bool_dyn, FullSpace_bool_dyn>(DG, k-1, lstintersect, nop_f, count_only);
    //else if (subspace_type == 1) count = KCliqueDir<FullSpace_csv_hash_dyn, FullSpace_csv_hash_dyn>(DG, k-1, lstintersect, nop_f, count_only);
    else count = KCliqueDir<FullSpace_csv_dyn, FullSpace_csv_dyn>(DG, k-1, lstintersect, nop_f, count_only);
  }
  else if (gen_type && space_type == 0) {
    if (inter_type == 0){
      count = assemble_induced_KCliqueDirGen(DG, k, lstintersect_par_struct{}, subspace_type, count_only);
    }
    else if (inter_type == 1) {
      assert (DG.n < INT_MAX);
      count = assemble_induced_KCliqueDirGen(DG, k, lstintersect_set_struct{}, subspace_type, count_only);
    }
    else {
      count = assemble_induced_KCliqueDirGen(DG, k, lstintersect_vec_struct{}, subspace_type, count_only);
    }
  }*/
  /*if (!induced && !gen) count = KCliqueIndDir_alloc(DG, k-1, lstintersect_par_struct{}, nop_f, true); //count = KCliqueDir(DG, k-1);
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
  }*/
  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";
  return count;
}


/*
template <class IN, class I, class Graph, class F, class H>
inline size_t KCliqueDirGen(Graph& DG, size_t k, F intersect_op, H base_op, bool count_only) {
IN::init();
 sequence<uintE> base = sequence<uintE>();
 if (!count_only) base = sequence<uintE>::no_init(k);
 switch (k) {
 case 2:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = I(DG, k, a);
 auto sizea = induceda.num_induced;
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto inducedb = IN();
 if (count_only) {
 sizeb = intersect_op(DG, 1, b, induceda, base, false, inducedb);
 } else {
 sizeb = intersect_op(DG, 1, b, induceda, base, true, inducedb);
 if (sizeb >= k - 1) {
 for (size_t xx = 0; xx < sizeb; xx++) {
 base[2] = inducedb.induced[xx];
 base_op(base);
 }
 }
 inducedb.del();
 }
 storeb[b] = sizeb;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 IN::finish();
 return pbbslib::reduce_add(storea);
 break; }
 case 3:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = I(DG, k, a);
 auto sizea = induceda.num_induced;
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto inducedb = IN();
 sizeb = intersect_op(DG, 1, b, induceda, base, true, inducedb);
 if (sizeb >= k - 1) {
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto inducedc = IN();
 if (count_only) {
 sizec = intersect_op(DG, 2, c, inducedb, base, false, inducedc);
 } else {
 sizec = intersect_op(DG, 2, c, inducedb, base, true, inducedc);
 if (sizec >= k - 2) {
 for (size_t xx = 0; xx < sizec; xx++) {
 base[3] = inducedc.induced[xx];
 base_op(base);
 }
 }
 inducedc.del();
 }
 storec[c] = sizec;
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 inducedb.del();
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 IN::finish();
 return pbbslib::reduce_add(storea);
 break; }
 case 4:  {
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = I(DG, k, a);
 auto sizea = induceda.num_induced;
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 size_t sizeb = 0;
 auto inducedb = IN();
 sizeb = intersect_op(DG, 1, b, induceda, base, true, inducedb);
 if (sizeb >= k - 1) {
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 size_t sizec = 0;
 auto inducedc = IN();
 sizec = intersect_op(DG, 2, c, inducedb, base, true, inducedc);
 if (sizec >= k - 2) {
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 size_t sized = 0;
 auto inducedd = IN();
 if (count_only) {
 sized = intersect_op(DG, 3, d, inducedc, base, false, inducedd);
 } else {
 sized = intersect_op(DG, 3, d, inducedc, base, true, inducedd);
 if (sized >= k - 3) {
 for (size_t xx = 0; xx < sized; xx++) {
 base[4] = inducedd.induced[xx];
 base_op(base);
 }
 }
 inducedd.del();
 }
 stored[d] = sized;
 });
 storec[c] = pbbslib::reduce_add(stored);} else storec[c] = 0;
 inducedc.del();
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 inducedb.del();
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 IN::finish();
 return pbbslib::reduce_add(storea);
 break; }
 default:
 auto storea = sequence<size_t>::no_init(DG.n);
 parallel_for (0, DG.n, [&] (size_t a) {
 auto induceda = I(DG, k, a);
 auto sizea = induceda.num_induced;
 if (sizea >= k) {
 auto storeb = sequence<size_t>::no_init(sizea);
 parallel_for (0, sizea, [&] (size_t b) {
 auto inducedb = IN();
 size_t sizeb = intersect_op(DG, 1, b, induceda, base, true, inducedb);
 if (sizeb >= k - 1) {
 auto storec = sequence<size_t>::no_init(sizeb);
 parallel_for (0, sizeb, [&] (size_t c) {
 auto inducedc = IN();
 size_t sizec = intersect_op(DG, 2, c, inducedb, base, true, inducedc);
 if (sizec >= k - 2) {
 auto stored = sequence<size_t>::no_init(sizec);
 parallel_for (0, sizec, [&] (size_t d) {
 auto inducedd = IN();
 size_t sized = intersect_op(DG, 3, d, inducedc, base, true, inducedd);
 if (sized >= k - 3) {
 stored[d] = KCliqueDir_rec<IN>(DG, 4, k, inducedd, intersect_op, base, base_op, count_only);} else stored[d] = 0;
 });
 storec[c] = pbbslib::reduce_add(stored);} else storec[c] = 0;
 });
 storeb[b] = pbbslib::reduce_add(storec);} else storeb[b] = 0;
 });
 storea[a] = pbbslib::reduce_add(storeb);
 } else storea[a] = 0;
 });
 IN::finish();
 return pbbslib::reduce_add(storea); 
 }
}
*/





/*

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


//GENERATED

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
*/