#include <algorithm>
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "gbbs/gbbs.h"
#include "gbbs/graph.h"
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"

// ****************************************************************************************
// This is the ESCAPE five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_ESCAPE(Graph& GA, long order_type = 0, size_t block_size = 500000, double epsilon = 0.1) {
  using W = typename Graph::weight_type;
  timer tp; tp.start();
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
                { rank[order_to_vertex[i]] = GA.n-1-i; });
  } else if (order_type == 3) {
    rank = DegeneracyOrder(GA).to_seq();
    order_to_vertex = sequence<uintT>(GA.n, [&](size_t i){return 0;} );
    par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t v)
               { order_to_vertex[rank[v]] = v; });
  }
  debug(cout << "Rank abd Order done\n"; fflush(stdout););
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  //auto GDO = GA;
  debug(cout << "Relabel done\n"; fflush(stdout););

  auto out_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u < v;
  };
  auto OUTG = GDO.filterGraph(GDO, out_direction); // only keeps track of out edges

  auto in_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u > v;
  };
  auto ING = GDO.filterGraph(GDO, in_direction); // only keeps track of in edges

  debug(cout << "Filter done\n"; fflush(stdout););
  double ttp = tp.stop();
  GDO.del();

  timer tc; tc.start();
  ulong cycleCount = 0;
  auto U = sequence<uintE>(GA.n, [&](size_t s) { return 0; });
  auto distinct = sequence<uintE>(GA.n, [&](size_t s){ return 0; } ); 
  uintE num_distinct = 0;
  for (uintE i = 0; i < GA.n; i++) {
    ulong tmp = 0;
    //uintE* nghs = (uintE*) vi.getOutNeighbors();
    uintE viOutDegree  = OUTG.get_vertex(i).getOutDegree();
    uintE viInDegree  = ING.get_vertex(i).getOutDegree();
    if (viOutDegree + viInDegree == 0) continue; //return;

    // auto nghs_seq = sequence<uintE>(nghs, degree);
    uintE* outnghs_vi = (uintE*) OUTG.get_vertex(i).getOutNeighbors();
    uintE* innghs_vi = (uintE*) ING.get_vertex(i).getOutNeighbors();

    uintE vj, vk, vl;
    // i -> j -> k
    for (uintE j = 0; j < viOutDegree; j++){
      vj = outnghs_vi[j];
      uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
      for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
        vk = outnghs_vj[k];
        if (U[vk] == 0) { distinct[num_distinct] = vk; num_distinct++; }
        U[vk] += 1;
      }
    }

    // i <- j -> k or i <- j <- k
    for (uintE j = 0; j < viInDegree; j++){
        vj = innghs_vi[j];
        uintE* nghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
          vk = nghs_vj[k];
          if (U[vk] == 0) { distinct[num_distinct] = vk; num_distinct++; }
          U[vk] += 1;
        }
        nghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < ING.get_vertex(vj).getOutDegree();  k++) {
          vk = nghs_vj[k];
          if (U[vk] == 0) { distinct[num_distinct] = vk; num_distinct++; }
          U[vk] += 1;
        }
    }

    U[i] = 0;


    for (uintE j = 0; j < viInDegree; j++) {
      vj = innghs_vi[j];
      uintE* innghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
      uintE vjInDegree = ING.get_vertex(vj).getOutDegree();
      uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
      uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

      if (vjInDegree + vjOutDegree == 0) continue;

      //auto nghs_vj_seq = sequence<uintE>(nghs_vj, vj_degree);
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
          if (((index_leq_vk = _binary_search(outnghs_vi, custom_less_k, viOutDegree)) < viOutDegree )
              && (outnghs_vi[index_leq_vk] == vk))  tmp--;
          else if (((index_leq_vk = _binary_search(innghs_vi, custom_less_k, viInDegree)) < viInDegree )
              && (innghs_vi[index_leq_vk] == vk))  tmp--;

          if (((index_leq_vl = _binary_search(outnghs_vj, custom_less_l, vjOutDegree)) < vjOutDegree)
            && (outnghs_vj[index_leq_vl] == vl)) tmp--;
          else if  (((index_leq_vl = _binary_search(innghs_vj, custom_less_l, vjInDegree)) < vjInDegree)
            && (innghs_vj[index_leq_vl] == vl)) tmp--;
        }

      }
      //nghs_vj = nghs_vj_seq.to_array();
    }
    cycleCount += tmp;

    for (int z = 0; z < num_distinct; z++){
      U[distinct[z]] = 0;
      distinct[z] = 0;
    }
    num_distinct = 0; 
   //nghs = nghs_seq.to_array();
  }

  double ttc = tc.stop();
  debug(std::cout << "##### Actual counting: " << ttc << std::endl;);
  std::cout << "##### Result: " << cycleCount << "," << ttp << "," << 0 << "," << ttc << std::endl;
  //GDO.del();
  ING.del();
  OUTG.del();
  return cycleCount;
}

// ****************************************************************************************
// This is the ESCAPE five-cycle counting algorithm.
// ****************************************************************************************
template <class Graph>
inline ulong Count5Cycle_ESCAPE_par(Graph& GA, long order_type = 0, size_t block_size = 500000, double epsilon = 0.1) {
  timer tp; tp.start();
  using W = typename Graph::weight_type;
  sequence<uintE> rank ;
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
                { rank[order_to_vertex[i]] = GA.n-1-i; });
  } else if (order_type == 3) {
    rank = DegeneracyOrder(GA).to_seq();
    order_to_vertex = sequence<uintT>(GA.n, [&](size_t i){return 0;} );
    par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t v)
               { order_to_vertex[rank[v]] = v; });
  }
  debug(cout << "Rank abd Order done\n"; fflush(stdout););
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  debug(cout << "Relabel done\n"; fflush(stdout););

  auto out_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u < v;
  };
  auto OUTG = GDO.filterGraph(GDO, out_direction); // only keeps track of out edges
  auto in_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u > v;
  };
  auto ING = GDO.filterGraph(GDO, in_direction); // only keeps track of in edges


  debug(cout << "Filter done\n"; fflush(stdout););
  double ttp = tp.stop();

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

  // size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
  size_t work_per_block = total_work / n_blocks;
  n_blocks = (total_work/work_per_block) + 1;

  double tts = ts.stop();
  GDO.del();
  debug(std::cout << "##### Scheduling: " << tts << std::endl;);

  timer tc; tc.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  sequence<ulong> cycleCounts = sequence<ulong>(n_blocks * eltsPerCacheLine, [&](size_t s) { return 0; });

  auto run_intersection = [&](size_t start_ind, size_t end_ind, size_t block_index, U_FastReset* V) {
  
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto U = V->U;
      ulong tmp = 0;

      //uintE* nghs = (uintE*) vi.getOutNeighbors();
      uintE viOutDegree  = OUTG.get_vertex(i).getOutDegree();
      uintE viInDegree  = ING.get_vertex(i).getOutDegree();
      if (viOutDegree + viInDegree == 0) continue; //return;
      // auto nghs_seq = sequence<uintE>(nghs, degree);

      uintE* outnghs_vi = (uintE*) OUTG.get_vertex(i).getOutNeighbors();
      uintE* innghs_vi = (uintE*) ING.get_vertex(i).getOutNeighbors();

      uintE vj, vk, vl;
      // i -> j -> k
      for (uintE j = 0; j < viOutDegree; j++){
        vj = outnghs_vi[j];
        uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
          vk = outnghs_vj[k];
          if (U[vk] == 0) { V->distinct[V->num_distinct] = vk; V->num_distinct++; }
          U[vk] += 1;
        }
      }
      // i <- j -> k or i <- j <- k
      for (uintE j = 0; j < viInDegree; j++){
        vj = innghs_vi[j];
        uintE* nghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
          vk = nghs_vj[k];
          if (U[vk] == 0) { V->distinct[V->num_distinct] = vk; V->num_distinct++; }
          U[vk] += 1;
        }
        nghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < ING.get_vertex(vj).getOutDegree();  k++) {
          vk = nghs_vj[k];
          if (U[vk] == 0) { V->distinct[V->num_distinct] = vk; V->num_distinct++; }
          U[vk] += 1;
        }
      }
      U[i] = 0;

      for (uintE j = 0; j < viInDegree; j++) {
        vj = innghs_vi[j];
        uintE* innghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
        uintE vjInDegree = ING.get_vertex(vj).getOutDegree();
        uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

        if (vjInDegree + vjOutDegree == 0) continue;

        //auto nghs_vj_seq = sequence<uintE>(nghs_vj, vj_degree);
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
            if (((index_leq_vk = _binary_search(outnghs_vi, custom_less_k, viOutDegree)) < viOutDegree )
                && (outnghs_vi[index_leq_vk] == vk))  tmp--;
            else if (((index_leq_vk = _binary_search(innghs_vi, custom_less_k, viInDegree)) < viInDegree )
                && (innghs_vi[index_leq_vk] == vk))  tmp--;

            if (((index_leq_vl = _binary_search(outnghs_vj, custom_less_l, vjOutDegree)) < vjOutDegree)
              && (outnghs_vj[index_leq_vl] == vl)) tmp--;
            else if  (((index_leq_vl = _binary_search(innghs_vj, custom_less_l, vjInDegree)) < vjInDegree)
              && (innghs_vj[index_leq_vl] == vl)) tmp--;
          }
        }
        // nghs_vj = nghs_vj_seq.to_array();
      }
      cycleCounts[block_index * eltsPerCacheLine] += tmp;
      V->reset();
      // for (uintE j = 0; j < viOutDegree; j++){
      //   vj = outnghs_vi[j];
      //   uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
      //   uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

      //   for (uintE k = 0; k < vjOutDegree;  k++) {
      //     vk = outnghs_vj[k];
      //     U[vk] = 0;
      //   }

      // }
      // // i -> j <- k or i -> j -> k
      // for (uintE j = 0; j < viInDegree; j++){
      //   vj = innghs_vi[j];
      //   uintE* nghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
      //   for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
      //     U[nghs_vj[k]] = 0;
      //   }
      //   nghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
      //   for (uintE k = 0; k < ING.get_vertex(vj).getOutDegree();  k++) {
      //     U[nghs_vj[k]] = 0;
      //   }
      // }
      // nghs = nghs_seq.to_array();
    }
  };
  auto init_V = [&](U_FastReset* V){ V->alloc(GA.n); };
  auto finish_V = [&](U_FastReset* V){ delete V; };
  parallel_for_alloc<U_FastReset>(init_V, finish_V, 0, n_blocks, [&] (size_t i, U_FastReset* V) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = pbbslib::binary_search(parallel_work, start, less_fn);
    size_t end_ind = pbbslib::binary_search(parallel_work, end, less_fn);
    run_intersection(start_ind, end_ind, i, V);
  });
  ulong total = pbbslib::reduce_add(cycleCounts);

  double ttc = tc.stop();
  debug(std::cout << "##### Actual counting: " << ttc << std::endl;);
  std::cout << "##### Result: " << total << "," << ttp << "," << tts << "," << ttc << std::endl;
  ING.del();
  OUTG.del();
  return total;
}


template <class Graph>
inline ulong Count5Cycle_ESCAPE_no_scheduling(Graph& GA, long order_type = 0, size_t block_size = 500000, double epsilon = 0.1) {
  timer tp; tp.start();
  using W = typename Graph::weight_type;
  sequence<uintE> rank ;
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
                { rank[order_to_vertex[i]] = GA.n-1-i; });
  } else if (order_type == 3) {
    rank = DegeneracyOrder(GA).to_seq();
    order_to_vertex = sequence<uintT>(GA.n, [&](size_t i){return 0;} );
    par_for(0, GA.n, pbbslib::kSequentialForThreshold, [&] (size_t v)
               { order_to_vertex[rank[v]] = v; });
  }
  debug(cout << "Rank and Order done\n"; fflush(stdout););
  auto GDO = relabel_graph(GA, order_to_vertex); // graph by degree ordering

  debug(cout << "Relabel done\n"; fflush(stdout););

  auto out_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u < v;
  };
  auto OUTG = GDO.filterGraph(GDO, out_direction); // only keeps track of out edges
  auto in_direction = [&](const uintE& u, const uintE& v, const W& wgh) {
    return u > v;
  };
  auto ING = GDO.filterGraph(GDO, in_direction); // only keeps track of in edges


  debug(cout << "Filter done\n"; fflush(stdout););
  double ttp = tp.stop();
  GDO.del();
  debug(std::cout << "##### Scheduling: " << tts << std::endl;);

  timer tc; tc.start();
  const size_t eltsPerCacheLine = 64/sizeof(ulong);
  sequence<ulong> cycleCounts = sequence<ulong>(GA.n * eltsPerCacheLine, [&](size_t s) { return 0; });

  auto run_intersection = [&](size_t i, U_FastReset* V) {
      auto U = V->U;
      ulong tmp = 0;

      //uintE* nghs = (uintE*) vi.getOutNeighbors();
      uintE viOutDegree  = OUTG.get_vertex(i).getOutDegree();
      uintE viInDegree  = ING.get_vertex(i).getOutDegree();
      if (viOutDegree + viInDegree == 0) return; //return;
      // auto nghs_seq = sequence<uintE>(nghs, degree);

      uintE* outnghs_vi = (uintE*) OUTG.get_vertex(i).getOutNeighbors();
      uintE* innghs_vi = (uintE*) ING.get_vertex(i).getOutNeighbors();

      uintE vj, vk, vl;
      // i -> j -> k
      for (uintE j = 0; j < viOutDegree; j++){
        vj = outnghs_vi[j];
        uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
          vk = outnghs_vj[k];
          if (U[vk] == 0) { V->distinct[V->num_distinct] = vk; V->num_distinct++; }
          U[vk] += 1;
        }
      }
      // i <- j -> k or i <- j <- k
      for (uintE j = 0; j < viInDegree; j++){
        vj = innghs_vi[j];
        uintE* nghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < OUTG.get_vertex(vj).getOutDegree();  k++) {
          vk = nghs_vj[k];
          if (U[vk] == 0) { V->distinct[V->num_distinct] = vk; V->num_distinct++; }
          U[vk] += 1;
        }
        nghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
        for (uintE k = 0; k < ING.get_vertex(vj).getOutDegree();  k++) {
          vk = nghs_vj[k];
          if (U[vk] == 0) { V->distinct[V->num_distinct] = vk; V->num_distinct++; }
          U[vk] += 1;
        }
      }
      U[i] = 0;

      for (uintE j = 0; j < viInDegree; j++) {
        vj = innghs_vi[j];
        uintE* innghs_vj = (uintE*) ING.get_vertex(vj).getOutNeighbors();
        uintE vjInDegree = ING.get_vertex(vj).getOutDegree();
        uintE* outnghs_vj = (uintE*) OUTG.get_vertex(vj).getOutNeighbors();
        uintE vjOutDegree = OUTG.get_vertex(vj).getOutDegree();

        if (vjInDegree + vjOutDegree == 0) continue;

        //auto nghs_vj_seq = sequence<uintE>(nghs_vj, vj_degree);
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
            if (((index_leq_vk = _binary_search(outnghs_vi, custom_less_k, viOutDegree)) < viOutDegree )
                && (outnghs_vi[index_leq_vk] == vk))  tmp--;
            else if (((index_leq_vk = _binary_search(innghs_vi, custom_less_k, viInDegree)) < viInDegree )
                && (innghs_vi[index_leq_vk] == vk))  tmp--;

            if (((index_leq_vl = _binary_search(outnghs_vj, custom_less_l, vjOutDegree)) < vjOutDegree)
              && (outnghs_vj[index_leq_vl] == vl)) tmp--;
            else if  (((index_leq_vl = _binary_search(innghs_vj, custom_less_l, vjInDegree)) < vjInDegree)
              && (innghs_vj[index_leq_vl] == vl)) tmp--;
          }
        }
        // nghs_vj = nghs_vj_seq.to_array();
      }
      cycleCounts[i * eltsPerCacheLine] = tmp;
      V->reset();
  };

  auto init_V = [&](U_FastReset* V){ V->alloc(GA.n); };
  auto finish_V = [&](U_FastReset* V){ delete V; };
  parallel_for_alloc<U_FastReset>(init_V, finish_V, 0, GA.n, [&] (size_t i, U_FastReset* V) {
    run_intersection(i, V);
  });
  ulong total = pbbslib::reduce_add(cycleCounts);

  double ttc = tc.stop();
  debug(std::cout << "##### Actual counting: " << ttc << std::endl;);
  std::cout << "##### Result: " << total << "," << ttp << "," << 0 << "," << ttc << std::endl;
  ING.del();
  OUTG.del();
  return total;
}