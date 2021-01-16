#pragma once

#include <algorithm>
#include <cmath>
#include "gbbs/gbbs.h"
#include <vector>

//#include "../Makkar/dynamic_graph.h"
#include "makkar_graph.h"

namespace gbbs {

template <class Graph>
inline size_t Makkar_Dynamic_Triangle(
    Graph& G,
    const vector<gbbs::gbbs_io::Edge<int>>& updates,
    // const F& f,
    int batch_size,
    int weight,
    commandLine& P) {
  size_t n = P.getOptionLongValue("-n", 1000000);
  size_t num_batches = (updates.size() + batch_size - 1) / batch_size;

  // Just convert to sequence for convenince.
  using Edge = gbbs::gbbs_io::Edge<int>;
  pbbs::sequence<Edge> U(updates.size());
  timer t;

  if(weight == 1){ ////////////////////// insertions

  parallel_for(0, updates.size(), [&] (size_t i) {
    U[i].from = updates[i].from;
    U[i].to = updates[i].to;
    U[i].weight = updates[i].weight;
  });

  t.start();


  auto DG = gbbs::DynamicGraph(n);

  for (size_t i=0; i<num_batches; i++) {
    // process batch i
    size_t batch_start = i*batch_size;
    size_t batch_end = std::min(updates.size(), (i+1)*batch_size);
    auto batch = U.slice(batch_start, batch_end);
    timer bt; bt.start();
    DG.process_batch(batch);
    bt.stop(); bt.reportTotal("batch time");
    std::cout << "### Batch " << i << " [" << batch_start << " " << batch_end << "]" << std::endl;
    std::cout << std::endl;
  }
  DG.report_stats();


  } else if(weight == 2){ ////////////////////// deletions

    if (G.n != 0) {
      parallel_for(0, updates.size(), [&] (size_t i) {
        U[i].from = updates[i].from;
        U[i].to = updates[i].to;
        U[i].weight = 0;
      });

      t.start();
      auto DG = gbbs::DynamicGraph(G);
      t.next("graph initialized");

      for (size_t i=0; i<=num_batches; i++) {
        // process batch i
        size_t batch_start = (num_batches-i) * batch_size;
        size_t batch_end = std::min(updates.size(), batch_start + batch_size);
        if(batch_end <= batch_start) continue;
        auto batch = U.slice(batch_start, batch_end);
        timer bt; bt.start();
        DG.process_batch(batch);
        bt.stop(); bt.reportTotal("batch time");
        std::cout << "### Batch " << i << " [" << batch_start << " " << batch_end << "]" << std::endl;
        std::cout << std::endl;
      }
      DG.report_stats();
    } else {
      parallel_for(0, updates.size(), [&] (size_t i) {
        U[i].from = updates[i].from;
        U[i].to = updates[i].to;
        U[i].weight = 1;
      });

      t.start();
      auto inserts = U.slice();
      auto DG = gbbs::DynamicGraph(n);
      DG.process_batch(inserts);
      t.next("graph initialized");


      parallel_for(0, updates.size(), [&] (size_t i) {
        U[i].weight = 0;
      });

      for (size_t i=0; i<=num_batches; i++) {
        // process batch i
        size_t batch_start = (num_batches-i) * batch_size;
        size_t batch_end = std::min(updates.size(), batch_start + batch_size);
        if(batch_end <= batch_start) continue;
        auto batch = U.slice(batch_start, batch_end);
        timer bt; bt.start();
        DG.process_batch(batch);
        bt.stop(); bt.reportTotal("batch time");
        std::cout << "### Batch " << i << " [" << batch_start << " " << batch_end << "]" << std::endl;
        std::cout << std::endl;
      }
      DG.report_stats();
    }
  }




  return 0;

  t.stop();
  t.reportTotal("total processing time");
}


}
