#pragma once

#include <algorithm>
#include <cmath>
#include "gbbs/gbbs.h"
#include <vector>

#include "../Makkar/dynamic_graph.h"

namespace gbbs {

template <class Graph>
inline size_t Makkar_Dynamic_Triangle(
    Graph& G,
    const vector<gbbs::gbbs_io::Edge<int>>& updates,
    // const F& f,
    int batch_size,
    commandLine& P) {
  size_t n = P.getOptionLongValue("-n", 1000000);
  size_t num_batches = (updates.size() + batch_size - 1) / batch_size;

  // Just convert to sequence for convenince.
  using Edge = gbbs::gbbs_io::Edge<int>;
  pbbs::sequence<Edge> U(updates.size());
  parallel_for(0, updates.size(), [&] (size_t i) {
    U[i].from = updates[i].from;
    U[i].to = updates[i].to;
    U[i].weight = updates[i].weight;
  });

  timer t;
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

  return 0;

//  t.end();
//  t.reportTotal("total processing time");
}


}