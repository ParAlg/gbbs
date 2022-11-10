#include "PushRelabel.h"

namespace gbbs {

template <class Graph>
double PushRelabel_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  uintE target = static_cast<uintE>(P.getOptionLongValue(
      "-target", (G.num_vertices() == 0) ? 0 : (G.num_vertices() - 1)));
  std::cout << "### Application: PushRelabel" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src << " -target = " << target
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  auto parents = PushRelabel(G, src, target);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::PushRelabel_runner, false);
