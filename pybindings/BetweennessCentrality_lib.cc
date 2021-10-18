#include "BetweennessCentrality_lib.h"

#include "benchmarks/SSBetweenessCentrality/Brandes/SSBetweennessCentrality.h"

#include "gbbs/gbbs.h"

namespace gbbs {
  namespace compiled {

    sequence<double>
    BetweennessCentrality(symmetric_unweighted_graph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_unweighted_graph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_uint32_graph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_uint32_graph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_float_graph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_float_graph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_double_graph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_double_graph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_unweighted_cgraph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_unweighted_cgraph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_uint32_cgraph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_uint32_cgraph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_float_cgraph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_float_cgraph>(G, start);
    }

    sequence<double>
    BetweennessCentrality(symmetric_double_cgraph& G,
			  size_t start)
    {
      return gbbs::bc_bfs::SSBetweennessCentrality_BFS<symmetric_double_cgraph>(G, start);
    }

  }  // namespace compiled
}  // namespace gbbs
