#include "StronglyConnectedComponents_lib.h"

#include "benchmarks/StronglyConnectedComponents/RandomGreedyBGSS16/StronglyConnectedComponents.h"

#include "gbbs/gbbs.h"

namespace gbbs {
  namespace compiled {

    sequence<size_t>
    StronglyConnectedComponents(symmetric_unweighted_graph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_unweighted_graph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_uint32_graph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_uint32_graph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_float_graph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_float_graph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_double_graph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_double_graph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_unweighted_cgraph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_unweighted_cgraph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_uint32_cgraph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_uint32_cgraph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_float_cgraph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_float_cgraph>(G, beta);
    }

    sequence<size_t>
    StronglyConnectedComponents(symmetric_double_cgraph& G,
				double beta)
    {
      return gbbs::StronglyConnectedComponents<symmetric_double_cgraph>(G, beta);
    }

  }  // namespace compiled
}  // namespace gbbs
