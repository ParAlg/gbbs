#include "CliqueCounting_lib.h"

#include "benchmarks/CliqueCounting/Clique.h"

#include "gbbs/gbbs.h"

namespace gbbs {
  namespace compiled {

    size_t
    CliqueCounting(symmetric_unweighted_graph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      return gbbs::Clique<symmetric_unweighted_graph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_uint32_graph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      return gbbs::Clique<symmetric_uint32_graph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_float_graph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      return gbbs::Clique<symmetric_float_graph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_double_graph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      return gbbs::Clique<symmetric_double_graph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_unweighted_cgraph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      return gbbs::Clique<symmetric_unweighted_cgraph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_uint32_cgraph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      throw std::runtime_error("Graph type not supported");
      //return gbbs::Clique<symmetric_uint32_cgraph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_float_cgraph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      return gbbs::Clique<symmetric_float_cgraph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

    size_t
    CliqueCounting(symmetric_double_cgraph& G, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base, long recursive_level, bool approx_peel, double approx_eps)
    {
      throw std::runtime_error("Graph type not supported");
      //return gbbs::Clique<symmetric_double_cgraph>(G, k, order_type, epsilon, space_type, label, filter, use_base, recursive_level, approx_peel, approx_eps);
    }

  }  // namespace compiled
}  // namespace gbbs
