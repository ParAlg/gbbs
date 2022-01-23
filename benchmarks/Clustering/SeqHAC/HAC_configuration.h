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

#include "gbbs/macros.h"

namespace gbbs {

struct EmptyToLogW {
  struct data {};
  template <class Graph, class WeightType = uintE>
  struct GetWeight {
    using weight_type = WeightType;
    using underlying_weight_type = gbbs::empty;
    Graph& G;

    GetWeight(Graph& G) : G(G) {}

    // Convert an underlying weight to an initial edge weight for this edge.
    weight_type get_weight(const uintE& u, const uintE& v,
                           const underlying_weight_type& wgh) const {
      auto v_u = G.get_vertex(u);
      auto v_v = G.get_vertex(v);
      return weight_type(
          1 + parlay::log2_up(1 + v_u.out_degree() +
                              v_v.out_degree()));  // [1, log(max_deg))
    }
  };
};

struct ActualWeight {
  struct data {};
  template <class Graph, class WeightType = float>
  struct GetWeight {
    using weight_type = typename Graph::weight_type;
    using underlying_weight_type = typename Graph::weight_type;
    Graph& G;

    GetWeight(Graph& G) : G(G) {}

    // Convert an underlying weight to an initial edge weight for this edge.
    weight_type get_weight(const uintE& u, const uintE& v,
                           const underlying_weight_type& wgh) const {
      return wgh;
    }
  };
};

struct DissimilarityClustering {
  struct data {};

  template <class Graph, class GetWeight = EmptyToLogW>
  struct Clustering : GetWeight::template GetWeight<Graph> {
    using base = typename GetWeight::template GetWeight<Graph>;
    using weight_type = typename base::weight_type;
    using base::base;  // import constructors

    // Used to specify whether we are doing similarity of dissimilarity
    // clustering. Similarity means taking max (heavier weights are more
    // similar)
    // and dissimilarity means taking min (smaller edges are "closer")
    static weight_type augmented_combine(const weight_type& lhs,
                                         const weight_type& rhs) {
      return std::min(lhs, rhs);  // similarity
    }

    static weight_type id() { return std::numeric_limits<weight_type>::max(); }
  };

  template <class Graph, class WeightType, class GetWeight = EmptyToLogW>
  struct WeightedClustering : GetWeight::template GetWeight<Graph, WeightType> {
    using base = typename GetWeight::template GetWeight<Graph, WeightType>;
    using weight_type = WeightType;
    using base::base;  // import constructors

    // Used to specify whether we are doing similarity of dissimilarity
    // clustering. Similarity means taking max (heavier weights are more
    // similar)
    // and dissimilarity means taking min (smaller edges are "closer")
    static weight_type augmented_combine(const weight_type& lhs,
                                         const weight_type& rhs) {
      return std::min(lhs, rhs);  // similarity
    }

    static weight_type id() { return std::numeric_limits<weight_type>::max(); }
  };
};

struct SimilarityClustering {
  struct data {};

  template <class Graph, class GetWeight = EmptyToLogW>
  struct Clustering : GetWeight::template GetWeight<Graph> {
    using base = typename GetWeight::template GetWeight<Graph>;
    using weight_type = typename base::weight_type;
    using base::base;  // import constructors

    // Used to specify whether we are doing similarity of dissimilarity
    // clustering. Similarity means taking max (heavier weights are more
    // similar)
    // and dissimilarity means taking min (smaller edges are "closer")
    static weight_type augmented_combine(const weight_type& lhs,
                                         const weight_type& rhs) {
      return std::max(lhs, rhs);  // similarity
    }

    static weight_type id() { return (weight_type)0; }
  };

  template <class Graph, class WeightType, class GetWeight = EmptyToLogW>
  struct WeightedClustering : GetWeight::template GetWeight<Graph, WeightType> {
    using base = typename GetWeight::template GetWeight<Graph, WeightType>;
    using weight_type = WeightType;
    using base::base;  // import constructors

    // Used to specify whether we are doing similarity of dissimilarity
    // clustering. Similarity means taking max (heavier weights are more
    // similar)
    // and dissimilarity means taking min (smaller edges are "closer")
    static weight_type augmented_combine(const weight_type& lhs,
                                         const weight_type& rhs) {
      return std::max(lhs, rhs);  // similarity
    }

    static weight_type id() { return (weight_type)0; }
  };
};

template <class Graph, class ClusteringType = DissimilarityClustering,
          class GetWeight = EmptyToLogW>
struct WeightedAverageLinkage
    : ClusteringType::template WeightedClustering<Graph, float, GetWeight> {
  using base =
      typename ClusteringType::template WeightedClustering<Graph, float,
                                                           GetWeight>;
  using weight_type = float;
  using base::base;

  // The linkage function.
  static weight_type linkage(const weight_type& lhs, const weight_type& rhs) {
    return (lhs + rhs) / static_cast<double>(2);
  }

  static std::string AsString(const weight_type& wgh) {
    return std::to_string(wgh);
  }
};

template <class Graph, class ClusteringType = DissimilarityClustering,
          class GetWeight = EmptyToLogW>
struct MaxLinkage : ClusteringType::template Clustering<Graph, GetWeight> {
  using base = typename ClusteringType::template Clustering<Graph, GetWeight>;
  using weight_type = typename base::weight_type;
  using base::base;

  // The linkage function.
  static weight_type linkage(const weight_type& lhs, const weight_type& rhs) {
    return std::max(lhs, rhs);
  }

  static std::string AsString(const weight_type& wgh) {
    return std::to_string(wgh);
  }
};

template <class Graph, class ClusteringType = DissimilarityClustering,
          class GetWeight = EmptyToLogW>
struct MinLinkage : ClusteringType::template Clustering<Graph, GetWeight> {
  using base = typename ClusteringType::template Clustering<Graph, GetWeight>;
  using weight_type = typename base::weight_type;
  using base::base;

  // The linkage function.
  static weight_type linkage(const weight_type& lhs, const weight_type& rhs) {
    return std::min(lhs, rhs);
  }

  static std::string AsString(const weight_type& wgh) {
    return std::to_string(wgh);
  }
};

struct NormAvgLinkWeight {
  uintE bundle_size;
  double total_weight;
  NormAvgLinkWeight() : bundle_size(0), total_weight(0) {}
  NormAvgLinkWeight(uintE single_edge_weight)
      : bundle_size(1), total_weight(single_edge_weight) {}
  NormAvgLinkWeight(uintE bundle_size, double total_weight)
      : bundle_size(bundle_size), total_weight(total_weight) {}
  double get_weight() const {
    return total_weight / static_cast<double>(bundle_size);
  }
  void print() const {
    std::cout << "{" << bundle_size << ", " << total_weight << "}" << std::endl;
  }
};
bool operator<(const NormAvgLinkWeight& l, const NormAvgLinkWeight& r) {
  return l.get_weight() < r.get_weight();
}
bool operator<=(const NormAvgLinkWeight& l, const NormAvgLinkWeight& r) {
  return l.get_weight() <= r.get_weight();
}
bool operator>(const NormAvgLinkWeight& l, const NormAvgLinkWeight& r) {
  return l.get_weight() > r.get_weight();
}
bool operator>=(const NormAvgLinkWeight& l, const NormAvgLinkWeight& r) {
  return l.get_weight() >= r.get_weight();
}
bool operator==(const NormAvgLinkWeight& l, const NormAvgLinkWeight& r) {
  return l.get_weight() == r.get_weight();
}
bool operator!=(const NormAvgLinkWeight& l, const NormAvgLinkWeight& r) {
  return l.get_weight() != r.get_weight();
}

template <class Graph, class ClusteringType = DissimilarityClustering,
          class GetWeight = EmptyToLogW>
struct NormAverageLinkage
    : ClusteringType::template WeightedClustering<Graph, NormAvgLinkWeight,
                                                  GetWeight> {
  using base = typename ClusteringType::template WeightedClustering<
      Graph, NormAvgLinkWeight, GetWeight>;
  using weight_type = NormAvgLinkWeight;
  using base::base;

  // The linkage function.
  static weight_type linkage(const weight_type& lhs, const weight_type& rhs) {
    size_t bundle_size = lhs.bundle_size + rhs.bundle_size;
    double total_weight = lhs.total_weight + rhs.total_weight;
    return weight_type(bundle_size, total_weight);
  }

  static std::string AsString(const weight_type& wgh) {
    return std::to_string(wgh.get_weight());
  }
};

struct AvgLinkWeight {
  double total_weight;    // weight going across this cut
  double current_weight;  // total_weight / (|A| * |B|)
  AvgLinkWeight() : total_weight(0), current_weight(0) {}
  AvgLinkWeight(uintE single_edge_weight)
      : total_weight(single_edge_weight), current_weight(single_edge_weight) {}
  AvgLinkWeight(double total_weight, double current_weight)
      : total_weight(total_weight), current_weight(current_weight) {}
  double get_weight() const { return current_weight; }
  void print() const {
    std::cout << "{" << total_weight << ", " << current_weight << "}"
              << std::endl;
  }
};
bool operator<(const AvgLinkWeight& l, const AvgLinkWeight& r) {
  return l.get_weight() < r.get_weight();
}
bool operator<=(const AvgLinkWeight& l, const AvgLinkWeight& r) {
  return l.get_weight() <= r.get_weight();
}
bool operator>(const AvgLinkWeight& l, const AvgLinkWeight& r) {
  return l.get_weight() > r.get_weight();
}
bool operator>=(const AvgLinkWeight& l, const AvgLinkWeight& r) {
  return l.get_weight() >= r.get_weight();
}
bool operator==(const AvgLinkWeight& l, const AvgLinkWeight& r) {
  return l.get_weight() == r.get_weight();
}
bool operator!=(const AvgLinkWeight& l, const AvgLinkWeight& r) {
  return l.get_weight() != r.get_weight();
}

template <class Graph, class ClusteringType = DissimilarityClustering,
          class GetWeight = EmptyToLogW>
struct ApproxAverageLinkage
    : ClusteringType::template WeightedClustering<Graph, AvgLinkWeight,
                                                  GetWeight> {
  using base =
      typename ClusteringType::template WeightedClustering<Graph, AvgLinkWeight,
                                                           GetWeight>;
  using weight_type = AvgLinkWeight;
  using base::base;
  using underlying_weight_type = typename Graph::weight_type;

  using value = std::pair<uintE, weight_type>;  // type of the value in the ngh
                                                // tree (neighbor_map)

  // Convert an underlying weight to an initial edge weight for this edge.
  static weight_type get_weight(const uintE& u, const uintE& v,
                                const underlying_weight_type& wgh) {
    double total_weight = wgh;  // divided by 1
    double current_weight = total_weight;
    return weight_type(total_weight, current_weight);
  }

  template <class Clusters>
  static auto GetLinkage(Clusters& clusters, const uintE& our_size) {
    return [&, our_size](const value& lhs, const value& rhs) -> value {
      uintE id = lhs.first;
      double ngh_size = clusters[id].size();
      double total_weight = lhs.second.total_weight + rhs.second.total_weight;
      double sizes_product = our_size * ngh_size;
      double current_weight = total_weight / sizes_product;
      return value(id, weight_type(total_weight, current_weight));
    };
  }

  template <class Clusters>
  static value UpdateWeight(Clusters& clusters, const value& wgh,
                            const uintE& our_size) {
    uintE ngh_id = wgh.first;
    uintE ngh_size = clusters[ngh_id].size();
    double total_weight = wgh.second.total_weight;
    double sizes_product = our_size * ngh_size;
    double current_weight = total_weight / sizes_product;
    return value(ngh_id, weight_type(total_weight, current_weight));
  }

  static std::string AsString(const weight_type& wgh) {
    return std::to_string(wgh.get_weight());
  }
};

}  // namespace gbbs
