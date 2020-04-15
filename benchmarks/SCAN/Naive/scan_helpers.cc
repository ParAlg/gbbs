#include "benchmarks/SCAN/Naive/scan_helpers.h"

#include <limits>

namespace naive_scan {

namespace internal {

CoreBFSEdgeMapFunctions::CoreBFSEdgeMapFunctions(
    const StructuralSimilarities& similarities,
    const Clustering& current_clustering,
    const float epsilon)
  : similarities_{similarities}
  , current_clustering_{current_clustering}
  , epsilon_{epsilon} {}

bool
CoreBFSEdgeMapFunctions::update(const uintE u, const uintE v, NoWeight) const {
  constexpr float
    kDefaultSimilarity{std::numeric_limits<float>::signaling_NaN()};
  return current_clustering_[v].empty()
    && similarities_.find({u, v}, kDefaultSimilarity) >= epsilon_;
}

bool CoreBFSEdgeMapFunctions::updateAtomic(
    const uintE u, const uintE v, NoWeight) const {
  return update(u, v, NoWeight{});
}

bool CoreBFSEdgeMapFunctions::cond(const uintE v) const {
  return current_clustering_[v].empty();
}

}  // namespace internal

}  // namespace naive_scan
