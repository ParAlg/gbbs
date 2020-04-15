#include "benchmarks/SCAN/Naive/scan_helpers.h"

namespace naive_scan {

namespace internal {

CoreBFSEdgeMapFunctions::CoreBFSEdgeMapFunctions(
    const StructuralSimilarities& similarities,
    const Clustering& current_clustering) 
  : similarities_{similarities}
  , current_clustering_{current_clustering} {}

bool 
CoreBFSEdgeMapFunctions::update(const uintE u, const uintE v, NoWeight) const {
  return current_clustering_[v].empty();
}

bool CoreBFSEdgeMapFunctions::updateAtomic(
    const uintE u, const uintE v, NoWeight) const {
  return current_clustering_[v].empty();
}

bool CoreBFSEdgeMapFunctions::cond(const uintE v) const {
  return current_clustering_[v].empty();
}

}  // namespace internal

}  // namespace naive_scan
