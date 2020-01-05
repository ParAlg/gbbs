// Copyright 2019 Thomas Tseng
#include "benchmarks/SCAN/scan.h"

#include <utility>

#include "utils/assert.h"

float ScanIndex::GetSimilarity(uintE u, uintE v) {
  ASSERT(0 <= u && u < num_vertices_);
  ASSERT(0 <= v && v < num_vertices_);
  if (u > v) {
    std::swap(u, v);
  }
  return similarities_[u * (u + 1) / 2 + v];
}

void ScanIndex::SetSimilarity(uintE u, uintE v, const float similarity) {
  ASSERT(0 <= u && u < num_vertices_);
  ASSERT(0 <= v && v < num_vertices_);
  if (u > v) {
    std::swap(u, v);
  }
  similarities_[u * (u + 1) / 2 + v] = similarity;
}
