// Copyright 2019 Thomas Tseng
#include "benchmarks/SCAN/scan.h"

#include <utility>

#include "utils/assert.h"

ScanIndex::ScanIndex(uint64_t num_vertices)
  : num_vertices_{num_vertices}
  , similarities_(num_vertices_ * (num_vertices_ + 1) / 2, 0.) {}

float ScanIndex::GetSimilarity(uint64_t u, uint64_t v) {
  ASSERT(0 <= u && u < num_vertices_);
  ASSERT(0 <= v && v < num_vertices_);
  if (u > v) {
    std::swap(u, v);
  }
  return similarities_[u * (u + 1) / 2 + v];
}

void ScanIndex::SetSimilarity(uint64_t u, uint64_t v, float similarity) {
  ASSERT(0 <= u && u < num_vertices_);
  ASSERT(0 <= v && v < num_vertices_);
  if (u > v) {
    std::swap(u, v);
  }
  similarities_[u * (u + 1) / 2 + v] = similarity;
}
