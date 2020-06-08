#pragma once

#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "gmock/gmock.h"

// Googletest-like matcher for EdgeSimilarity.
inline auto EdgeSimilarityEq(
    const uintE expected_source,
    const uintE expected_neighbor,
    const float expected_similarity) {
  return ::testing::AllOf(
   ::testing::Field(
     &scan::EdgeSimilarity::source, ::testing::Eq(expected_source)),
   ::testing::Field(
     &scan::EdgeSimilarity::neighbor, ::testing::Eq(expected_neighbor)),
   ::testing::Field(
     &scan::EdgeSimilarity::similarity,
     ::testing::FloatEq(expected_similarity)));
}
