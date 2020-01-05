// Copyright 2019 Thomas Tseng
#pragma once

#include <vector>

// Index for an undirected graph from which clustering the graph with SCAN
// quick.
//
// Based off of SCAN index presented in "Efficient Structural Graph Clustering:
// An Index-Based Approach" by Wen et al.
//
// Vertices are 0-indexed.
class ScanIndex {
 public:
  ScanIndex(uint64_t num_vertices);

 private:
  // Gets the structural similarity between vertices `u` and `v`.
  float GetSimilarity(uint64_t u, uint64_t v);
  // Sets the structural similarity between vertices `u` and `v` to be
  // `similarity`.
  void SetSimilarity(uint64_t u, uint64_t v, float similarity);

  uint64_t num_vertices_;
  // Stores structural similarities between each pair of vertices.
  std::vector<float> similarities_;
};
