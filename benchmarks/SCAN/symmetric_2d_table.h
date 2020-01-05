#pragma once

#include <vector>

// Holds an n-by-n two-dimensional table where table[i][j] = table[j][i].
class Symmetric2DTable {
 public:
  // Creates a `dimension_length`-by-`dimension_length` table.
  explicit Symmetric2DTable(size_t dimension_length);

  // Get entry at index (i, j).
  float GetEntry(size_t i, size_t j);
  // Set entry at index (i, j) to be the specified value.
  void SetEntry(size_t i, size_t j, float value);

  const size_t dimension_length_;

 private:
  std::vector<float> entries_;
};
