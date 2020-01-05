#include "benchmarks/SCAN/symmetric_2d_table.h"

#include <utility>

#include "utils/assert.h"

Symmetric2DTable::Symmetric2DTable(const size_t dimension_length)
  : dimension_length_{dimension_length}
  , entries_(dimension_length * (dimension_length + 1) / 2, 0.) {}

float Symmetric2DTable::GetEntry(size_t i, size_t j) {
  ASSERT(0 <= i && i < dimension_length_);
  ASSERT(0 <= j && j < dimension_length_);
  if (i > j) {
    std::swap(i, j);
  }
  return entries_[i * (i + 1) / 2 + j];
}

void Symmetric2DTable::SetEntry(size_t i, size_t j, const float similarity) {
  ASSERT(0 <= i && i < dimension_length_);
  ASSERT(0 <= j && j < dimension_length_);
  if (i > j) {
    std::swap(i, j);
  }
  entries_[i * (i + 1) / 2 + j] = similarity;
}
