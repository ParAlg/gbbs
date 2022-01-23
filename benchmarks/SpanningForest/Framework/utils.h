#pragma once

#include <limits.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_map>

namespace gbbs {
// From gapbs/cc.c, Thanks to S. Beamer + M. Sutton for the original code this
// snippet is based on.
template <class Seq>
std::pair<typename Seq::value_type, double> sample_frequent_element(
    Seq& S, uint64_t num_samples = 1024) {
  using VT = typename Seq::value_type;
  using T = typename std::remove_volatile<VT>::type;
  std::unordered_map<T, int> sample_counts(32);
  using kvp_type = typename std::unordered_map<T, int>::value_type;
  // Sample elements from 'S'
  std::mt19937 gen;
  std::uniform_int_distribution<T> distribution(0, S.size() - 1);
  for (T i = 0; i < num_samples; i++) {
    T n = distribution(gen);
    sample_counts[(T)(S[n])]++;
  }
  // Find most frequent element in samples (estimate of most frequent overall)
  auto most_frequent = std::max_element(
      sample_counts.begin(), sample_counts.end(),
      [](const kvp_type& a, const kvp_type& b) { return a.second < b.second; });

  double frac_of_graph =
      static_cast<double>(most_frequent->second) / num_samples;
  //  std::cout
  //    << "# Skipping largest intermediate component (ID: " <<
  //    most_frequent->first
  //    << ", approx. " << (frac_of_graph * 100)
  //    << "% of the graph)" << std::endl;
  return std::make_pair(most_frequent->first, frac_of_graph);
}
}  // namespace gbbs
