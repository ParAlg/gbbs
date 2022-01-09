#pragma once

#include "gbbs/helpers/atomic_max_counter.h"
#include "gbbs/helpers/atomic_sum_counter.h"

#include <limits.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

namespace gbbs {

/**************************** Framework options ****************************/
enum SamplingOption { sample_kout, sample_bfs, sample_ldd, no_sampling };

/* Union-Find options */
enum FindOption {
  find_compress,
  find_naive,
  find_split,
  find_halve,
  find_atomic_split,
  find_atomic_halve
};
enum UniteOption {
  unite,
  unite_early,
  unite_nd,
  unite_rem_lock,
  unite_rem_cas
};

/* RemCAS-specific options */
enum SpliceOption {
  split_atomic_one,
  halve_atomic_one,
  splice_simple,
  splice_atomic
};

/* Jayanti-specific options */
enum JayantiFindOption { find_twotrysplit, find_simple };

/* LiuTarjan-specific options */
enum LiuTarjanConnectOption {
  simple_connect,
  parent_connect,
  extended_connect
};
enum LiuTarjanUpdateOption { simple_update, root_update };
enum LiuTarjanShortcutOption { shortcut, full_shortcut };
enum LiuTarjanAlterOption { alter, no_alter };

/* Finish Algorithm Types */
enum AlgorithmType {
  union_find_type,
  liu_tarjan_type,
  shiloach_vishkin_type,
  label_prop_type
};

enum UpdateType { insertion_type, query_type };

namespace connectit {

// Converts enum to string.
std::string find_to_string(FindOption);
std::string splice_to_string(SpliceOption);
std::string unite_to_string(UniteOption);
std::string sampling_to_string(SamplingOption);
std::string jayanti_find_to_string(JayantiFindOption);
std::string connect_to_string(LiuTarjanConnectOption);
std::string update_to_string(LiuTarjanUpdateOption);
std::string shortcut_to_string(LiuTarjanShortcutOption);
std::string alter_to_string(LiuTarjanAlterOption);

std::string uf_options_to_string(SamplingOption, FindOption, UniteOption);
std::string uf_options_to_string(SamplingOption, FindOption, UniteOption,
                                 SpliceOption);
std::string jayanti_options_to_string(SamplingOption sampling_option,
                                      JayantiFindOption find_option);
std::string liu_tarjan_options_to_string(SamplingOption, LiuTarjanConnectOption,
                                         LiuTarjanUpdateOption,
                                         LiuTarjanShortcutOption,
                                         LiuTarjanAlterOption);

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

}  // namespace connectit
}  // namespace gbbs
