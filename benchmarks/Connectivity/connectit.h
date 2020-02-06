#pragma once

#include "ligra/pbbslib/atomic_max_counter.h"
#include "ligra/pbbslib/atomic_sum_counter.h"

#include <iostream>
#include <limits.h>
#include <algorithm>
#include <unordered_map>
#include <random>

/**************************** Framework options ****************************/
enum SamplingOption {
  sample_kout, sample_bfs, sample_ldd, no_sampling
};

/* Union-Find options */
enum FindOption {
  find_compress, find_naive, find_split, find_halve, find_atomic_split, find_atomic_halve
};
enum UniteOption {
  unite, unite_early, unite_nd, unite_rem_lock, unite_rem_cas
};

/* RemCAS-specific options */
enum SpliceOption {
  split_atomic_one, halve_atomic_one, splice_simple, splice_atomic
};

/* Jayanti-specific options */
enum JayantiFindOption {
  find_twotrysplit, find_simple
};

/* LiuTarjan-specific options */
enum LiuTarjanConnectOption {
  simple_connect, parent_connect, extended_connect
};
enum LiuTarjanUpdateOption {
  simple_update, root_update
};
enum LiuTarjanShortcutOption {
  shortcut, full_shortcut
};
enum LiuTarjanAlterOption {
  alter, no_alter
};

/* Finish Algorithm Types */
enum AlgorithmType {
  union_find_type, liu_tarjan_type, shiloach_vishkin_type, label_prop_type
};

enum UpdateType {
  insertion_type,
  query_type
};

namespace connectit {

  template <FindOption find_option>
  std::string find_to_string() {
    if constexpr (find_option == find_compress) {
      return "find_compress";
    } else if constexpr (find_option == find_naive) {
      return "find_naive";
    } else if constexpr (find_option == find_split) {
      return "find_split";
    } else if constexpr (find_option == find_halve) {
      return "find_halve";
    } else if constexpr (find_option == find_atomic_split) {
      return "find_atomic_split";
    } else if constexpr (find_option == find_atomic_halve) {
      return "find_atomic_halve";
    } else {
      abort();
      return "";
    }
  }

  template <SpliceOption splice_option>
  auto splice_to_string() {
    if constexpr (splice_option == split_atomic_one) {
      return "split_atomic_one";
    } else if constexpr (splice_option == halve_atomic_one) {
      return "halve_atomic_one";
    } else {
      return "splice_atomic";
    }
  };

  template <UniteOption unite_option>
  std::string unite_to_string() {
    if constexpr (unite_option == unite) {
      return "unite";
    } else if constexpr (unite_option == unite_early) {
      return "unite_early";
    } else if constexpr (unite_option == unite_nd) {
      return "unite_nd";
    } else if constexpr (unite_option == unite_rem_lock) {
      return "unite_rem_lock";
    } else if constexpr (unite_option == unite_rem_cas) {
      return "unite_rem_cas";
    } else {
      abort();
    }
    return "";
  }

  template <SamplingOption sampling_option>
  std::string sampling_to_string() {
    if constexpr (sampling_option == sample_kout) {
      return "kout";
    } else if constexpr (sampling_option == sample_bfs) {
      return "bfs";
    } else if constexpr (sampling_option == sample_ldd) {
      return "ldd";
    } else {
      return "no_sampling";
    }
  }

  template <SamplingOption sampling_option, FindOption find_option, UniteOption unite_option>
  std::string uf_options_to_string() {
    return "uf; sample="
      + sampling_to_string<sampling_option>()
      + "; unite=" + unite_to_string<unite_option>()
      + "; find=" + find_to_string<find_option>();
  };

  template <SamplingOption sampling_option, FindOption find_option, UniteOption unite_option, SpliceOption splice_option>
  std::string uf_options_to_string() {
    return "uf; sample="
      + sampling_to_string<sampling_option>()
      + "; unite=" + unite_to_string<unite_option>()
      + "; find=" + find_to_string<find_option>() +
      + "; splice=" + splice_to_string<splice_option>();
  };

  template <JayantiFindOption find_option>
  auto jayanti_find_to_string() {
    if constexpr (find_option == find_twotrysplit) {
      return "find_twotrysplitting";
    } else {
      return "find_simple";
    }
  }

  template <SamplingOption sampling_option, JayantiFindOption find_option>
  std::string jayanti_options_to_string() {
    return "jayanti; sample="
      + sampling_to_string<sampling_option>()
      + "; find=" + jayanti_find_to_string<find_option>();
  }

  template <LiuTarjanConnectOption connect_option>
  auto connect_to_string() {
    if constexpr (connect_option == simple_connect) {
      return "connect";
    } else if constexpr (connect_option == parent_connect) {
      return "parent_connect";
    } else if constexpr (connect_option == extended_connect) {
      return "extended_connect";
    } else {
      abort();
    }
  }

  template <LiuTarjanUpdateOption update_option>
  auto update_to_string() {
    if constexpr (update_option == simple_update) {
      return "simple_update";
    } else if constexpr (update_option == root_update) {
      return "root_update";
    } else {
      abort();
    }
  }

  template <LiuTarjanShortcutOption shortcut_option>
  auto shortcut_to_string() {
    if constexpr (shortcut_option == shortcut) {
      return "shortcut";
    } else if constexpr (shortcut_option == full_shortcut) {
      return "full_shortcut";
    } else {
      abort();
    }
  }

  template <LiuTarjanAlterOption alter_option>
  auto alter_to_string() {
    if constexpr (alter_option == alter) {
      return "alter";
    } else {
      return "no_alter";
    }
  }

  template <
    SamplingOption sampling_option,
    LiuTarjanConnectOption connect_option,
    LiuTarjanUpdateOption update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption alter_option>
  std::string liu_tarjan_options_to_string() {
    return "liu_tarjan; sample=" + sampling_to_string<sampling_option>()
      + "; connect=" + connect_to_string<connect_option>()
      + "; update=" + update_to_string<update_option>()
      + "; shortcut=" + shortcut_to_string<shortcut_option>()
      + "; alter=" + alter_to_string<alter_option>();
  }


  // From gapbs/cc.c, Thanks to S. Beamer + M. Sutton for the original code this
  // snippet is based on.
  template <class Seq>
  std::pair<typename Seq::value_type, double>
  sample_frequent_element(Seq& S, uint64_t num_samples=1024) {
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

    double frac_of_graph = static_cast<double>(most_frequent->second) / num_samples;
  //  std::cout
  //    << "# Skipping largest intermediate component (ID: " << most_frequent->first
  //    << ", approx. " << (frac_of_graph * 100)
  //    << "% of the graph)" << std::endl;
    return std::make_pair(most_frequent->first, frac_of_graph);
  }

} // namespace connectit
