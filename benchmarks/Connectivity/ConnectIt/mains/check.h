#pragma once

/* ************************* Utils *************************** */

namespace gbbs {
template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  parlay::scan_inplace(flags);
  std::cout << "# n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = parlay::reduce_max(flags);
  std::cout << "# largest_cc has size: " << sz << "\n";
  return sz;
}

template <class Seq>
inline void RelabelDet(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto component_map = sequence<T>(n + 1, (T)0);
  T cur_comp = 0;
  for (size_t i = 0; i < n; i++) {
    T comp = ids[i];
    if (component_map[comp] == 0) {
      component_map[comp] = cur_comp;
      ids[i] = cur_comp;
      cur_comp++;
    } else {
      ids[i] = component_map[comp];
    }
  }
}

template <class S1, class S2>
inline void cc_check(S1& correct, S2& check) {
  RelabelDet(check);

  bool is_correct = true;
  parent max_cor = 0;
  parent max_chk = 0;
  parallel_for(0, correct.size(), [&](size_t i) {
    if ((correct[i] != check[i])) {
      is_correct = false;
      std::cout << "# at i = " << i << " cor = " << correct[i]
                << " got: " << check[i] << std::endl;
      std::cout.flush();
      assert(correct[i] == check[i]);
      abort();
    }
    if (correct[i] > max_cor) {
      gbbs::write_max(&max_cor, correct[i], std::less<parent>());
    }
    if (check[i] > max_chk) {
      gbbs::write_max(&max_chk, check[i], std::less<parent>());
    }
  });
  std::cout << "# correctness check: " << is_correct << std::endl;
  std::cout << "# max_cor = " << max_cor << " max_chk = " << max_chk
            << std::endl;
}
}  // namespace gbbs
