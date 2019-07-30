#pragma once

template <typename G,
          typename GW,
          typename std::enable_if<std::is_same<typename G::weight_type, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto gw(GW& unweighted_weights) {
  return unweighted_weights;
}

template <typename G,
          typename GW,
          typename std::enable_if<!std::is_same<typename G::weight_type, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto gw(GW& unweighted_weights) {
  using W = typename G::weight_type;
  // default weights
  return [](const uintE& u, const uintE& v, const W& wgh) __attribute__((always_inline)) -> uintE {
    return wgh;
  };
}



