#pragma once

namespace gbbs {
size_t get_edge_pri(uintE u, uintE v) {
  uintE min_u = std::min(u, v);
  uintE max_u = std::max(u, v);
  size_t key =
      (static_cast<size_t>(min_u) << 32UL) | static_cast<size_t>(max_u);
  return parlay::hash64(key);
}

struct sorted_it {
  uintE* ngh_a;
  uintE* ngh_b;
  size_t deg_a;
  size_t deg_b;
  uintE id_a;
  uintE id_b;
  size_t ct_a;
  size_t ct_b;
  size_t tot;
  sorted_it(uintE* ngh_a, uintE* ngh_b, uintE deg_a, uintE deg_b, uintE id_a,
            uintE id_b)
      : ngh_a(ngh_a),
        ngh_b(ngh_b),
        deg_a(deg_a),
        deg_b(deg_b),
        id_a(id_a),
        id_b(id_b) {
    ct_a = 0;
    ct_b = 0;
    tot = deg_a + deg_b;
  }

  /* specialize for now */
  auto lt(const uintE& a, const uintE& n_a, const uintE& b, const uintE& n_b) {
    assert(a < n_a);
    assert(b < n_b);
    auto a_p = get_edge_pri(a, n_a);
    auto b_p = get_edge_pri(b, n_b);
    if (a_p < b_p) {
      return std::make_pair(a, n_a);
    } else if (b_p < a_p) {
      return std::make_pair(b, n_b);
    } else {
      assert(a != n_a && b != n_b);
      if (std::make_pair(a, n_a) < std::make_pair(b, n_b)) {
        return std::make_pair(a, n_a);
      } else {
        return std::make_pair(b, n_b);
      }
    }
  }

  std::pair<uintE, uintE> get_next() {
    if (ct_a < deg_a) {
      if (ct_b < deg_b) {
        uintE a = ngh_a[ct_a];
        uintE b = ngh_b[ct_b];
        auto pair_a = std::make_pair(std::min(id_a, a), std::max(id_a, a));
        auto pair_b = std::make_pair(std::min(id_b, b), std::max(id_b, b));
        auto ret = lt(pair_a.first, pair_a.second, pair_b.first, pair_b.second);
        if (ret == pair_a) {
          ct_a++;
        } else {
          ct_b++;
        }
        return ret;
      } else { /* exhausted b */
        assert(ct_a < deg_a);
        auto ret = std::make_pair(id_a, ngh_a[ct_a++]);
        return std::make_pair(std::min(ret.first, ret.second),
                              std::max(ret.first, ret.second));
      }
    } else { /* exhausted a */
      assert(ct_b < deg_b);
      auto ret = std::make_pair(id_b, ngh_b[ct_b++]);
      return std::make_pair(std::min(ret.first, ret.second),
                            std::max(ret.first, ret.second));
    }
  }

  bool has_next() { return ct_a + ct_b < tot; }
};

}  // namespace gbbs
