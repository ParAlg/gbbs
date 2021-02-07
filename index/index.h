#pragma once

#include <vector>

#ifdef USE_PAM
#include <pam/pam.h>
#else
#include <cpam/pam.h>
using namespace cpam;
#endif

auto add = [] (float a, float b) -> float {
    return a + b;};

using token = parlay::sequence<char>;

struct inv_index {
  using doc_id = int;
  using weight = float;

  using post_elt = pair<doc_id, weight>;

  struct doc_entry {
    using key_t = doc_id;
    using val_t = weight;
    static inline bool comp(key_t a, key_t b) { return a < b;}
    using aug_t = weight;
    static aug_t get_empty() {return 0;}
    static aug_t from_entry(key_t k, val_t v) {return v;}
    static aug_t combine(aug_t a, aug_t b) {
      return (b > a) ? b : a;}
    using entry_t = std::pair<key_t,val_t>;
  };

  using post_list = aug_map<doc_entry>;
  //using post_list = pam_map<doc_entry>;

  struct token_entry {
    using key_t = token;
    using val_t = post_list;
    static inline bool comp(const key_t& a, const key_t& b) {
      size_t m = std::min(a.size(), b.size());
      for (size_t i = 0; i < m; i++)
	if (a[i] != b[i]) return a[i] < b[i];
      return a.size() < b.size();
    }
  };

  using index_elt = pair<token, post_elt>;
  using index_list = pair<token, post_list>;
  using index = pam_map<token_entry>;

  index idx;

  inv_index(parlay::sequence<index_elt> const &S) {
    timer t("build index", false);
    //size_t n = S.size();
    auto reduce = [&] (parlay::slice<post_elt*,post_elt*> R) {
      // should be optimized to build sequentially in place
      //return post_list(R, add, true); };
      return post_list(R, add); };
    idx = index::multi_insert_reduce(index(), S, reduce);
    t.next("build");
  }

  inv_index() {}

  post_list get_list(const token w) {
    std::optional<post_list> p = idx.find(w);
    if (p) return *p;
    else return post_list();
  }

  static post_list And(post_list a, post_list b) {
    return post_list::map_intersect(a,b,add);}

  static post_list Or(post_list a, post_list b) {
    return post_list::map_union(a,b,add);}

  static post_list And_Not(post_list a, post_list b) {
    return post_list::map_difference(a,b);}

// Optional(todo).
//  vector<post_elt> top_k(post_list a, int k) {
//    int l = min<int>(k,a.size());
//    vector<post_elt> vec(l);
//    post_list b = a;
//    for (int i=0; i < l; i++) {
//      weight m = b.aug_val();
//      auto f = [m] (weight v) {return v < m;};
//      vec[i] = *b.aug_select(f);
//      b = post_list::remove(move(b),vec[i].first);
//    }
//    return vec;
//  }
};
