#pragma once
#include "utils.h"
#include "map_ops.h"

// *******************************************
//   AUGMENTED MAP OPERATIONS
// *******************************************

namespace cpam {

template<class Map>
struct augmented_ops : Map {
  using Entry = typename Map::Entry;
  using node = typename Map::node;
  using Seq = typename Map::_Seq;
  using ET = typename Map::ET;
  using GC = typename Map::GC;
  using K = typename Map::K;
  using aug_t = typename Entry::aug_t;
  using ptr = typename GC::ptr;

  static inline aug_t aug_val(node* b) {
    return Seq::aug_val(b);
  }

  struct aug_sum_t {
    aug_t result;
    aug_sum_t() : result(Entry::get_empty()) {}
    void add_entry(ET e) {
      result = Entry::combine(result, Entry::from_entry(e));
    }
    void add_aug_val(aug_t av) {
      result = Entry::combine(result, av);
    }
  };

// TODO
//  template<class aug>
//  // the sum right of or at key
//  static void aug_sum_right(node* b, const K& key, aug& a) {
//    while (b) {
//      if (!Map::comp(Map::get_key(b), key)) {
//	a.add_entry(Map::get_entry(b));
//	//if (b->rc) a.add_aug_val(aug_val(b->rc));
//	if (b->rc) a.add_aug_val(b->rc->entry.second);
//	b = b->lc;
//      } else b = b->rc;
//    }
//  }

  template<class aug>
  // the sum left of or at key
  static void aug_sum_left(node* b, const K& key, aug& a) {
    while (b) {
      if (Map::is_compressed(b)) {
        bool lt = true;
        auto fn = [&] (const auto& et) {
          if (lt && Map::comp(Entry::get_key(et), key)) a.add_entry(et);
          else lt = false;
        };
        Map::iterate_seq(b, fn);
        return;
      }
      auto rb = Map::cast_to_regular(b);
      if (!Map::comp(key, Map::get_key(rb))) {
	a.add_entry(Map::get_entry(rb));
	if (rb->lc) a.add_aug_val(Map::aug_val(rb->lc));
	b = rb->rc;
      } else b = rb->lc;
    }
  }

//  template<class aug>
//  // the sum left of or at key
//  static void aug_sum_left(ptr p_b, const K& key, aug& a) {
//    if (p_b.empty()) return;
//
//    if (p_b.is_compressed()) {
//      aug_sum_left_compressed(p_b.node_ptr(), key, a);
//      return;
//    }
//
//    auto b = Map::cast_to_regular(p_b.node_ptr());
//
//
//
//
//    if (!b.empty()) {
//      auto [lc, e, rc, m] = Seq::expose(std::move(b));
//      if (!Map::comp(key, Entry::get_key(e))) {
//        a.add_entry(e);
//        if (!lc.empty()) {  // extract aug_val from ptr
//          auto l = lc.node_ptr();
//          a.add_aug_val(Map::aug_val(l));  // TODO
//        }
//        GC::decrement(m);
//        aug_sum_left(std::move(rc), key, a);
//      } else {
//        GC::decrement(m);
//        return aug_sum_left(std::move(lc), key, a);
//      }
//    }
//  }



//  template<class aug>
//  // the sum left of or at key
//  static void aug_sum_left(ptr b, const K& key, aug& a) {
//    if (!b.empty()) {
//      auto [lc, e, rc, m] = Seq::expose(std::move(b));
//      if (!Map::comp(key, Entry::get_key(e))) {
//        a.add_entry(e);
//        if (!lc.empty()) {  // extract aug_val from ptr
//          auto l = lc.node_ptr();
//          a.add_aug_val(Map::aug_val(l));  // TODO
//        }
//        GC::decrement(m);
//        aug_sum_left(std::move(rc), key, a);
//      } else {
//        GC::decrement(m);
//        return aug_sum_left(std::move(lc), key, a);
//      }
//    }
//  }


//  template<class aug>
//  static void aug_sum_range(node* b, const K& key_left, const K& key_right, aug& a) {
//    node* r = Map::range_root(b, key_left, key_right);
//    if (r) {
//      // add in left side (right of or at key_left)
//      aug_sum_right(r->lc, key_left, a);
//      a.add_entry(Map::get_entry(r));   // add in middle
//      // add in right side (left of or at key_right)
//      aug_sum_left(r->rc, key_right, a);
//    }
//  }
//
//  template<typename Func>
//  static node* aug_select(node* b, const Func& f) {
//    if (b == NULL) return NULL;
//    if (f(aug_val(b->lc))) {
//      if (f(Entry::from_entry(Map::get_entry(b)))) {
//	return aug_select(b->rc, f);
//      }
//      return b;
//    }
//    return aug_select(b->lc, f);
//  }
//
//  static node* aug_eq(node* b, aug_t query) {
//    if (b == NULL) return NULL;
//    assert(aug_val(b) == query);
//    if (Entry::from_entry(Map::get_entry(b)) == query) return b;
//    else if (aug_val(b->lc) == query) return aug_eq(b->lc, query);
//    return aug_eq(b->rc, query);
//  }
//
//  template<class Func>
//  static node* aug_filter(node* b, const Func& f, bool extra_ptr = false) {
//    if (!b) return NULL;
//    if (!f(aug_val(b))) return NULL;
//    bool copy = extra_ptr || (b->ref_cnt > 1);
//
//    auto P = utils::fork<node*>(Map::size(b) >= utils::node_limit,
//	        [&]() {return aug_filter(b->lc, f, copy);},
//		[&]() {return aug_filter(b->rc, f, copy);});
//
//    if (f(Entry::from_entry(Map::get_entry(b)))) {
//	return Map::node_join(P.first, P.second, GC::copy_if(b, copy, extra_ptr));
//    } else {
//      GC::dec_if(b, copy, extra_ptr);
//      return Map::join2(P.first, P.second);
//    }
//  }
//
//  template <class Func>
//  static node* insert_lazy(node* b, const ET& e, const Func& f) {
//    aug_t av = Entry::from_entry(e);
//    auto g = [&] (const aug_t& a) { return Entry::combine(av,a);};
//
//    auto lazy_join = [&] (node* l, node* r, node* m) {
//      m->rc = r; m->lc = l;
//      if (Map::is_balanced(m)) {
//	Map::lazy_update(m,g);
//	return m;
//      } else return Map::node_join(l,r,m);
//    };
//
//    return Map::insert_j(b, e, f, lazy_join, false);
//  }

};

}  // namespace cpam
