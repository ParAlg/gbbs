#pragma once
#include "gc.h"
#include "utils.h"

// *******************************************
//   SEQUENCES
// *******************************************

namespace cpam {

template<class Tree>
struct sequence_ops : Tree {
  using node = typename Tree::node;
  using regular_node = typename Tree::regular_node;
  using compressed_node = typename Tree::compressed_node;
  using ET = typename Tree::ET;
  using GC = gc<Tree>;
  using ptr = typename GC::ptr;

  using expose_tuple = std::tuple<ptr, ET&, ptr, regular_node*>;

  // Takes a compressed node and returns the node's elements as a tree
  static regular_node* to_tree_impl(ET* A, uint32_t n, uint32_t depth = 0) {
    if (n <= 0) return (regular_node*)Tree::empty();
    if (n == 1) return Tree::single(A[0]);
    if (n == B && depth > 0) {
      return (regular_node*)(Tree::make_single_compressed_node(A, n));
    }

    size_t mid = n/2;
    regular_node* m = Tree::make_regular_node(A[mid]);

    size_t l_size = mid, r_size = n-mid-1;
    auto P = utils::fork<regular_node*>(n >= utils::node_limit,
      [&]() {return to_tree_impl(A, l_size, depth + 1);},
      [&]() {return to_tree_impl(A+mid+1, r_size, depth + 1);});
    return Tree::regular_node_join(P.first, P.second, m);
  }

  static regular_node* to_tree(compressed_node* c) {
    ET tmp_arr[2*B];
    ET* arr = Tree::compressed_node_elms(c, tmp_arr);
    size_t arr_size = c->s;
    assert(arr_size > 0);
    return to_tree_impl(arr, arr_size);
  }

  static expose_tuple expose(ptr a) {
    node* p = ptr::strip_flag(a.p);
    bool extra = a.extra();
    bool other_real_refs = (Tree::ref_cnt(p) > 1);
    bool copy = extra || other_real_refs;
    bool is_regular = Tree::is_regular(p);
    if (is_regular) {
      regular_node* ret = nullptr;
      regular_node* rp = (regular_node*)p;
      if (!copy) {
        a.p = nullptr;
        ret = Tree::cast_to_regular(p);
      }
      auto ret_tup = expose_tuple(ptr(rp->lc, copy), Tree::get_entry(rp), ptr(rp->rc, copy), ret);
      return ret_tup;
    } else {
      // convert to tree
      //if (print) {
      //  std::cout << "Hit a compressed node, converting to tree. Size = " << Tree::size(p) << std::endl;
      //  size_t sz = 0;
      //  auto fn = [&] (const auto& et) {
      //    std::cout << std::get<0>(et) << " " << sz << std::endl;
      //    sz++;
      //  };
      //  Tree::iterate_seq(p, fn);
      //}
      regular_node* root = to_tree((compressed_node*)p);
      //if (print) {
      //  auto fn = [&] (const auto& et) {
      //    std::cout << et.first << std::endl;
      //  };
      //  std::cout << "printing new root" << std::endl;
      //  Tree::iterate_seq(root, fn);
      //}
      return expose_tuple(ptr(root->lc, false), Tree::get_entry(root), ptr(root->rc, false), root);
    }
  }

  static node* join(node* l, ET e, node* r, regular_node* root) {
    if (root == nullptr) {
      root = Tree::make_regular_node(e);
    }
    return Tree::node_join(l, r, root);
  }


  static node_size_t depth(node* a) {
    if (a == NULL) return 0;
    auto P = utils::fork<node_size_t>(Tree::size(a) >= utils::node_limit,
		[&]() {return depth(a->lc);},
		[&]() {return depth(a->rc);});
    return std::max(P.first, P.second) + 1;
  }

  static bool check_balance(node* a) {
    if (a == NULL) return true;
    auto P = utils::fork<bool>(Tree::size(a) >= utils::node_limit,
		[&]() {return check_balance(a->lc);},
		[&]() {return check_balance(a->rc);});
    return Tree::is_balanced(a) && P.first && P.second;
  }

  static std::optional<ET> select(ptr b, size_t rank) {
    if (b.empty()) return {};
    auto [lc, e, rc, m] = expose(std::move(b));
    size_t left_size = lc.size();
    std::optional<ET> ret;
    if (rank > left_size) {
      ret = select(std::move(rc), rank - left_size - 1);
    } else if (rank < left_size) {
      ret = select(std::move(lc), rank);
    } else {
      ret = e;
    }
    GC::decrement(m);
    return ret;
  }

// TODO: fix.
//  static node* take(node* b, size_t rank) {
//    if (!b) return NULL;
//    size_t left_size = Tree::size(b->lc);
//    if (rank < left_size) return take(b->lc, rank);
//    node* join = Tree::make_node(Tree::get_entry(b));
//    node* r = take(b->rc, rank - left_size - 1);
//    GC::increment(b->lc);
//    return Tree::node_join(b->lc, r, join);
//  }

  static node* join2_i(ptr b1, ptr b2) {
    if (b1.empty()) return b2.node_ptr();
    if (b2.empty()) return b1.node_ptr();

    if (b1.size() > b2.size()) {
      if (b1.is_compressed() || Tree::will_be_compressed(b1.size() + b2.size())) {
        auto ret = Tree::make_compressed(b1.node_ptr(), b2.node_ptr());
#ifdef DEBUG
        Tree::check_structure(ret);
#endif
        return ret;
      }
      auto[lc, e, rc, root] = expose(std::move(b1));
      node* l = lc.node_ptr();
      node* r = join2_i(std::move(rc), std::move(b2));
      auto ret = join(l, e, r, root);
#ifdef DEBUG
      Tree::check_structure(ret);
#endif
      return ret;
    } else {
      if (b2.is_compressed() || Tree::will_be_compressed(b1.size() + b2.size())) {
        return Tree::make_compressed(b1.node_ptr(), b2.node_ptr());
      }
      auto[lc, e, rc, root] = expose(std::move(b2));
      node* l = join2_i(std::move(b1), std::move(lc));
      assert(l != root);
      node* r = rc.node_ptr();
      auto ret = join(l, e, r, root);
#ifdef DEBUG
      Tree::check_structure(ret);
#endif
      return ret;
    }
  }

  static node* join2(node* b1, node* b2) {
    return join2_i(ptr(b1), ptr(b2));
  }

// TODO:
//  template<class InTree, class Func>
//  static node* map(typename InTree::node* b, const Func& f) {
//    if (!b) return NULL;
//    size_t n = InTree::size(b);
//    auto P = utils::fork<node*>(n >= utils::node_limit,
//       [&] () {return map<InTree>(b->lc, f);},
//       [&] () {return map<InTree>(b->rc, f);});
//    auto y = f(InTree::get_entry(b));
//	node* m = Tree::make_node(y);
//	return Tree::node_join(P.first, P.second, m);
//  }

  // Written without any customization for the compressed case.
  template <typename F>
  static void foreach_index(ptr a, size_t start, const F& f,
                            size_t granularity = utils::node_limit) {
    if (a.empty()) return;
    if (a.is_compressed()) {
      auto c = a.node_ptr();
      size_t i = 0;
      auto fn = [&] (ET a) {
        f(a, start + i);
        i++;
      };
      Tree::iterate_seq(c, fn);
      return;
    }
    auto[lc, e, rc, root] = expose(std::move(a));
    size_t lsize = lc.size();
    f(e, start + lsize);
    utils::fork_no_result(
        lsize >= granularity,
        [&]() { foreach_index(std::move(lc), start, f, granularity); },
        [&]() {
          foreach_index(std::move(rc), start + lsize + 1, f, granularity);
        });
    GC::decrement(root);
  }

  // F : entry x index -> bool
  template <typename F>
  static bool foreach_cond(ptr a, size_t start, const F& f,
                           size_t granularity = utils::node_limit) {
    if (a.empty()) return true;
    if (a.is_compressed()) {
      auto c = a.node_ptr();
      size_t i = 0;
      auto fn = [&] (ET a) -> bool {
        bool ret = f(a, start + i);
        i++;
        return ret;
      };
      return Tree::iterate_cond(c, fn);
    }
    auto [lc, e, rc, root] = expose(std::move(a));
    size_t lsize = lc.size();
    bool ret = f(e, start + lsize);
    if (!ret) {
      GC::decrement(root);
      return false;
    }
    auto P = utils::fork<bool>(
        false,
        [&]() { return foreach_cond(std::move(lc), start, f, granularity); },
        [&]() {
          return foreach_cond(std::move(rc), start + lsize + 1, f, granularity);
    });
    GC::decrement(root);
    return P.first && P.second;
  }

  // F : entry x index -> bool
  template <typename F, typename C>
  static bool foreach_cond_par(ptr a, size_t start, const F& f, const C& cond,
                               size_t granularity = utils::node_limit) {
    if (a.empty()) return true;
    if (a.is_compressed()) {
      auto c = a.node_ptr();
      size_t i = 0;
      auto fn = [&] (const ET& a) -> bool {
        bool ret = f(a, start + i);
        i++;
        return ret;
      };
      return Tree::iterate_cond(c, fn);
    }
    auto [lc, e, rc, root] = expose(std::move(a));
    size_t lsize = lc.size();
    bool ret = cond();
    if (ret) {
      ret = f(e, start + lsize);
    }
    if (!ret) {
      GC::decrement(root);
      return false;
    }
    auto P = utils::fork<bool>(
        true,
        [&]() { return foreach_cond_par(std::move(lc), start, f, cond, granularity); },
        [&]() {
          return foreach_cond_par(std::move(rc), start + lsize + 1, f, cond, granularity);
    });
    GC::decrement(root);
    return P.first && P.second;
  }

  // similar to above but sequential using in-order traversal
  // usefull if using 20th century constructs such as iterators
  template<typename F>
  static void foreach_seq(ptr a, const F& f) {
    if (a.empty()) return;
    auto[lc, e, rc, root] = expose(std::move(a));
    foreach_seq(std::move(lc), f);
    f(e);
    foreach_seq(std::move(rc), f);
    GC::decrement(root);
  }

  template <class Func>
  static node* filter_bc(ptr b1, const Func& f) {
    assert(b1.size() > 0);
    ET stack[utils::kBaseCaseSize + 1];

    auto b1_node = b1.node_ptr();
    size_t offset = 0;
    auto copy_f = [&] (const ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Tree::iterate_seq(b1_node, copy_f);
    assert(offset <= utils::kBaseCaseSize);

    Tree::decrement_recursive(b1_node);

    size_t k = 0;
    for (size_t i=0; i<offset; i++) {
      if (f(stack[i])) {
        if (i > k) {
          stack[k++] = std::move(stack[i]);
        }
      }
    }

    if (k < utils::compression_block_size) {
      return to_tree_impl((ET*)stack, k);
    } else {
      // need to refactor
      return Tree::make_compressed(stack, k);
    }
  }

  template<class Func>
  static node* filter(ptr b, const Func& f, size_t granularity=utils::node_limit) {
    if (b.empty()) return NULL;
    size_t n = b.size();

    if (n <= utils::kBaseCaseSize) {
      return filter_bc(std::move(b), f);
    }

    auto[lc, e, rc, root] = expose(std::move(b));

    auto [l, r] = utils::fork<node*>(n >= granularity,
      [&]() {return filter(std::move(lc), f, granularity);},
      [&]() {return filter(std::move(rc), f, granularity);});

    if (f(e)) {
      if (!root) root = Tree::single(e);
      return join(l, e, r, root);
    } else {
      GC::decrement(root);
      return join2(l, r);
    }
  }

// TODO
//  template<class Func>
//  static bool if_exist(node* b, const Func& f, bool* indicator) {
//    if (!b) return false;
//	if (*indicator) return true;
//	if (f(Tree::get_entry(b))) {
//		utils::atomic_compare_and_swap(indicator, false, true);
//		//*indicator = true;
//		return true;
//	}
//
//    auto P = utils::fork<bool>(Tree::size(b) >= utils::node_limit,
//      [&]() {return if_exist(b->lc, f, indicator);},
//      [&]() {return if_exist(b->rc, f, indicator);});
//	return P.first || P.second;
//  }

  // Assumes the input is sorted and there are no duplicate keys
  static node* from_array(ET* A, size_t n) {
	  //cout << "from_array: " << n << endl;
    if (n <= 0) return Tree::empty();
    if (n == 1) return Tree::single(A[0]);
    if (n >= B && n <= 2*B) {
      return Tree::make_compressed(A, n);
    }
    size_t mid = n/2;

    regular_node* m = Tree::make_regular_node(A[mid]);

    auto P = utils::fork<node*>(n >= utils::node_limit,
      [&]() {return from_array(A, mid);},
      [&]() {return from_array(A+mid+1, n-mid-1);});

    return Tree::node_join(P.first, P.second, m);
  }

// TODO
//  template<class Seq1, class Func>
//  static node* map_filter(typename Seq1::node* b, const Func& f,
//			  size_t granularity=utils::node_limit) {
//    if (!b) return NULL;
//
//    auto P = utils::fork<node*>(Seq1::size(b) >= granularity,
//      [&]() {return map_filter<Seq1>(b->lc, f, granularity);},
//      [&]() {return map_filter<Seq1>(b->rc, f, granularity);});
//
//    std::optional<ET> me = f(Seq1::get_entry(b));
//    if (me.has_value()) {
//      node* r = Tree::make_node(*me);
//      return Tree::node_join(P.first, P.second, r);
//    } else return join2(P.first, P.second);
//  }

// TODO
//  template<class R, class F>
//  static typename R::T map_reduce(node* b, F f, R r,
//				  size_t grain=utils::node_limit) {
//    using T = typename R::T;
//    if (!b) return r.identity();
//
//    auto P = utils::fork<T>(Tree::size(b) >= grain,
//      [&]() {return map_reduce<R>(b->lc, f, r, grain);},
//      [&]() {return map_reduce<R>(b->rc, f, r, grain);});
//
//    T v = f(Tree::get_entry(b));
//    return R::add(P.first, r.add(v, P.second));
//  }

// TODO
//  template<class F, class T>
//  static void semi_map_reduce_seq(node* b, T& v, F& f) {
//    if (!b) return;
//    semi_map_reduce_seq(b->lc, v, f);
//    f(v, Tree::get_entry(b));
//    semi_map_reduce_seq(b->rc, v, f);
//  }

// TODO
//  // Class R must be a structure representing a monoid with:
//  //   T = a type
//  //   T identity();
//  //   T add(T, T);  an associative operation on T
//  // Function f must be of type (T&, E) -> void
//  //   it side affects its first argument by adding to it
//  template<class R, class F>
//  static typename R::T semi_map_reduce(node* b, F& f, R reduce, size_t grain) {
//    using T = typename R::T;
//    if (Tree::size(b) >= grain) {
//      T r, l;
//      auto left = [&] () -> void {
//	r = semi_map_reduce<R,F>(b->rc, f, reduce, grain);};
//      auto right = [&] () -> void {
//	l = semi_map_reduce<R,F>(b->lc, f, reduce, grain);};
//      parlay::par_do(left,right);
//      f(l, Tree::get_entry(b));
//      return reduce.add(l,r);
//    } else {
//      // when going sequential, create identity and then update it
//      T id = reduce.identity();
//      semi_map_reduce_seq(b, id, f);
//      return id;
//    }
//  }

};

}  // namespace cpam
