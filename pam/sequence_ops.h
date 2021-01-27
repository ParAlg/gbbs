#pragma once
#include "gc.h"
#include "utils.h"

// *******************************************
//   SEQUENCES
// *******************************************

template<class Tree>
struct sequence_ops : Tree {
  using node = typename Tree::node;
  using ET = typename Tree::ET;
  using GC = gc<Tree>;

  static node* join(node* l, ET e, node* r) {
    node *x = Tree::make_node(e);
    return Tree::node_join(l, r, x);
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

  static node* select(node* b, size_t rank) {
    size_t lrank = rank;
    while (b) {
      size_t left_size = Tree::size(b->lc);
      if (lrank > left_size) {
	lrank -= left_size + 1;
	b = b->rc;
      }
      else if (lrank < left_size) b = b->lc;
      else return b;
    }
    return NULL;
  }

  static node* take(node* b, size_t rank) {
    if (!b) return NULL;
    size_t left_size = Tree::size(b->lc);
    if (rank < left_size) return take(b->lc, rank);
    node* join = Tree::make_node(Tree::get_entry(b));
    node* r = take(b->rc, rank - left_size - 1);
    GC::increment(b->lc);
    return Tree::node_join(b->lc, r, join);
  }

  static node* join2_i(node* b1, node* b2, bool extra_b1, bool extra_b2) {
    if (!b1) return GC::inc_if(b2, extra_b2);
    if (!b2) return GC::inc_if(b1, extra_b1);

    if (Tree::size(b1) > Tree::size(b2)) {
      bool copy_b1 = extra_b1 || b1->ref_cnt > 1;
      node* l = GC::inc_if(b1->lc, copy_b1);
      node* r = join2_i(b1->rc, b2, copy_b1, extra_b2);
      return Tree::node_join(l, r, GC::copy_if(b1, copy_b1, extra_b1));
    }
    else {
      bool copy_b2 = extra_b2 || b2->ref_cnt > 1;
      node* l = join2_i(b1, b2->lc, extra_b1, copy_b2);
      node* r = GC::inc_if(b2->rc, copy_b2);
      return Tree::node_join(l, r, GC::copy_if(b2, copy_b2, extra_b2));
    }
  }

  static node* join2(node* b1, node* b2) {
    return join2_i(b1, b2, false, false);
  }

  template<class InTree, class Func>
  static node* map(typename InTree::node* b, const Func& f) {
    if (!b) return NULL;
    size_t n = InTree::size(b);
    auto P = utils::fork<node*>(n >= utils::node_limit,
       [&] () {return map<InTree>(b->lc, f);},
       [&] () {return map<InTree>(b->rc, f);});
    auto y = f(InTree::get_entry(b));
	node* m = Tree::make_node(y);
    //node* r = join(P.first, y , P.second);
	return Tree::node_join(P.first, P.second, m);
    //return r;
  }

  template<typename F>
  static void foreach_index(node* a, size_t start, const F& f,
			    size_t granularity=utils::node_limit,
			    bool extra_ptr = false) {

    if (!a) return;
    bool copy = extra_ptr || (a->ref_cnt > 1);
    size_t lsize = Tree::size(a->lc);
    f(Tree::get_entry(a), start+lsize);
    utils::fork_no_result(lsize >= granularity,
      [&] () {foreach_index(a->lc, start, f, granularity, copy);},
      [&] () {foreach_index(a->rc, start + lsize + 1,f, granularity, copy);});
    if (!extra_ptr) GC::decrement(a);
  }

  // similar to above but sequential using in-order traversal
  // usefull if using 20th century constructs such as iterators
  template<typename F>
  static void foreach_seq(node* a, const F& f, bool extra_ptr = false) {
    if (!a) return;
    bool copy = extra_ptr || (a->ref_cnt > 1);
    foreach_seq(a->lc, f, copy);
    f(Tree::get_entry(a));
    foreach_seq(a->rc, f, copy);
    if (!extra_ptr) GC::decrement(a);
  }

  template<class Func>
  static node* filter(node* b, const Func& f, size_t granularity=utils::node_limit, bool extra_ptr = false) {
    if (!b) return NULL;
    bool copy = extra_ptr || (b->ref_cnt > 1);

    auto P = utils::fork<node*>(Tree::size(b) >= granularity,
      [&]() {return filter(b->lc, f, granularity, copy);},
      [&]() {return filter(b->rc, f, granularity, copy);});

    if (f(Tree::get_entry(b))) {
      return Tree::node_join(P.first, P.second, GC::copy_if(b, copy, extra_ptr));
    } else {
      GC::dec_if(b, copy, extra_ptr);
      return join2(P.first, P.second);
    }
  }


  template<class Func>
  static bool if_exist(node* b, const Func& f, bool* indicator) {
    if (!b) return false;
	if (*indicator) return true;
	if (f(Tree::get_entry(b))) {
		utils::atomic_compare_and_swap(indicator, false, true);
		//*indicator = true;
		return true;
	}

    auto P = utils::fork<bool>(Tree::size(b) >= utils::node_limit,
      [&]() {return if_exist(b->lc, f, indicator);},
      [&]() {return if_exist(b->rc, f, indicator);});
	return P.first || P.second;
  }

  /*
  template<class Func>
  static bool if_exist(node* b, const Func& f) {
    if (!b) return false;
	if (f(Tree::get_entry(b))) {
		return true;
	}

    auto P = utils::fork<bool>(Tree::size(b) >= utils::node_limit,
      [&]() {return if_exist(b->lc, f);},
      [&]() {return if_exist(b->rc, f);});
	return P.first || P.second;
  }*/

  // Assumes the input is sorted and there are no duplicate keys
  static node* from_array(ET* A, size_t n) {
	  //cout << "from_array: " << n << endl;
    if (n <= 0) return Tree::empty();
    if (n == 1) return Tree::single(A[0]);
    size_t mid = n/2;

    node* m = Tree::make_node(A[mid]);

    auto P = utils::fork<node*>(n >= utils::node_limit,
      [&]() {return from_array(A, mid);},
      [&]() {return from_array(A+mid+1, n-mid-1);});

    return Tree::node_join(P.first, P.second, m);
  }

  template<class Seq1, class Func>
  static node* map_filter(typename Seq1::node* b, const Func& f,
			  size_t granularity=utils::node_limit) {
    if (!b) return NULL;

    auto P = utils::fork<node*>(Seq1::size(b) >= granularity,
      [&]() {return map_filter<Seq1>(b->lc, f, granularity);},
      [&]() {return map_filter<Seq1>(b->rc, f, granularity);});

    std::optional<ET> me = f(Seq1::get_entry(b));
    if (me.has_value()) {
      node* r = Tree::make_node(*me);
      return Tree::node_join(P.first, P.second, r);
    } else return join2(P.first, P.second);
  }

  template<class R, class F>
  static typename R::T map_reduce(node* b, F f, R r,
				  size_t grain=utils::node_limit) {
    using T = typename R::T;
    if (!b) return r.identity();

    auto P = utils::fork<T>(Tree::size(b) >= grain,
      [&]() {return map_reduce<R>(b->lc, f, r, grain);},
      [&]() {return map_reduce<R>(b->rc, f, r, grain);});

    T v = f(Tree::get_entry(b));
    return R::add(P.first, r.add(v, P.second));
  }

  // template<class T, class Map, class Reduce>
  // static T map_reduce_index(node* a, Map m, Reduce r, T identity,
  // 			    size_t start, size_t grain=utils::node_limit) {
  //   if (!a) return identity;
  //   size_t lsize = Tree::size(a->lc);
  //   auto P = utils::fork<T>(Tree::size(a) >= grain,
  //     [&] () {return map_reduce_index(a->lc, m, r, identity, start, grain);},
  //     [&] () {return map_reduce_index(a->rc, m, r, identity, start + lsize + 1, grain);});
  //   T v = m(Tree::get_entry(a), start+lsize);
  //   return r(P.first, r(v, P.second));
  // }

  template<class F, class T>
  static void semi_map_reduce_seq(node* b, T& v, F& f) {
    if (!b) return;
    semi_map_reduce_seq(b->lc, v, f);
    f(v, Tree::get_entry(b));
    semi_map_reduce_seq(b->rc, v, f);
  }

  // Class R must be a structure representing a monoid with:
  //   T = a type
  //   T identity();
  //   T add(T, T);  an associative operation on T
  // Function f must be of type (T&, E) -> void
  //   it side affects its first argument by adding to it
  template<class R, class F>
  static typename R::T semi_map_reduce(node* b, F& f, R reduce, size_t grain) {
    using T = typename R::T;
    if (Tree::size(b) >= grain) {
      T r, l;
      auto left = [&] () -> void {
	r = semi_map_reduce<R,F>(b->rc, f, reduce, grain);};
      auto right = [&] () -> void {
	l = semi_map_reduce<R,F>(b->lc, f, reduce, grain);};
      parlay::par_do(left,right);
      f(l, Tree::get_entry(b));
      return reduce.add(l,r);
    } else {
      // when going sequential, create identity and then update it
      T id = reduce.identity();
      semi_map_reduce_seq(b, id, f);
      return id;
    }
  }

};
