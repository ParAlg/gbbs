#pragma once
#include "utils.h"
#include "parlay/primitives.h"
#include "parlay/internal/binary_search.h"

// *******************************************
//   MAPS and SETS
// *******************************************

template<class Seq, class EntryT>
struct map_ops : Seq {
  using Entry = EntryT;
  using node = typename Seq::node;
  using ET = typename Seq::ET;
  using GC = typename Seq::GC;
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;

  static bool comp(K a, K b) { return Entry::comp(a,b);}
  static K get_key(node *s) { return Entry::get_key(*Seq::get_entry_p(s));}
  static V get_val(node *s) { return Entry::get_val(*Seq::get_entry_p(s));}
  //static K get_key(node *s) { return Entry::get_key(Seq::get_entry(s));}
  //static V get_val(node *s) { return Entry::get_val(Seq::get_entry(s));}

  static node* find(node* b, const K& key) {
    while (b) {
      if (Entry::comp(key, get_key(b))) b = b->lc;
      else if (Entry::comp(get_key(b), key)) b = b->rc;
      else return b;
    }
    return NULL;
  }

  template <class BinaryOp>
  static inline void update_value(node* a, const BinaryOp& op) {
    ET re = Seq::get_entry(a);
    Entry::set_val(re, op(re));
    Seq::set_entry(a, re);
  }


  /*template<class Func>
  static node* update(node* b, const K& key, Func& f) {
    node* nd = find(b, key);
	if (!nd) return nd;
	ET e = Seq::get_entry(nd);
	V new_value = f(e);
	ET new_entry = make_pair(key, new_value);
	auto replace = [](V a, V b) {return b;};
	return insert(b, new_entry, replace);
  }*/

  static node* previous(node* b, const K& key) {
    node* r = NULL;
    while (b) {
      if (Entry::comp(get_key(b), key)) {r = b; b = b->rc;}
      else b = b->lc;
    }
    return r;
  }

  static node* next(node* b, const K& key) {
    node* r = NULL;
    while (b) {
      if (Entry::comp(key, get_key(b)) ) {r = b; b = b->lc;}
      else b = b->rc;
    }
    return r;
  }

  static size_t rank(node* b, const K& key) {
    size_t ret = 0;
    while (b) {
      if (Entry::comp(get_key(b), key) ) {
	ret += 1 + Seq::size(b->lc);
	b = b->rc;
      } else b = b->lc;
    }
    return ret;
  }

  struct split_info {
    split_info(node* first, node* second, bool removed)
      : first(first), second(second), removed(removed) {};
    split_info(node* first, ET entry, node* second)
      : first(first), second(second), entry(entry), removed(true) {};
    node* first;  node* second; ET entry; bool removed;
  };

  static split_info split_copy(node* bst, const K& e) {
    if (!bst) return split_info(NULL, NULL, false);

    else if (Entry::comp(get_key(bst), e)) {
      // following should be a copy ????
      node* join = Seq::make_node(Seq::get_entry(bst));
      split_info bstpair = split_copy(bst->rc, e);
      GC::increment(bst->lc);
      bstpair.first = Seq::node_join(bst->lc, bstpair.first, join);
      return bstpair;
    }
    else if (Entry::comp(e, get_key(bst))) {
      node* join = Seq::make_node(Seq::get_entry(bst));
      split_info bstpair = split_copy(bst->lc, e);
      GC::increment(bst->rc);
      bstpair.second = Seq::node_join(bstpair.second, bst->rc, join);
      return bstpair;
    }
    else {
      GC::increment(bst->lc); GC::increment(bst->rc);
      split_info ret(bst->lc, bst->rc, true);
      ret.entry = Seq::get_entry(bst);
      return ret;
    }
  }

  // A version that will reuse a node if ref count is 1
  // Will decrement ref count of root if it is copied
  static split_info split_inplace(node* bst, const K& e) {
    if (!bst) return split_info(NULL, NULL, false);
    else if (bst->ref_cnt > 1) {
		//std::cout << "rc>1" << std::endl;
      split_info ret = split_copy(bst, e);
      GC::decrement_recursive(bst);
      return ret;
    }
    else if (Entry::comp(get_key(bst), e)) {
      split_info bstpair = split_inplace(bst->rc, e);
      bstpair.first = Seq::node_join(bst->lc, bstpair.first, bst);
      return bstpair;
    }
    else if (Entry::comp(e, get_key(bst))) {
      split_info bstpair = split_inplace(bst->lc, e);
      bstpair.second = Seq::node_join(bstpair.second, bst->rc, bst);
      return bstpair;
    }
    else {
      split_info ret(bst->lc, bst->rc, true);
      ret.entry = Seq::get_entry(bst);
      GC::decrement(bst);
      return ret;
    }
  }

  static inline split_info split(node* bst, const K& e) {
    return split_inplace(bst, e);
  }

  template <class BinaryOp>
  static inline void combine_values(node* a, ET e, bool reverse, const BinaryOp& op) {
    ET re = Seq::get_entry(a);
	//std::cout << Entry::get_val(re) << " " << Entry::get_val(e) << std::endl;
    if (reverse) Entry::set_val(re, op(Entry::get_val(e), Entry::get_val(re)));
    else Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(e)));
    Seq::set_entry(a, re);
  }

  template <class VE, class BinaryOp>
  static inline void update_valuev(node* a, VE v0, const BinaryOp& op) {
    ET re = Seq::get_entry(a);
    Entry::set_val(re, op(Entry::get_val(re), v0));
    Seq::set_entry(a, re);
  }

  //valid only if key type is printable
  static void output_r(node* t) {
	  if (!t) return;
	  output_r(t->lc);
	  std::cout << get_key(t) << " ";
	  output_r(t->rc);
  }

  //valid only if key type is printable
  static void output(node* t) {
	  output_r(t);
	  std::cout << std::endl;
  }

  // Works in-place when possible.
  // extra_b2 means there is an extra pointer to b2 not included
  // in the reference count.  It is used as an optimization to reduce
  // the number of ref_cnt updates.
  template <class BinaryOp>
  static node* uniont(node* b1, node* b2, const BinaryOp op,
		      bool extra_b2 = false) {
    if (!b1) return GC::inc_if(b2, extra_b2);
    if (!b2) return b1;

    size_t n1 = Seq::size(b1);   size_t n2 = Seq::size(b2);

    bool copy = extra_b2 || (b2->ref_cnt > 1);
    node* r = copy ? Seq::make_node(Seq::get_entry(b2)) : b2;
    split_info bsts = split(b1, get_key(b2));

    auto P = utils::fork<node*>(utils::do_parallel(n1, n2),
      [&] () {return uniont(bsts.first, b2->lc, op, copy);},
      [&] () {return uniont(bsts.second, b2->rc, op, copy);});

    if (copy && !extra_b2) GC::decrement_recursive(b2);
    if (bsts.removed) combine_values(r, bsts.entry, true, op);
    return Seq::node_join(P.first, P.second, r);
  }

  // extra_b2 is for b2
  template <class Seq1, class Seq2, class BinaryOp>
  static node* intersect(typename Seq1::node* b1, typename Seq2::node* b2,
			 const BinaryOp& op,
			 bool extra_b2 = false) {
	//std::cout << "size: " << n1 << " " << n2 << std::endl;
    if (!b1) {if (!extra_b2) Seq2::GC::decrement_recursive(b2); return NULL;}
    if (!b2) {Seq1::GC::decrement_recursive(b1); return NULL;}
    size_t n1 = Seq1::size(b1);   size_t n2 = Seq2::size(b2);
    bool copy = extra_b2 || (b2->ref_cnt > 1);
    typename Seq1::split_info bsts = Seq1::split(b1, Seq2::get_key(b2));

    auto P = utils::fork<node*>(utils::do_parallel(n1, n2),
	[&]() {return intersect<Seq1,Seq2>(bsts.first, b2->lc, op, copy);},
	[&]() {return intersect<Seq1,Seq2>(bsts.second, b2->rc, op, copy);}
    );

    if (bsts.removed) {
      ET e(Seq2::get_key(b2),
	   op(Seq1::Entry::get_val(bsts.entry),
	      Seq2::Entry::get_val(Seq2::get_entry(b2))));
      node* r = Seq::make_node(e);
      //if (copy && !extra_b2) Seq2::GC::decrement_recursive(b2);
      Seq2::GC::dec_if(b2, copy, extra_b2);
      return Seq::node_join(P.first, P.second, r);
    } else {
      //if (copy && !extra_b2) Seq2::GC::decrement_recursive(b2);
      Seq2::GC::dec_if(b2, copy, extra_b2);
      return Seq::join2(P.first, P.second);
    }
  }

  static node* difference(node* b1, node* b2, bool extra_b1 = false) {
    if (!b1) {GC::decrement_recursive(b2); return NULL;}
    if (!b2) return GC::inc_if(b1, extra_b1);
    size_t n1 = Seq::size(b1);   size_t n2 = Seq::size(b2);
    bool copy = extra_b1 || (b1->ref_cnt > 1);
    split_info bsts = split(b2, get_key(b1));

    auto P = utils::fork<node*>(utils::do_parallel(n1, n2),
      [&]() {return difference(b1->lc, bsts.first, copy);},
      [&]() {return difference(b1->rc, bsts.second, copy);});

    if (bsts.removed) {
      GC::dec_if(b1, copy, extra_b1);
      return Seq::join2(P.first, P.second);
    } else {
      return Seq::node_join(P.first, P.second, GC::copy_if(b1, copy, extra_b1));
    }
  }

  static node* range_root(node* b, const K& key_left, const K& key_right) {
    while (b) {
      if (Entry::comp(key_right, get_key(b))) { b = b->lc; continue; }
      if (Entry::comp(get_key(b), key_left)) { b = b->rc; continue; }
      break;
    }
    return b;
  }

  static node* left(node* b, const K& e) {
    if (!b) return NULL;
    if (Entry::comp(e, get_key(b))) return left(b->lc, e);
    node* r = Seq::make_node(Seq::get_entry(b));
    GC::increment(b->lc);
    return Seq::node_join(b->lc, left(b->rc, e), r);
  }

  static node* right(node* b, const K& e) {
    if (!b) return NULL;
    if (Entry::comp(get_key(b), e)) return right(b->rc, e);
    node* r = Seq::make_node(Seq::get_entry(b));
    GC::increment(b->rc);
    return Seq::node_join(right(b->lc, e), b->rc, r);
  }

  static node* range(node* b, const K& low, const K& high) {
    node* r = range_root(b, low, high);
    if (!r) return NULL;
    node* rr = Seq::make_node(Seq::get_entry(r));
    return Seq::node_join(right(r->lc, low), left(r->rc, high), rr);
  }

  static node* left_number(node* b, size_t rg) {
	  if (!b) return NULL;
	  if (rg == 0) return NULL;
	  if (Seq::size(b) <= rg) {
		  GC::increment(b); return b;
	  }//?
	  size_t x = Seq::size(b->lc);
	  if (x < rg) {
		  node* rr = Seq::make_node(Seq::get_entry(b));
		  GC::increment(b->lc);
		  node* rtree = left_number(b->rc, rg-x-1);
		  return Seq::node_join(b->lc, rtree, rr);
	  } else {
		  return left_number(b->lc, rg);
	  }
  }

  static node* range_num(node* b, const K& low, size_t rg) {
    while (b) {
      if (Entry::comp(get_key(b), low)) { b = b->rc; continue; }
      break;
    }
	if (!b) return NULL;
	node* ltree = range_num(b->lc, low, rg);
	size_t x = Seq::size(ltree);
	if (x == rg) return ltree;
	node* rtree = left_number(b->rc, rg-x-1);
	node* rr = Seq::make_node(Seq::get_entry(b));
	return Seq::node_join(ltree, rtree, rr);;
  }

  template<class Map, class Reduce>
  static std::pair<typename Reduce::T, size_t> left_number_mr(node* b, size_t rg, const Map& mp, const Reduce& rdc,
													size_t grain=utils::node_limit) {
	  using T = typename Reduce::T;
	  if (!b) return std::make_pair(rdc.identity(), 0);
	  if (rg == 0) return std::make_pair(rdc.identity(), 0);
	  size_t sz = Seq::size(b);
	  if (sz <= rg) {
		  //T res = Seq::map_reducef(b, mp, rdc, I);
		  T res = Seq::template map_reduce<Reduce>(b, mp, rdc, grain);
		  return std::make_pair(res, sz);
	  }
	  size_t x = Seq::size(b->lc);
	  if (x < rg) {
		  //T lres = Seq::map_reducef(b->lc, mp, rdc, I);
		  T lres = Seq::template map_reduce<Reduce>(b->lc, mp, rdc, grain);
		  T res = mp(Seq::get_entry(b));
		  res = rdc.add(lres, res);
		  std::pair<T, size_t> rres = left_number_mr(b->rc, rg-x-1, mp, rdc);
		  size_t s = x+1+rres.second;
		  res = rdc.add(res, rres.first);
		  return std::make_pair(res, s);
	  } else {
		  return left_number_mr(b->lc, rg, mp, rdc);
	  }
  }

  template<class Map, class Reduce>
  static std::pair<typename Reduce::T, size_t> range_num_mr(node* b, const K& low, size_t rg, const Map& mp, const Reduce& rdc,
												size_t grain=utils::node_limit) {
	using T = typename Reduce::T;
    while (b) {
      if (Entry::comp(get_key(b), low)) { b = b->rc; continue; }
      break;
    }
	if (!b) return std::make_pair(rdc.identity(), 0);
	using ptype = std::pair<T, size_t>;

	ptype ltree = range_num_mr(b->lc, low, rg, mp, rdc);
	size_t x = ltree.second;
	if (x == rg) return ltree;
	ptype rtree = left_number_mr(b->rc, rg-x-1, mp, rdc);
	T rr = mp(Seq::get_entry(b));
	size_t s = x + 1 + rtree.second;
	rr = rdc.add(ltree.first, rr);
	rr = rdc.add(rr, rtree.first);
	return std::make_pair(rr, s);
  }


  template <class Func, class J>
  static node* insert_j(node* b, const ET& e, const Func& f, const J& join,
			bool extra_ptr=false){
    if (!b) return Seq::single(e);
    bool copy = extra_ptr || (b->ref_cnt > 1);

    if (Entry::comp(get_key(b), Entry::get_key(e))) {
      node* l = GC::inc_if(b->lc, copy);
      node* r = insert_j(b->rc, e, f, join, copy);
      node* o = GC::copy_if(b, copy, extra_ptr);
      return join(l, r, o);
    }
    else if (Entry::comp(Entry::get_key(e), get_key(b))) {
      node* l = insert_j(b->lc, e, f, join, copy);
      node* r = GC::inc_if(b->rc, copy);
      node* o = GC::copy_if(b, copy, extra_ptr);
      return join(l, r, o);
    }
    else {
      node* l = GC::inc_if(b->lc, copy);
      node* r = GC::inc_if(b->rc, copy);
      node* o = GC::copy_if(b, copy, extra_ptr);
	  ET be = Seq::get_entry(b);
	  Seq::set_entry(o, e);
      //combine_values(o, e, false, f);
	  combine_values(o, be, true, f);
      return join(l, r, o);
    }
  }

  template <class Func>
  static node* insert(node* b, const ET& e, const Func& f, bool extra_ptr=false){
    auto join = [] (node* l, node* r, node* m) {return Seq::node_join(l,r,m);};
    return insert_j(b, e, f, join, extra_ptr);
  }

  template <class Func, class J>
  static node* update_j(node* b, const K& k, const Func& f, const J& join,
			bool extra_ptr=false){
    if (!b) return b;
    bool copy = extra_ptr || (b->ref_cnt > 1);

    if (Entry::comp(get_key(b), k)) {
      node* l = GC::inc_if(b->lc, copy);
      node* r = update_j(b->rc, k, f, join, copy);
      node* o = GC::copy_if(b, copy, extra_ptr);
      return join(l, r, o);
    }
    else if (Entry::comp(k, get_key(b))) {
      node* l = update_j(b->lc, k, f, join, copy);
      node* r = GC::inc_if(b->rc, copy);
      node* o = GC::copy_if(b, copy, extra_ptr);
      return join(l, r, o);
    }
    else {
      node* l = GC::inc_if(b->lc, copy);
      node* r = GC::inc_if(b->rc, copy);
      node* o = GC::copy_if(b, copy, extra_ptr);

	  update_value(o, f);
      return join(l, r, o);
    }
  }

  template <class Func>
  static node* update(node* b, const K& k, const Func& f, bool extra_ptr=false){
    auto join = [] (node* l, node* r, node* m) {return Seq::node_join(l,r,m);};
    return update_j(b, k, f, join, extra_ptr);
  }

  static node* deletet(node* b, const K& k, bool extra_ptr = false) {
    if (!b) return Seq::empty();

    bool copy = extra_ptr || (b->ref_cnt > 1);
    if (Entry::comp(get_key(b), k)) {
      node* l = GC::inc_if(b->lc, copy);
      node* r = deletet(b->rc, k, copy);
      return Seq::node_join(l, r, GC::copy_if(b, copy, extra_ptr));
    }
    else if (Entry::comp(k, get_key(b))) {
      node* r = GC::inc_if(b->rc, copy);
      node* l = deletet(b->lc, k, copy);
      return Seq::node_join(l, r, GC::copy_if(b, copy, extra_ptr));
    }
    else {
      node* o = Seq::join2_i(b->lc, b->rc, copy, copy);
      GC::dec_if(b, copy, extra_ptr);
      return o;
    }
  }
  template <typename F, typename AT>
  static size_t PAM_linear_search(AT* A, size_t size, const F& less) {
    for (size_t i = 0; i < size; i++)
      if (!less(A[i])) return i;
    return size;
  }

  // ?? why does PAM have its own binary search?
  template <typename F, typename AT>
  static size_t PAM_binary_search(AT* A, size_t n, const F& less) {
    size_t start = 0;
    size_t end = n;
    while (end-start > parlay::internal::_binary_search_base) {
      size_t mid = (end+start)/2;
      if (!less(A[mid])) end = mid;
      else start = mid + 1;
    }
    size_t x = start + PAM_linear_search(A+start, end-start, less);
	return x;
  }

  static node* multi_delete_sorted(node* b, K* A, size_t n,
				   bool extra_ptr = false) {
    if (!b) {
		return NULL;
	}
    if (n == 0) return GC::inc_if(b, extra_ptr);

    bool copy = extra_ptr || (b->ref_cnt > 1);
    K bk = get_key(b);

    auto less_val = [&] (K& a) -> bool {return Entry::comp(a,bk);};

	//---
	//parlay::sequence<ET> stemp(A, n);

    //size_t mid = parlay::binary_search(stemp, less_val);
	size_t mid = PAM_binary_search(A, n, less_val);

    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));

    auto P = utils::fork<node*>(utils::do_parallel(Seq::size(b), n),
	       [&] () {return multi_delete_sorted(b->lc, A, mid, copy);},
	       [&] () {return multi_delete_sorted(b->rc, A+mid+dup,
						  n-mid-dup, copy);});

    if (!dup) {
	  node* r = GC::copy_if(b, copy, extra_ptr);
      return Seq::node_join(P.first, P.second, r);
	} else {
	  GC::dec_if(b, copy, extra_ptr);
	  return Seq::join2(P.first, P.second);
	}
  }

  // assumes array A is of length n and is sorted with no duplicates
  template <class BinaryOp>
  static node* multi_insert_sorted(node* b, ET* A, size_t n,
				   const BinaryOp& op, bool extra_ptr = false) {
    if (!b) {
      node* x = Seq::from_array(A,n);
      return x;
    }
    if (n == 0) return GC::inc_if(b, extra_ptr);

    bool copy = extra_ptr || (b->ref_cnt > 1);
    K bk = get_key(b);

    auto less_val = [&] (ET& a) -> bool {return Entry::comp(Entry::get_key(a),bk);};

    	//---
    	//parlay::sequence<ET> stemp(A, n);

    //size_t mid = parlay::binary_search(stemp, less_val);
    	size_t mid = PAM_binary_search(A, n, less_val);

    bool dup = (mid < n) && (!Entry::comp(bk, Entry::get_key(A[mid])));

    auto P = utils::fork<node*>(utils::do_parallel(Seq::size(b), n),
    	       [&] () {return multi_insert_sorted(b->lc, A, mid, op, copy);},
    	       [&] () {return multi_insert_sorted(b->rc, A+mid+dup,
    						  n-mid-dup, op, copy);});

    node* r = GC::copy_if(b, copy, extra_ptr);
    if (dup) combine_values(r, A[mid], false, op);
    return Seq::node_join(P.first, P.second, r);
  }

  // ?? why are we using pointers
  template <class VE, class BinaryOp>
  static node* multi_update_sorted(node* b, std::pair<K, VE>* A, size_t n,
				   const BinaryOp& op, bool extra_ptr = false) {
	if (!b) return NULL;
    if (n == 0) return GC::inc_if(b, extra_ptr);
    bool copy = extra_ptr || (b->ref_cnt > 1);
    K bk = get_key(b);
    auto less_val = [&] (std::pair<K, VE>& a) -> bool {return Entry::comp(a.first,bk);};
    auto AS = parlay::make_slice(A,A+n);
    size_t mid = parlay::internal::binary_search(AS, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid].first));

    auto P = utils::fork<node*>(utils::do_parallel(Seq::size(b), n),
	       [&] () {return multi_update_sorted(b->lc, A, mid, op, copy);},
	       [&] () {return multi_update_sorted(b->rc, A+mid+dup,
						  n-mid-dup, op, copy);});

    node* r = GC::copy_if(b, copy, extra_ptr);
    if (dup) update_valuev(r, A[mid].second, op);
    return Seq::node_join(P.first, P.second, r);
  }

  static bool multi_find_sorted(node* b, K* A, size_t n, V* ret, size_t offset) {
    if (!b) return true;
    if (n == 0) return true;
    K bk = get_key(b);
    auto less_val = [&] (K a) -> bool {return Entry::comp(a,bk);};
    size_t mid = parlay::internal::binary_search(parlay::make_slice(A, A+n), less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));
    utils::fork<bool>(utils::do_parallel(Seq::size(b), n),
	       [&] () {return multi_find_sorted(b->lc, A, mid, ret, offset);},
	       [&] () {return multi_find_sorted(b->rc, A+mid+dup,
						  n-mid-dup, ret, offset+mid+dup);});

    if (dup) ret[offset+mid] = get_val(b);
	return true;
  }

  template<class InTree, class Func>
  static node* map(typename InTree::node* b, const Func& f) {
    auto g = [&] (typename InTree::ET& a) {
      return ET(InTree::Entry::get_key(a),f(a));};
    return Seq::template map<InTree>(b, g);
  }

  template<class InTree, class Func>
  static node* map_set(typename InTree::node* b, const Func& f) {
    return Seq::template map<InTree>(b, f);
  }

  template<class Seq1, class Func>
    static node* map_filter(typename Seq1::node* b, const Func& f,
			    size_t granularity=utils::node_limit) {
    auto g = [&] (typename Seq1::ET& a) {
      std::optional<V> v = f(a);
      if (v) return std::optional<ET>(ET(Seq1::Entry::get_key(a),*v));
      else return std::optional<ET>();
    };

    return Seq::template map_filter<Seq1>(b, g, granularity);
  }

};
