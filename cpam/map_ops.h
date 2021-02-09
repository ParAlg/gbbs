#pragma once
#include "utils.h"
#include "parlay/primitives.h"
#include "parlay/internal/binary_search.h"

// *******************************************
//   MAPS and SETS
// *******************************************

namespace cpam {

template<class Seq, class EntryT>
struct map_ops : Seq {
  using Entry = EntryT;
  using node = typename Seq::node;
  using regular_node = typename Seq::regular_node;
  using compressed_node = typename Seq::compressed_node;
  using ET = typename Seq::ET;
  using GC = typename Seq::GC;
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using _Seq = Seq;
  using ptr = typename GC::ptr;

  static bool comp(K a, K b) { return Entry::comp(a,b);}
  static K get_key(node *s) { return Entry::get_key(*Seq::get_entry_p(s));}
  static V get_val(node *s) { return Entry::get_val(*Seq::get_entry_p(s));}
  //static K get_key(node *s) { return Entry::get_key(Seq::get_entry(s));}
  //static V get_val(node *s) { return Entry::get_val(Seq::get_entry(s));}

  static std::optional<ET> find_compressed(ptr b, const K& key) {
    std::optional<ET> ret;
    auto find_f = [&] (const ET& et) {
      if (!(Entry::comp(key, Entry::get_key(et)) || Entry::comp(Entry::get_key(et), key))) {
        ret = et;
      }
    };
    Seq::iterate_seq(b.node_ptr(), find_f);
    // TODO: should probably use decode_cond or something like this.
    return ret;
  }

  static std::optional<ET> find(ptr b, const K& key) {
    if (b.empty()) return {};
    if (b.is_compressed()) {
      return find_compressed(std::move(b), key);
    }
    auto [lc, e, rc, m] = Seq::expose(std::move(b));
    std::optional<ET> ret = {};
    if (Entry::comp(key, Entry::get_key(e))) {
      ret = find(std::move(lc), key);
    } else if (Entry::comp(Entry::get_key(e), key)) {
      ret = find(std::move(rc), key);
    } else {
      ret = {e};
    }
    GC::decrement(m);
    return ret;
  }

// TODO
//  template <class BinaryOp>
//  static inline void update_value(node* a, const BinaryOp& op) {
//    ET re = Seq::get_entry(a);
//    Entry::set_val(re, op(re));
//    Seq::set_entry(a, re);
//  }

// TODO
//  static node* previous(node* b, const K& key) {
//    node* r = NULL;
//    while (b) {
//      if (Entry::comp(get_key(b), key)) {r = b; b = b->rc;}
//      else b = b->lc;
//    }
//    return r;
//  }
//
//  static node* next(node* b, const K& key) {
//    node* r = NULL;
//    while (b) {
//      if (Entry::comp(key, get_key(b)) ) {r = b; b = b->lc;}
//      else b = b->rc;
//    }
//    return r;
//  }

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
    node* l;
    std::optional<ET> mid;
    node* r;
    split_info(node* l, std::optional<ET> mid, node* r)
      : l(l), mid(mid), r(r) {};
  };

  static split_info split(ptr a, const K& k) {
    if (a.empty()) return split_info(NULL, std::nullopt, NULL);
    auto [lc, e, rc, m] = Seq::expose(std::move(a));
    const K& kmid = Entry::get_key(e);

    if (Entry::comp(kmid, k)) {
      split_info bstpair = split(std::move(rc), k);
      bstpair.l = Seq::join(lc.node_ptr(), e, bstpair.l, m);
      return bstpair;
    } else if (Entry::comp(k, kmid)) {
      split_info bstpair = split(std::move(lc), k);
      bstpair.r = Seq::join(bstpair.r, e, rc.node_ptr(), m);
      return bstpair;
    } else {
      GC::decrement(m);
      return split_info(lc.node_ptr(), e, rc.node_ptr());
    }
  }

// TODO
//  static split_info split_copy(node* bst, const K& e) {
//    if (!bst) return split_info(NULL, NULL, false);
//
//    else if (Entry::comp(get_key(bst), e)) {
//      // following should be a copy ????
//      node* join = Seq::make_node(Seq::get_entry(bst));
//      split_info bstpair = split_copy(bst->rc, e);
//      GC::increment(bst->lc);
//      bstpair.first = Seq::node_join(bst->lc, bstpair.first, join);
//      return bstpair;
//    }
//    else if (Entry::comp(e, get_key(bst))) {
//      node* join = Seq::make_node(Seq::get_entry(bst));
//      split_info bstpair = split_copy(bst->lc, e);
//      GC::increment(bst->rc);
//      bstpair.second = Seq::node_join(bstpair.second, bst->rc, join);
//      return bstpair;
//    }
//    else {
//      GC::increment(bst->lc); GC::increment(bst->rc);
//      split_info ret(bst->lc, bst->rc, true);
//      ret.entry = Seq::get_entry(bst);
//      return ret;
//    }
//  }

// TODO: deprecate?
//  // A version that will reuse a node if ref count is 1
//  // Will decrement ref count of root if it is copied
//  static split_info split_inplace(node* bst, const K& e) {
//    if (!bst) return split_info(NULL, NULL, false);
//    else if (bst->ref_cnt > 1) {
//		//std::cout << "rc>1" << std::endl;
//      split_info ret = split_copy(bst, e);
//      GC::decrement_recursive(bst);
//      return ret;
//    }
//    else if (Entry::comp(get_key(bst), e)) {
//      split_info bstpair = split_inplace(bst->rc, e);
//      bstpair.first = Seq::node_join(bst->lc, bstpair.first, bst);
//      return bstpair;
//    }
//    else if (Entry::comp(e, get_key(bst))) {
//      split_info bstpair = split_inplace(bst->lc, e);
//      bstpair.second = Seq::node_join(bstpair.second, bst->rc, bst);
//      return bstpair;
//    }
//    else {
//      split_info ret(bst->lc, bst->rc, true);
//      ret.entry = Seq::get_entry(bst);
//      GC::decrement(bst);
//      return ret;
//    }
//  }
//
//  static inline split_info split(node* bst, const K& e) {
//    return split_inplace(bst, e);
//  }

  template <class BinaryOp>
  static inline void combine_values(ET& re, ET e, bool reverse, const BinaryOp& op) {
    if (reverse) Entry::set_val(re, op(Entry::get_val(e), Entry::get_val(re)));
    else Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(e)));
  }

  template <class BinaryOp>
  static inline void combine_values(node* a, ET e, bool reverse, const BinaryOp& op) {
//    ET& re = Seq::get_entry(a);
//    combine_values(re, e, reverse, op);
    ET re = Seq::get_entry(a);
    if (reverse) Entry::set_val(re, op(Entry::get_val(e), Entry::get_val(re)));
    else Entry::set_val(re, op(Entry::get_val(re), Entry::get_val(e)));
    Seq::set_entry(a, re);
  }


// TODO
//  template <class VE, class BinaryOp>
//  static inline void update_valuev(node* a, VE v0, const BinaryOp& op) {
//    ET re = Seq::get_entry(a);
//    Entry::set_val(re, op(Entry::get_val(re), v0));
//    Seq::set_entry(a, re);
//  }

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

// TODO: remove
//  // Works in-place when possible.
//  // extra_b2 means there is an extra pointer to b2 not included
//  // in the reference count.  It is used as an optimization to reduce
//  // the number of ref_cnt updates.
//  template <class BinaryOp>
//  static node* uniont(node* b1, node* b2, const BinaryOp op,
//		      bool extra_b2 = false) {
//    if (!b1) return GC::inc_if(b2, extra_b2);
//    if (!b2) return b1;
//
//    size_t n1 = Seq::size(b1);   size_t n2 = Seq::size(b2);
//
//    bool copy = extra_b2 || (b2->ref_cnt > 1);
//    node* r = copy ? Seq::make_node(Seq::get_entry(b2)) : b2;
//    split_info bsts = split(b1, get_key(b2));
//
//    auto P = utils::fork<node*>(utils::do_parallel(n1, n2),
//      [&] () {return uniont(bsts.first, b2->lc, op, copy);},
//      [&] () {return uniont(bsts.second, b2->rc, op, copy);});
//
//    if (copy && !extra_b2) GC::decrement_recursive(b2);
//    if (bsts.removed) combine_values(r, bsts.entry, true, op);
//    return Seq::node_join(P.first, P.second, r);
//  }

// TODO: remove
//  // extra_b2 is for b2
//  template <class Seq1, class Seq2, class BinaryOp>
//  static node* intersect(typename Seq1::node* b1, typename Seq2::node* b2,
//			 const BinaryOp& op,
//			 bool extra_b2 = false) {
//	//std::cout << "size: " << n1 << " " << n2 << std::endl;
//    if (!b1) {if (!extra_b2) Seq2::GC::decrement_recursive(b2); return NULL;}
//    if (!b2) {Seq1::GC::decrement_recursive(b1); return NULL;}
//    size_t n1 = Seq1::size(b1);   size_t n2 = Seq2::size(b2);
//    bool copy = extra_b2 || (b2->ref_cnt > 1);
//    typename Seq1::split_info bsts = Seq1::split(b1, Seq2::get_key(b2));
//
//    auto P = utils::fork<node*>(utils::do_parallel(n1, n2),
//	[&]() {return intersect<Seq1,Seq2>(bsts.first, b2->lc, op, copy);},
//	[&]() {return intersect<Seq1,Seq2>(bsts.second, b2->rc, op, copy);}
//    );
//
//    if (bsts.removed) {
//      ET e(Seq2::get_key(b2),
//	   op(Seq1::Entry::get_val(bsts.entry),
//	      Seq2::Entry::get_val(Seq2::get_entry(b2))));
//      node* r = Seq::make_node(e);
//      //if (copy && !extra_b2) Seq2::GC::decrement_recursive(b2);
//      Seq2::GC::dec_if(b2, copy, extra_b2);
//      return Seq::node_join(P.first, P.second, r);
//    } else {
//      //if (copy && !extra_b2) Seq2::GC::decrement_recursive(b2);
//      Seq2::GC::dec_if(b2, copy, extra_b2);
//      return Seq::join2(P.first, P.second);
//    }
//  }

  // reuses node
  template <class BinaryOp>
  static node* uniont(ptr b1, ptr b2, const BinaryOp& op) {
    if (b1.empty()) return b2.node_ptr();
    if (b2.empty()) return b1.node_ptr();

    size_t n1 = b1.size(), n2 = b2.size();
    if (n1 + n2 <= utils::kBaseCaseSize) {
      return union_bc(std::move(b1), std::move(b2), op);
    }

    auto [l2, e, r2, m2] = Seq::expose(std::move(b2));
    if (!m2) m2 = Seq::single(e);

    auto sp1 = split(std::move(b1), Entry::get_key(e));
#ifdef DEBUG
    Seq::check_structure(sp1.l); Seq::check_structure(sp1.r);
#endif

    auto [l, r] = utils::fork<node*>(utils::do_parallel(n1, n2),
      [&] () {return uniont(sp1.l, std::move(l2), op);},
      [&] () {return uniont(sp1.r, std::move(r2), op);});
#ifdef DEBUG
    Seq::check_structure(l); Seq::check_structure(r);
#endif

    if (sp1.mid) combine_values(m2, *sp1.mid, true, op);
    return Seq::node_join(l, r, m2);
  }


  template <class Seq1, class Seq2, class BinaryOp>
  static node* intersect(typename Seq1::node* b1, typename Seq2::ptr b2,
                         const BinaryOp& op) {
    if (!b1) return NULL;
    if (b2.empty()) {Seq1::GC::decrement_recursive(b1); return NULL;}

    size_t n1 = Seq1::size(b1);
    size_t n2 = b2.size();

    if (n1 + n2 <= utils::kBaseCaseSize) {
      return intersect_bc(std::move(b1), std::move(b2), op);
    }

    auto [l2, e2, r2, m2] = Seq2::expose(std::move(b2));
    auto key = Seq2::Entry::get_key(e2);
    auto sp1 = Seq1::split(b1, key);

#ifdef DEBUG
    Seq::check_structure(sp1.l); Seq::check_structure(sp1.r);
#endif

    auto [l, r] = utils::fork<node*>(utils::do_parallel(n1, n2),
        [&]() {return intersect<Seq1,Seq2>(sp1.l, std::move(l2), op);},
        [&]() {return intersect<Seq1,Seq2>(sp1.r, std::move(r2), op);}
    );

#ifdef DEBUG
    Seq::check_structure(l); Seq::check_structure(r);
#endif

    Seq2::GC::decrement(m2);
    if (sp1.mid) {
      ET e(key,
           op(Seq1::Entry::get_val(*sp1.mid),
              Seq2::Entry::get_val(e)));
      return Seq::join(l, e, r, nullptr);
    } else return Seq::join2(l, r);
  }


  static node* difference(ptr b1, node* b2) {
    if (b1.empty()) {GC::decrement_recursive(b2); return NULL;}
    if (!b2) return b1.node_ptr();

    size_t n1 = b1.size(), n2 = Seq::size(b2);
    if (n1 + n2 <= utils::kBaseCaseSize) {
      return difference_bc(std::move(b1), std::move(b2));
    }

    auto [l1, e, r1, m1] = Seq::expose(std::move(b1));
    auto sp2 = split(b2, Entry::get_key(e));

    auto [l, r] = utils::fork<node*>(utils::do_parallel(n1, n2),
        [&]() {return difference(std::move(l1), sp2.l);},
        [&]() {return difference(std::move(r1), sp2.r);});

    if (sp2.mid) {
      GC::decrement(m1);
      return Seq::join2(l, r);
    } else {
      if (!m1) m1 = Seq::single(e);
      return Seq::node_join(l, r, m1);
    }
  }


// TODO
//  static node* range_root(node* b, const K& key_left, const K& key_right) {
//    while (b) {
//      if (Entry::comp(key_right, get_key(b))) { b = b->lc; continue; }
//      if (Entry::comp(get_key(b), key_left)) { b = b->rc; continue; }
//      break;
//    }
//    return b;
//  }

// TODO
//  static node* left(node* b, const K& e) {
//    if (!b) return NULL;
//    if (Entry::comp(e, get_key(b))) return left(b->lc, e);
//    node* r = Seq::make_node(Seq::get_entry(b));
//    GC::increment(b->lc);
//    return Seq::node_join(b->lc, left(b->rc, e), r);
//  }

// TODO
//  static node* right(node* b, const K& e) {
//    if (!b) return NULL;
//    if (Entry::comp(get_key(b), e)) return right(b->rc, e);
//    node* r = Seq::make_node(Seq::get_entry(b));
//    GC::increment(b->rc);
//    return Seq::node_join(right(b->lc, e), b->rc, r);
//  }

// TODO
//  static node* range(node* b, const K& low, const K& high) {
//    node* r = range_root(b, low, high);
//    if (!r) return NULL;
//    node* rr = Seq::make_node(Seq::get_entry(r));
//    return Seq::node_join(right(r->lc, low), left(r->rc, high), rr);
//  }

// TODO
//  static node* left_number(node* b, size_t rg) {
//	  if (!b) return NULL;
//	  if (rg == 0) return NULL;
//	  if (Seq::size(b) <= rg) {
//		  GC::increment(b); return b;
//	  }//?
//	  size_t x = Seq::size(b->lc);
//	  if (x < rg) {
//		  node* rr = Seq::make_node(Seq::get_entry(b));
//		  GC::increment(b->lc);
//		  node* rtree = left_number(b->rc, rg-x-1);
//		  return Seq::node_join(b->lc, rtree, rr);
//	  } else {
//		  return left_number(b->lc, rg);
//	  }
//  }

// TODO
//  static node* range_num(node* b, const K& low, size_t rg) {
//    while (b) {
//      if (Entry::comp(get_key(b), low)) { b = b->rc; continue; }
//      break;
//    }
//	if (!b) return NULL;
//	node* ltree = range_num(b->lc, low, rg);
//	size_t x = Seq::size(ltree);
//	if (x == rg) return ltree;
//	node* rtree = left_number(b->rc, rg-x-1);
//	node* rr = Seq::make_node(Seq::get_entry(b));
//	return Seq::node_join(ltree, rtree, rr);;
//  }

// TODO
//  template<class Map, class Reduce>
//  static std::pair<typename Reduce::T, size_t> left_number_mr(node* b, size_t rg, const Map& mp, const Reduce& rdc,
//													size_t grain=utils::node_limit) {
//	  using T = typename Reduce::T;
//	  if (!b) return std::make_pair(rdc.identity(), 0);
//	  if (rg == 0) return std::make_pair(rdc.identity(), 0);
//	  size_t sz = Seq::size(b);
//	  if (sz <= rg) {
//		  //T res = Seq::map_reducef(b, mp, rdc, I);
//		  T res = Seq::template map_reduce<Reduce>(b, mp, rdc, grain);
//		  return std::make_pair(res, sz);
//	  }
//	  size_t x = Seq::size(b->lc);
//	  if (x < rg) {
//		  //T lres = Seq::map_reducef(b->lc, mp, rdc, I);
//		  T lres = Seq::template map_reduce<Reduce>(b->lc, mp, rdc, grain);
//		  T res = mp(Seq::get_entry(b));
//		  res = rdc.add(lres, res);
//		  std::pair<T, size_t> rres = left_number_mr(b->rc, rg-x-1, mp, rdc);
//		  size_t s = x+1+rres.second;
//		  res = rdc.add(res, rres.first);
//		  return std::make_pair(res, s);
//	  } else {
//		  return left_number_mr(b->lc, rg, mp, rdc);
//	  }
//  }

// TODO
//  template<class Map, class Reduce>
//  static std::pair<typename Reduce::T, size_t> range_num_mr(node* b, const K& low, size_t rg, const Map& mp, const Reduce& rdc,
//												size_t grain=utils::node_limit) {
//	using T = typename Reduce::T;
//    while (b) {
//      if (Entry::comp(get_key(b), low)) { b = b->rc; continue; }
//      break;
//    }
//	if (!b) return std::make_pair(rdc.identity(), 0);
//	using ptype = std::pair<T, size_t>;
//
//	ptype ltree = range_num_mr(b->lc, low, rg, mp, rdc);
//	size_t x = ltree.second;
//	if (x == rg) return ltree;
//	ptype rtree = left_number_mr(b->rc, rg-x-1, mp, rdc);
//	T rr = mp(Seq::get_entry(b));
//	size_t s = x + 1 + rtree.second;
//	rr = rdc.add(ltree.first, rr);
//	rr = rdc.add(rr, rtree.first);
//	return std::make_pair(rr, s);
//  }

// TODO: remove
//  template <class Func, class J>
//  static node* insert_j(node* b, const ET& e, const Func& f, const J& join,
//			bool extra_ptr=false){
//    if (!b) return Seq::single(e);
//    bool copy = extra_ptr || (b->ref_cnt > 1);
//
//    if (Entry::comp(get_key(b), Entry::get_key(e))) {
//      node* l = GC::inc_if(b->lc, copy);
//      node* r = insert_j(b->rc, e, f, join, copy);
//      node* o = GC::copy_if(b, copy, extra_ptr);
//      return join(l, r, o);
//    }
//    else if (Entry::comp(Entry::get_key(e), get_key(b))) {
//      node* l = insert_j(b->lc, e, f, join, copy);
//      node* r = GC::inc_if(b->rc, copy);
//      node* o = GC::copy_if(b, copy, extra_ptr);
//      return join(l, r, o);
//    }
//    else {
//      node* l = GC::inc_if(b->lc, copy);
//      node* r = GC::inc_if(b->rc, copy);
//      node* o = GC::copy_if(b, copy, extra_ptr);
//	  ET be = Seq::get_entry(b);
//	  Seq::set_entry(o, e);
//      //combine_values(o, e, false, f);
//	  combine_values(o, be, true, f);
//      return join(l, r, o);
//    }
//  }


  template <bool copy=false>
  static node* inc_tmpl(node* x) {
    if constexpr (copy) {
      GC::increment(x);
    }
    return x;
  }

  // TODO: move this to a better location. Doesn't belong in map_ops.
  template <bool copy=false>
  static node* make_node_tmpl(node* x) {
    if constexpr (copy) {
      return Seq::make_regular_node(Seq::get_entry(x));
    }
    return x;
  }

  template <class Func, class J>
  static node* insert_compressed(compressed_node* b, const ET& e, const Func& f) {
    ET arr[2*B + 1];
    Seq::compressed_node_elms(b, arr);
    size_t n = Seq::size(b);
    assert(n <= 2*B);
    ET merged[2*B + 1];
    size_t out_off = 0;
    size_t k = 0;
    K key = Entry::get_key(e);
    bool placed = false;
    while (k < n) {
      if (Entry::comp(Entry::get_key(arr[k]), key)) {
        parlay::move_uninitialized(merged[out_off++], arr[k++]);
      } else if (Entry::comp(key, Entry::get_key(arr[k]))) {
        parlay::assign_uninitialized(merged[out_off++], e);
        placed = true;
        break;
      } else {  // arr[k] == key
        parlay::assign_uninitialized(merged[out_off], arr[k++]);
        combine_values(merged[out_off], e, true, f);
        out_off++;
        placed = true;
        break;
      }
    }
    while (k < n) {
      parlay::move_uninitialized(merged[out_off++], arr[k++]);}
    if (!placed) {
      parlay::assign_uninitialized(merged[out_off++], e);
    }
    return Seq::make_compressed(merged, out_off);
  }

  // A specialized version of insert that will switch to the version with
  // copy=true once it hits a node with ref_cnt(node) > 1.
  template <class Func, class J, bool copy=false>
  static node* insert_tmpl(node* b, const ET& e, const Func& f, const J& join) {
    if (!b) return Seq::single(e);

    if constexpr (!copy) {
      if (Seq::ref_cnt(b) > 1) {
        auto r = insert_tmpl<Func, J, true>(b, e, f, join);
        GC::decrement_recursive(b);
        return r;
      }
    }
    if (Seq::is_compressed(b)) {
      auto r = insert_compressed<Func, J>(Seq::cast_to_compressed(b), e, f);
      if constexpr (!copy) {
        GC::decrement_recursive(b);}
      return r;
    }

    regular_node* br = Seq::cast_to_regular(b);
    regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(b));
    if (Entry::comp(get_key(br), Entry::get_key(e))) {
      node* r = insert_tmpl<Func, J, copy>(br->rc, e, f, join);
      return join(inc_tmpl<copy>(br->lc), r, o);
    } else if (Entry::comp(Entry::get_key(e), get_key(br))) {
      node* l = insert_tmpl<Func, J, copy>(br->lc, e, f, join);
      return join(l, inc_tmpl<copy>(br->rc), o);
    } else {
      ET be = Seq::get_entry(br);
      Seq::set_entry(o, e);
      combine_values(o, be, true, f);
      return join(inc_tmpl<copy>(br->lc), inc_tmpl<copy>(br->rc), o);
    }
  }

  template <class Func>
  static node* insert(node* b, const ET& e, const Func& f) {
    auto join = [] (node* l, node* r, regular_node* m) {
      auto ret = Seq::node_join(l,r,m);
      return ret;
    };
    return insert_tmpl<Func, decltype(join), false>(b, e, f, join);
  }

  static node* remove_compressed(compressed_node* b, const K& key) {
    ET arr[2*B];
    Seq::compressed_node_elms(b, arr);
    size_t n = Seq::size(b);
    assert(n <= 2*B);
    ET merged[2*B];
    size_t k = 0;
    size_t out_off = 0;
    while (k < n) {
      if (Entry::comp(Entry::get_key(arr[k]), key)) {
        parlay::move_uninitialized(merged[out_off++], arr[k++]);
      } else if (Entry::comp(key, Entry::get_key(arr[k]))) {
        break;
      } else {  // arr[k] == key
        arr[k].~ET();
        k++;
        out_off++;
        break;
      }
    }
    while (k < n) {
      parlay::move_uninitialized(merged[out_off++], arr[k++]);}
    return Seq::make_compressed(merged, out_off);
  }

  // A specialized version of insert that will switch to the version with
  // copy=true once it hits a node with ref_cnt(node) > 1.
  template <class J, bool copy=false>
  static node* remove_tmpl(node* b, const K& k, const J& join) {
    if (!b) return Seq::empty();

    if constexpr (!copy) {
      if (Seq::ref_cnt(b) > 1) {
        auto r = remove_tmpl<J, true>(b, k, join);
        GC::decrement_recursive(b);
        return r;
      }
    }

    if (Seq::is_compressed(b)) {
      auto r = remove_compressed(Seq::cast_to_compressed(b), k);
      if constexpr (!copy) {
        GC::decrement_recursive(b);}
      return r;
    }

    regular_node* br = Seq::cast_to_regular(b);
    if (Entry::comp(get_key(br), k)) {
      regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(b));
      node* r = remove_tmpl<J, copy>(br->rc, k, join);
      return join(inc_tmpl<copy>(br->lc), r, o);
    } else if (Entry::comp(k, get_key(br))) {
      regular_node* o = Seq::cast_to_regular(make_node_tmpl<copy>(b));
      node* l = remove_tmpl<J, copy>(br->lc, k, join);
      return join(l, inc_tmpl<copy>(br->rc), o);
    } else {
      auto lc = br->lc, rc = br->rc;
      GC::decrement(br);
      return Seq::join2(inc_tmpl<copy>(lc), inc_tmpl<copy>(rc));
    }
  }

  static node* deletet(node* b, const K& k) {
    auto join = [] (node* l, node* r, regular_node* m) {
      auto ret = Seq::node_join(l,r,m);
      return ret;
    };
    return remove_tmpl<decltype(join), false>(b, k, join);
  }



// TODO
//  template <class Func, class J>
//  static node* update_j(node* b, const K& k, const Func& f, const J& join,
//			bool extra_ptr=false){
//    if (!b) return b;
//    bool copy = extra_ptr || (b->ref_cnt > 1);
//
//    if (Entry::comp(get_key(b), k)) {
//      node* l = GC::inc_if(b->lc, copy);
//      node* r = update_j(b->rc, k, f, join, copy);
//      node* o = GC::copy_if(b, copy, extra_ptr);
//      return join(l, r, o);
//    }
//    else if (Entry::comp(k, get_key(b))) {
//      node* l = update_j(b->lc, k, f, join, copy);
//      node* r = GC::inc_if(b->rc, copy);
//      node* o = GC::copy_if(b, copy, extra_ptr);
//      return join(l, r, o);
//    }
//    else {
//      node* l = GC::inc_if(b->lc, copy);
//      node* r = GC::inc_if(b->rc, copy);
//      node* o = GC::copy_if(b, copy, extra_ptr);
//
//	  update_value(o, f);
//      return join(l, r, o);
//    }
//  }

// TODO
//  template <class Func>
//  static node* update(node* b, const K& k, const Func& f, bool extra_ptr=false){
//    auto join = [] (node* l, node* r, node* m) {return Seq::node_join(l,r,m);};
//    return update_j(b, k, f, join, extra_ptr);
//  }

// TODO
//  static node* deletet(node* b, const K& k, bool extra_ptr = false) {
//    if (!b) return Seq::empty();
//
//    bool copy = extra_ptr || (b->ref_cnt > 1);
//    if (Entry::comp(get_key(b), k)) {
//      node* l = GC::inc_if(b->lc, copy);
//      node* r = deletet(b->rc, k, copy);
//      return Seq::node_join(l, r, GC::copy_if(b, copy, extra_ptr));
//    }
//    else if (Entry::comp(k, get_key(b))) {
//      node* r = GC::inc_if(b->rc, copy);
//      node* l = deletet(b->lc, k, copy);
//      return Seq::node_join(l, r, GC::copy_if(b, copy, extra_ptr));
//    }
//    else {
//      node* o = Seq::join2_i(b->lc, b->rc, copy, copy);
//      GC::dec_if(b, copy, extra_ptr);
//      return o;
//    }
//  }

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

// TODO
//  static node* multi_delete_sorted(node* b, K* A, size_t n,
//				   bool extra_ptr = false) {
//    if (!b) {
//		return NULL;
//	}
//    if (n == 0) return GC::inc_if(b, extra_ptr);
//
//    bool copy = extra_ptr || (b->ref_cnt > 1);
//    K bk = get_key(b);
//
//    auto less_val = [&] (K& a) -> bool {return Entry::comp(a,bk);};
//
//	//---
//	//parlay::sequence<ET> stemp(A, n);
//
//    //size_t mid = parlay::binary_search(stemp, less_val);
//	size_t mid = PAM_binary_search(A, n, less_val);
//
//    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));
//
//    auto P = utils::fork<node*>(utils::do_parallel(Seq::size(b), n),
//	       [&] () {return multi_delete_sorted(b->lc, A, mid, copy);},
//	       [&] () {return multi_delete_sorted(b->rc, A+mid+dup,
//						  n-mid-dup, copy);});
//
//    if (!dup) {
//	  node* r = GC::copy_if(b, copy, extra_ptr);
//      return Seq::node_join(P.first, P.second, r);
//	} else {
//	  GC::dec_if(b, copy, extra_ptr);
//	  return Seq::join2(P.first, P.second);
//	}
//  }

  // assumes array A is of length n and is sorted with no duplicates
  template <class BinaryOp>
  static node* multi_insert_sorted(ptr b, ET* A, size_t n,
				   const BinaryOp& op) {
    if (b.empty()) {
      node* x = Seq::from_array(A,n);
      return x;
    }
    if (n == 0) return b.node_ptr();

    size_t tot = b.size() + n;
    if (tot <= utils::kBaseCaseSize) {
      return multiinsert_bc(std::move(b), A, n, op);
    }

    auto [lc, e, rc, root] = Seq::expose(std::move(b));
    if (!root) root = Seq::single(e);

    K bk = Entry::get_key(e);
    auto less_val = [&](ET& a) -> bool {
      return Entry::comp(Entry::get_key(a), bk);
    };

    size_t mid = PAM_binary_search(A, n, less_val);
    bool dup = (mid < n) && (!Entry::comp(bk, Entry::get_key(A[mid])));

    auto P = utils::fork<node*>(utils::do_parallel(b.size(), n),
        [&] () {return multi_insert_sorted(std::move(lc), A, mid, op);},
        [&] () {return multi_insert_sorted(std::move(rc), A+mid+dup, n-mid-dup, op);});

    if (dup) combine_values(root, A[mid], false, op);
    return Seq::node_join(P.first, P.second, root);
  }

// TODO
//  template <class VE, class BinaryOp>
//  static node* multi_update_sorted(node* b, std::pair<K, VE>* A, size_t n,
//				   const BinaryOp& op, bool extra_ptr = false) {
//	if (!b) return NULL;
//    if (n == 0) return GC::inc_if(b, extra_ptr);
//    bool copy = extra_ptr || (b->ref_cnt > 1);
//    K bk = get_key(b);
//    auto less_val = [&] (std::pair<K, VE>& a) -> bool {return Entry::comp(a.first,bk);};
//    auto AS = parlay::make_slice(A,A+n);
//    size_t mid = parlay::internal::binary_search(AS, less_val);
//    bool dup = (mid < n) && (!Entry::comp(bk, A[mid].first));
//
//    auto P = utils::fork<node*>(utils::do_parallel(Seq::size(b), n),
//	       [&] () {return multi_update_sorted(b->lc, A, mid, op, copy);},
//	       [&] () {return multi_update_sorted(b->rc, A+mid+dup,
//						  n-mid-dup, op, copy);});
//
//    node* r = GC::copy_if(b, copy, extra_ptr);
//    if (dup) update_valuev(r, A[mid].second, op);
//    return Seq::node_join(P.first, P.second, r);
//  }

// TODO
//  static bool multi_find_sorted(node* b, K* A, size_t n, V* ret, size_t offset) {
//    if (!b) return true;
//    if (n == 0) return true;
//    K bk = get_key(b);
//    auto less_val = [&] (K a) -> bool {return Entry::comp(a,bk);};
//    size_t mid = parlay::internal::binary_search(parlay::make_slice(A, A+n), less_val);
//    bool dup = (mid < n) && (!Entry::comp(bk, A[mid]));
//    utils::fork<bool>(utils::do_parallel(Seq::size(b), n),
//	       [&] () {return multi_find_sorted(b->lc, A, mid, ret, offset);},
//	       [&] () {return multi_find_sorted(b->rc, A+mid+dup,
//						  n-mid-dup, ret, offset+mid+dup);});
//
//    if (dup) ret[offset+mid] = get_val(b);
//	return true;
//  }

// TODO
//  template<class InTree, class Func>
//  static node* map(typename InTree::node* b, const Func& f) {
//    auto g = [&] (typename InTree::ET& a) {
//      return ET(InTree::Entry::get_key(a),f(a));};
//    return Seq::template map<InTree>(b, g);
//  }

// TODO
//  template<class InTree, class Func>
//  static node* map_set(typename InTree::node* b, const Func& f) {
//    return Seq::template map<InTree>(b, f);
//  }

// TODO
//  template<class Seq1, class Func>
//    static node* map_filter(typename Seq1::node* b, const Func& f,
//			    size_t granularity=utils::node_limit) {
//    auto g = [&] (typename Seq1::ET& a) {
//      std::optional<V> v = f(a);
//      if (v) return std::optional<ET>(ET(Seq1::Entry::get_key(a),*v));
//      else return std::optional<ET>();
//    };
//
//    return Seq::template map_filter<Seq1>(b, g, granularity);
//  }


  template <class BinaryOp>
  static node* multiinsert_bc(ptr b1, ET* A, size_t n, const BinaryOp& op) {
    auto n_b1 = b1.node_ptr();

    ET stack[utils::kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&] (const ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    for (size_t i=0; i<n; i++) {
      parlay::move_uninitialized(stack[offset++], A[i]);
    }
    assert(offset <= utils::kBaseCaseSize);

    Seq::decrement_recursive(n_b1);

    ET output[utils::kBaseCaseSize + 1];

    // TODO: refactor merge code to share between union_bc and this code.

    // merge
    size_t nA = offset_2; size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      const auto& k_a = Entry::get_key(stack[i]);
      const auto& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        parlay::move_uninitialized(output[out_off++], stack[j]);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, op(Entry::get_val(stack[j]), Entry::get_val(re)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      parlay::move_uninitialized(output[out_off++], stack[j]);
      j++;
    }

    // build tree
    if (out_off < utils::compression_block_size) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }


  template <class BinaryOp>
  static node* union_bc(ptr b1, ptr b2, const BinaryOp& op) {
    auto n_b1 = b1.node_ptr();
    auto n_b2 = b2.node_ptr();

    ET stack[utils::kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&] (const ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    Seq::iterate_seq(n_b2, copy_f);
    assert(offset <= utils::kBaseCaseSize);

    Seq::decrement_recursive(n_b1);
    Seq::decrement_recursive(n_b2);

    ET output[utils::kBaseCaseSize + 1];

    // merge
    size_t nA = offset_2; size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      const auto& k_a = Entry::get_key(stack[i]);
      const auto& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        parlay::move_uninitialized(output[out_off++], stack[j]);
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, op(Entry::get_val(stack[j]), Entry::get_val(re)));
        out_off++;
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }
    while (j < nB) {
      parlay::move_uninitialized(output[out_off++], stack[j]);
      j++;
    }

    // build tree
    if (out_off < utils::compression_block_size) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

  static node* difference_bc(ptr b1, ptr b2) {
    auto n_b1 = b1.node_ptr();
    auto n_b2 = b2.node_ptr();

    ET stack[utils::kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&] (const ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    Seq::iterate_seq(n_b2, copy_f);
    assert(offset <= utils::kBaseCaseSize);

    Seq::decrement_recursive(n_b1);
    Seq::decrement_recursive(n_b2);

    ET output[utils::kBaseCaseSize + 1];

    // merge
    size_t nA = offset_2; size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      const auto& k_a = Entry::get_key(stack[i]);
      const auto& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        parlay::move_uninitialized(output[out_off++], stack[i]);
        i++;
      } else if (comp(k_b, k_a)) {
        j++;
      } else {
        i++;
        j++;
      }
    }
    while (i < nA) {
      parlay::move_uninitialized(output[out_off++], stack[i]);
      i++;
    }

    // build tree
    if (out_off < utils::compression_block_size) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

  template <class BinaryOp>
  static node* intersect_bc(ptr b1, ptr b2, const BinaryOp& op) {
    auto n_b1 = b1.node_ptr();
    auto n_b2 = b2.node_ptr();

    ET stack[utils::kBaseCaseSize + 1];
    size_t offset = 0;
    auto copy_f = [&] (const ET& a) {
      parlay::move_uninitialized(stack[offset++], a);
    };
    Seq::iterate_seq(n_b1, copy_f);
    size_t offset_2 = offset;
    Seq::iterate_seq(n_b2, copy_f);
    assert(offset <= utils::kBaseCaseSize);

    Seq::decrement_recursive(n_b1);
    Seq::decrement_recursive(n_b2);

    ET output[utils::kBaseCaseSize + 1];

    // merge
    size_t nA = offset_2; size_t nB = offset;
    size_t i = 0, j = offset_2, out_off = 0;
    while (i < nA && j < nB) {
      const auto& k_a = Entry::get_key(stack[i]);
      const auto& k_b = Entry::get_key(stack[j]);
      if (comp(k_a, k_b)) {
        i++;
      } else if (comp(k_b, k_a)) {
        j++;
      } else {
        parlay::move_uninitialized(output[out_off], stack[i]);
        ET& re = output[out_off];
        Entry::set_val(re, op(Entry::get_val(stack[j]), Entry::get_val(re)));
        out_off++;
        i++;
        j++;
      }
    }

    // build tree
    if (out_off < utils::compression_block_size) {
      return Seq::to_tree_impl((ET*)output, out_off);
    } else {
      // need to refactor
      return Seq::make_compressed(output, out_off);
    }
  }

};

}  // namespace cpam
