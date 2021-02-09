#pragma once

using namespace std;

// *******************************************
//   AUG MAPS
// *******************************************

template <class _Entry, class Join_Tree>
struct aug_map_ : private map_<_Entry, Join_Tree> {
public:
  using Map = map_<_Entry, Join_Tree>;
  using Entry = typename Map::Entry;
  using Tree = augmented_ops<typename Map::Tree>;
  using node = typename Tree::node;
  using E = typename Entry::entry_t;
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using A = typename Entry::aug_t;
  using M = aug_map_;
  using GC = typename Map::GC;
  using maybe_V = std::optional<V>;
  using maybe_E = std::optional<E>;

  template<class F>
  static M aug_filter(M m, const F& f) {
    return M(Tree::aug_filter(m.get_root(), f)); }

  // extract the augmented values
  A aug_val() { return Tree::aug_val(Map::root); }

  A aug_left (const K& key) {
    typename Tree::aug_sum_t a;
    Tree::aug_sum_left(Map::root, key, a);
    return a.result;}

  A aug_right(const K& key) {
    typename Tree::aug_sum_t a;
    Tree::aug_sum_right(Map::root, key, a);
    return a.result;}

  A aug_range(const K& key_left, const K& key_right) {
    typename Tree::aug_sum_t a;
    Tree::aug_sum_range(Map::root, key_left, key_right, a);
    return a.result;}

  // just side effecting
  template <class restricted_sum>
  void range_sum(const K& key_left,
		 const K& key_right,
		 restricted_sum& rs) {
    Tree::aug_sum_range(Map::root, key_left, key_right, rs);
  }

  template <class Func>
  maybe_E aug_select(Func f) {
    return Map::node_to_entry(Tree::aug_select(Map::root, f));};

  maybe_E aug_eq(A aug_val) {
    return Map::node_to_entry(Tree::aug_eq(Map::root, aug_val));
  }

  static M insert_lazy(M m, const E& p) {
    auto replace = [] (const V& a, const V& b) {return b;};
    return M(Tree::insert_lazy(m.get_root(), p, replace)); }

  // for coercing a map to an aug_map, should be a better way
  static M to_aug(Map&& m) {
    aug_map_ x;
    x.root = m.root; m.root = NULL; return x;}
  aug_map_() : Map() { }
  // install Map's constructors
  using Map::Map;

  template <class Func>
  static M insert(M m, const E& p, const Func& f) {
    return to_aug(Map::insert(std::move(m),p,f));}
  static M insert(M m, const E& p) {//cout << "ins a " << endl;
    return to_aug(Map::insert(std::move(m),p));}
  static M remove(M m, const K& k) {return to_aug(Map::remove(std::move(m), k));}
  template<class F>
  static M filter(M m, const F& f) {return to_aug(Map::filter(std::move(m), f));}
  static M multi_insert(M m, parlay::sequence<E> const &SS) {
    return to_aug(Map::multi_insert(std::move(m), SS));}

  // ??
  //static M multi_insert(M m, parlay::sequence<E> &&SS, bool sequential) { // ?? should it be &
  //  return Map::multi_insert_xx(std::move(m), std::move(SS), sequential);}

  static M multi_delete(M m, parlay::sequence<K> SS, bool seq_inplace = false) {  // ?? should it be &
    return to_aug(Map::multi_delete(std::move(m), SS, seq_inplace));}

  template<class Bin_Op>
  static M multi_insert_combine(M m, parlay::sequence<E> S, Bin_Op f,  // ?? should it be &
				bool seq_inplace = false) {
    return to_aug(Map::multi_insert_combine(std::move(m), S, f, seq_inplace));}
  template<class Val, class Reduce>
  static M multi_insert_reduce(M m, parlay::sequence<pair<K,Val>> S, Reduce g) {  // ?? should it be &
    return to_aug(Map::multi_insert_reduce(std::move(m), S, g)); }
  template<class M1, class M2, class F>
  static M map_intersect(M1 a, M2 b, const F& op) {
    return to_aug(Map::map_intersect(std::move(a), std::move(b), op));}
  static M map_intersect(M a, M b) {return to_aug(Map::map_intersect(std::move(a), std::move(b)));}
  template<class F>
  static M map_union(M a, M b, const F& op) {return to_aug(Map::map_union(std::move(a), std::move(b), op));}
  static M map_union(M a, M b) {return to_aug(Map::map_union(std::move(a), std::move(b)));}
  static M map_difference(M a, M b) {return to_aug(Map::map_difference(std::move(a), std::move(b)));}
  static M join2(M a, M b) {return to_aug(Map::join2(std::move(a), std::move(b)));}
  static M range(M& a, K kl, K kr) {return to_aug(Map::range(a,kl,kr));}
  static M upTo(M& a, K kr) {return to_aug(Map::upTo(a,kr));}
  template<class Ma, class F>
  static M map(Ma a, const F f) {return to_aug(Map::map(a, f));}
  template<class F>
  static void map_void(M& a, const F& f,
		       size_t granularity=utils::node_limit) {
    return Map::map_void(a, f, granularity);
  }
  static void entries(M m, E* out) { Map::entries(std::move(m),out);}
  static parlay::sequence<E> entries(M m, size_t granularity=utils::node_limit) { return Map::entries(std::move(m));}
  template <class outItter>
  static void keys(M m, outItter out) {Map::keys(std::move(m),out);}
  static void keys_to_array(M m, K* out) {Map::keys_to_array(std::move(m),out);}
  static parlay::sequence<K> keys(M m, size_t granularity=utils::node_limit) {
	  return Map::keys(m, granularity);
  }
  bool operator == (const M& m) { return Map::operator==(m);}
  template<class R, class F>
  static typename R::T map_reduce(const M& m, const F& f, const R& r,
				  size_t grain=utils::node_limit) {
    return Map::template map_reduce<R>(m, f, r, grain);}
  template<class F>
  static void map_index(M m, const F& f, size_t granularity = utils::node_limit,
			size_t start=0) {
    Map::map_index(m, f, granularity, start); }
  template<class Ma, class F>
  static M map_filter(Ma a, const F& f) {return to_aug(Map::map_filter(a,f));}
  template<class F>
  static bool if_exist(M m, const F& f) {return Map::if_exist(m,f);}
  template<class Ma, class F>
  static M map_set(Ma a, const F& f) {return to_aug(Map::map_set(a, f));}
  template <class F>
  static void foreach_index(M m, F f, size_t start=0,
			    size_t granularity = utils::node_limit) {
    Map::foreach_index(m, f, start, granularity); }
  template <class F>
  static void foreach_seq(M m, F f) {
    Map::foreach_seq(m, f); }
public:
  using Map::from_sorted;
  using Map::size;
  using Map::is_empty;
  using Map::init;
  using Map::reserve;
  using Map::finish;
  using Map::clear;
  using Map::find;
  using Map::contains;
  using Map::next;
  using Map::previous;
  using Map::rank;
  using Map::select;
  using Map::root;
  using Map::get_root;
  using Map::insert;
  using Map::check_balance;
  using Map::size_in_bytes;
};

// creates a key-value pair for the entry, and redefines from_entry
template <class entry>
struct aug_map_full_entry : entry {
  using val_t = typename entry::val_t;
  using key_t = typename entry::key_t;
  using aug_t = typename entry::aug_t;
  using entry_t = std::tuple<key_t,val_t>;
  static inline key_t get_key(const entry_t& e) {return std::get<0>(e);}
  static inline val_t get_val(const entry_t& e) {return std::get<1>(e);}
  static inline void set_val(entry_t& e, const val_t& v) {std::get<1>(e) = v;}
  static inline aug_t from_entry(const entry_t& e) {
    return entry::from_entry(std::get<0>(e), std::get<1>(e));}
};

template <class _Entry, class Balance=weight_balanced_tree>
using aug_map =
  aug_map_<aug_map_full_entry<_Entry>,
  typename Balance::template
  balance<aug_node<typename Balance::data,
		   aug_map_full_entry<_Entry>>>>;

// creates a key-value pair for the entry, and redefines from_entry
template <class entry>
struct aug_set_full_entry : entry {
  using val_t = bool; // not used
  using key_t = typename entry::key_t;
  using aug_t = typename entry::aug_t;
  using entry_t = key_t;
  static inline key_t get_key(const entry_t& e) {return e;}
  static inline val_t get_val(const entry_t& e) {return 0;}
  static inline void set_val(entry_t& e, const val_t& v) {}
};

// _Entry needs:
//    key_t, aug_t,
//    comp(key_t, key_t) -> bool,
//    from_entry(key_t) -> aug_t,
//    get_empty() -> aug_tm,
//    combine(aug_t, aug_t) -> aug_t
template <class _Entry, class Balance=weight_balanced_tree>
using aug_set =
  aug_map_<aug_set_full_entry<_Entry>,
  typename Balance::template
  balance<aug_node<typename Balance::data,
		   aug_set_full_entry<_Entry>>>>;
