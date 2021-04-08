#pragma once

using namespace std;

// *******************************************
//   MAPS
// *******************************************

namespace cpam {

template <class _Entry, class Join_Tree>
class map_ {
public:
  using Entry = _Entry;
  using Seq_Tree = sequence_ops<Join_Tree>;
  using Tree = map_ops<Seq_Tree,Entry>;
  using node = typename Tree::node;
  using E = typename Entry::entry_t;
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using M = map_;
  using GC = typename Tree::GC;
  using Build = build<Entry>;
  using maybe_V = std::optional<V>;
  using maybe_E = std::optional<E>;
  using ptr = typename GC::ptr;

  // initializing, reserving and finishing
  static void init() { GC::init(); }
  static void reserve(size_t n, bool randomize=false) {
    GC::reserve(n/utils::compression_block_size, randomize);};
  static void finish() { GC::finish(); }

  void swap(map_& b) { std::swap(root, b.root); }

  // empty constructor
  map_() : root(NULL) { GC::init(); }

  // copy constructor, increment reference count
  map_(const M& m) {
    root = m.root; GC::increment(root);}

  // move constructor, clear the source, leave reference count as is
  map_(M&& m) {
    root = m.root; m.root = NULL;}

  // copy assignment, clear target, increment reference count,
  M& operator = (const M& m) {
    if (this != &m) {
      clear();
      root = m.root; GC::increment(root);
    }
    return *this;
  }

  // move assignment, clear target, leave reference count as is
  M& operator = (M&& m){
    if (this != &m) { clear(); root = m.root; m.root = NULL;}
    return *this;
  }

  // destruct.
  ~map_() { clear(); }

  // singleton (TODO: finalize?)
  map_(const E& e) { GC::init(); root = Tree::single(e);}

  // construct from an array keeping one of the equal keys
  map_(E* s, E* e) { //, bool seq_inplace = false) {
    M empty = M();
    root = Tree::finalize(multi_insert(empty, parlay::slice<E*,E*>(s,e)).get_root()); }

  // construct from an array keeping one of the equal keys
  template<class Bin_Op>
  map_(E* s, E* e, Bin_Op f) { //, bool seq_inplace = false) {
    M empty = M();
    root = Tree::finalize(multi_insert_combine(empty, parlay::slice<E*,E*>(s,e), f).get_root()); }

  // construct from sequence keeping one of the equal keys
  //template<class Range>
  // map_(Range const &S) {
  // due to overloading cannot use the generic version on range
  map_(parlay::sequence<E> const &S) {
    M empty = M();
    root = Tree::finalize(multi_insert(empty, parlay::make_slice(S)).get_root());
  }

  map_(parlay::slice<E*, E*> R) {
    M empty = M();
    root = Tree::finalize(multi_insert(empty, R).get_root());
  }

  // construct from sequence, combining the values corresponding to equal keys with f
  // f should have type V x V -> V
  template<class Range, class Bin_Op>
  map_(Range S, Bin_Op f) {
    M empty = M();
    root = Tree::finalize(multi_insert_combine(empty, parlay::make_slice(S), f).get_root()); }

  template<class Bin_Op>
  map_(std::initializer_list<E> L, Bin_Op f) {
    M empty = M();
    root = Tree::finalize(multi_insert_combine(empty, L, f).get_root()); }

  // clears contents, decrementing ref counts
  // has a read/destruct race if concurrent with a reader
  void clear() {
    node* t = root;
    if (__sync_bool_compare_and_swap(&(this->root), t, NULL)) {
      if (GC::initialized())
	GC::decrement_recursive(t);
    }
  }

  /* ======================================================================= */

  // some basic functions
  // tree size
  size_t size() const {
    return Tree::size(root); }

  bool is_empty() {return root == NULL;}

  // equality
  bool operator == (const M& m) const {
    return (size() == m.size()) && (size() == map_union(*this,m).size());
  }
  bool operator != (const M& m) const { return !(*this == m); }

  // apply function f on all entries
  template <class F>
  static void foreach_index(const M& m, const F& f, size_t start=0,
			    size_t granularity = utils::node_limit) {
    Tree::foreach_index(ptr(m.root, true), start, f, granularity);
  }

  // apply function f on all entries sequentially
  template <class F>
  static void foreach_seq(M m, const F& f) {
    Tree::foreach_seq(ptr(m.get_root(), true), f);
  }

  // apply function f on all entries sequentially. F returns a boolean
  // indicating whether to proceed further or not.
  template <class F>
  static bool foreach_cond(M m, const F& f, size_t start=0,
      size_t granularity = utils::node_limit) {
    return Tree::foreach_cond(ptr(m.get_root(), true), start, f, granularity);
  }

  // apply function f on all entries sequentially. F returns a boolean
  // indicating whether to proceed further or not.
  template <class F, class C>
  static bool foreach_cond_par(M m, const F& f, const C& cond, size_t start=0,
      size_t granularity = utils::node_limit) {
    return Tree::foreach_cond_par(ptr(m.get_root(), true), start, f, cond, granularity);
  }

  // apply function f to all entries in the tree and flatten them to a sequence
  template <class OT, class F>
  static parlay::sequence<OT> to_seq(M m, const F& f,
			     size_t granularity=utils::node_limit) {
    parlay::sequence<OT> out = parlay::sequence<OT>::uninitialized(m.size());
    auto g = [&] (E& e, size_t i) {
      parlay::assign_uninitialized(out[i],f(e));};
    Tree::foreach_index(m.get_root(), 0, g, granularity);
    return out;
  }

  // flatten all entries to an array
  static void entries(M m, E* out) {
    auto f = [&] (E& e, size_t i) {out[i] = e;};
    Tree::foreach_index(m.get_root(), 0, f);
  }

  // flatten all entries to a sequence
  static parlay::sequence<E> entries(M m, size_t granularity=utils::node_limit) {
    auto f = [] (E e) -> E {return e;};
    return to_seq<E>(std::move(m), f, granularity);
  }

  // flatten all values to an output random access iterator
  template <class outIter>
  static void values(M m, outIter out) {
    auto f = [&] (E& e, size_t i) {out[i] = Entry::get_val(e);};
    Tree::foreach_index(m.get_root(), 0, f);
  }

  // flatten all keys to an output random access iterator
  template <class outIter>
  static void keys(M m, outIter out) {
    auto f = [&] (const E& e, const size_t i) {out[i] = Entry::get_key(e);};
    Tree::foreach_index(m.get_root(), 0, f);
  }

  static parlay::sequence<K> keys(M m, size_t granularity=utils::node_limit) {
    auto f = [] (E e) -> K {return Entry::get_key(e);};
    return to_seq<K>(std::move(m), f, granularity);
  }

  static void keys_to_array(M m, K* out) {
    auto f = [&] (E& e, size_t i) {out[i] = Entry::get_key(e);};
    Tree::foreach_index(m.get_root(), 0, f);
  }


  // insertion, updates, and deletion
  template <class Func>
  static M insert(M m, const E& p, const Func& f) {
    return M(Tree::insert(m.get_root(), p, f)); }

  template <class Func>
  void insert(const E& p, const Func& f) {
    root = Tree::insert(root, p, f); }

  static M insert(M m, const E& p) {
    auto replace = [] (const V& a, const V& b) {return b;};
    return M(Tree::insert(m.get_root(), p, replace)); }

  void insert(const E& p) {
    auto replace = [] (const V& a, const V& b) {return b;};
    root = Tree::insert(root, p, replace); }

// TODO
//  template <class Func>
//  void update(const K& k, const Func& f) {
//    root = Tree::update(root, k, f); }
//
//  template <class Func>
//  static M update(M m, const K& k, const Func& f) {
//    return M(Tree::update(m.get_root(), k, f)); }

  static M remove(M m, const K& k) {
    return M(Tree::deletet(m.get_root(), k)); }

  // Remove entry with key k from the current map.
  void remove(const K& k) {
    root = Tree::deletet(root, k);
  }

  static M join2(M a, M b) {
    return M(Tree::join2(a.get_root(), b.get_root()));
  }


  // filters elements that satisfy the predicate when applied to the elements.
  template<class F>
  static M filter(M m, const F& f, size_t granularity=utils::node_limit) {
    return M(Tree::filter(m.get_root(), f, granularity)); }

  template<class Seq>
  //static M from_sorted(Seq const &S) {
  static M from_sorted(Seq &S) {
    return M(Tree::from_array(S.data(), S.size()));
  }


// TODO
//  // determines if there is any entry in the tree satisfying indicator f
//  template<class F>
//  static bool if_exist(M m, const F& f) {
//	bool* indicator = new bool;
//	*indicator = false;
//    bool ans = Tree::if_exist(m.get_root(), f, indicator);
//	delete indicator;
//	return ans;
//  }

// TODO
//  // insert multiple keys from an array
//  template<class Seq>
//  static M multi_delete(M m, Seq const &SS) {
//    //bool seq_inplace = false, bool inplace = false) {
//
//    using EE = typename Seq::value_type;
//
//    struct e_type {
//      using entry_t = EE;
//      using key_t = EE;
//      using val_t = bool;
//      static bool comp(const key_t& a, const key_t& b) {return a<b;}
//      static inline key_t get_key(const entry_t& e) {return e;}
//      static inline val_t get_val(const entry_t& e) {return false;}
//    };
//
//    using BD = build<e_type>;
//    parlay::sequence<K> A = BD::sort_remove_duplicates(SS);
//
//    auto x = M(Tree::multi_delete_sorted(m.get_root(), A.data(), A.size()));
//    return x;
//  }

  // insert multiple keys from a sequence
  template<class Seq>
  static M multi_insert(M m, Seq const &SS) {
    auto replace = [] (const V& a, const V& b) {return b;};
    parlay::sequence<E> A = Build::sort_remove_duplicates(SS);
    M A_m = Tree::multi_insert_sorted(nullptr, A.data(), A.size(), replace);
    auto x = M(Tree::uniont(m.get_root(), A_m.get_root(), replace));
    //auto x = M(Tree::multi_insert_sorted(m.get_root(), A.data(), A.size(), replace));
    return x;
  }

  // insert multiple keys from a sequence, destroying the sequence
  // ?? pass by rvalue reference instead ??
  // template<class Seq>
  // static M multi_insert_xx(M m, Seq &SS, bool sequential = false) {
  //   auto replace = [] (const V& a, const V& b) {return b;};
  //   parlay::sequence<E> A = Build::sort_remove_duplicates(std::move(SS), sequential);
  //   auto x = M(Tree::multi_insert_sorted(m.get_root(), A.data(), A.size(), replace));
  //   return x;
  // }

// TODO
//  // update multiple entries from a sorted array
//  //template<class Seq>
//  template<class Seq, class Bin_Op>
//  static M multi_update_sorted(M m, Seq const &SS, Bin_Op f) {
//    return M(Tree::multi_update_sorted(m.get_root(), SS.data(), SS.size(), f));
//  }

// TODO
//  // update multiple entries from an array
//  //template<class Seq>
//  template<class Seq, class Bin_Op>
//  static M multi_update(M m, Seq const &SS, Bin_Op f) {
//    using EE = typename Seq::value_type;
//    using VV = typename EE::second_type;
//
//    struct e_type {
//      using entry_t = EE;
//      using key_t = K;
//      using val_t = VV;
//      static bool comp(const key_t& a, const key_t& b) {return a<b;}
//      static inline key_t get_key(const entry_t& e) {return e.first;}
//      static inline val_t get_val(const entry_t& e) {return e.second;}
//      static inline void set_val(entry_t& e, const val_t& v) {e.second = v;}
//    };
//
//    using BD = build<e_type>;
//    parlay::sequence<pair<K, VV>> B = BD::sort_remove_duplicates(SS);
//    return M(Tree::multi_update_sorted(m.get_root(), B.data(), B.size(), f));
//  }

// TODO
//  template<class Seq>
//  static V* multi_find(M m, Seq const &SS) {
//    using K = typename Seq::value_type;
//    auto less = [&] (K& a, K &b) {return Entry::comp(a,b);};
//    parlay::sequence<K> B = parlay::internal::sample_sort(SS, less);
//    V* ret = new V[B.size()];
//    Tree::multi_find_sorted(m.get_root(), B.data(), B.size(), ret, 0);
//    return ret;
//  }

  // insert multiple keys from an array
  template<class Seq>
  static M multi_insert_sorted(M m, Seq const &SS) {
    auto replace = [] (const V& a, const V& b) {return b;};
    return M(Tree::multi_insert_sorted(m.get_root(), SS.data(),
				       SS.size(), replace));
  }

  // insert multiple keys from an array, combine duplicates with f
  // here f must have type: V x V -> V
  // if key in map, then also combined with f
  template<class Seq, class Bin_Op>
  static M multi_insert_combine(M m, Seq const &S, Bin_Op f) {
    auto A = Build::sort_combine_duplicates(S, f);
    return M(Tree::multi_insert_sorted(m.get_root(),
				       A.data(), A.size(), f));
  }

  // What is this function _xx? in-place?
  template<class Seq, class Bin_Op>
  static M multi_insert_combine_xx(M m, Seq &S, Bin_Op f) {
    auto A = Build::sort_combine_duplicates_inplace(S, f);
    return M(Tree::multi_insert_sorted(m.get_root(),
				       A.begin(), A.size(), f));
  }

  template <class Func>
  void iterate_seq(const Func& f) {
    Tree::iterate_seq(root, f);
  }

  // insert multiple keys from an array, reduce duplicates with g
  // here g must have type: sequence<Val> -> V
  // if key in map, then replaced
  template<class Seq, class Reduce>
  static M multi_insert_reduce(M m, Seq const &S, Reduce g) {
    auto replace = [] (const V& a, const V& b) {return b;};
    parlay::sequence<E> A = Build::sort_reduce_duplicates(S, g);
    auto x =  M(Tree::multi_insert_sorted(m.get_root(),
					  A.data(), A.size(), replace));
    //auto print_fn = [&] (const E& e) {
    //  std::cout << e.first << " " << e.second.size() << std::endl;
    //};
    //Tree::iterate_seq(x.root, print_fn);
    return x;
  }

  template<class Val, class Reduce>
  static M multi_insert_reduce(M m, pair<K,Val>* A, size_t n, Reduce g) {
    auto replace = [] (const V& a, const V& b) {return b;};
    auto B = Build::sort_reduce_duplicates(A, n, g);
    auto x =  M(Tree::multi_insert_sorted(m.get_root(),
					  B.first, B.second, replace));
    return x;
  }

  template<class Val, class Reduce, class Bin_Op>
  static M multi_insert_reduce(M m, pair<K,Val>* A, size_t n,
			       Reduce g,
			       Bin_Op& f = [] (const V& a, const V& b) {
				 return b;}) {
    auto B = Build::sort_reduce_duplicates(A, n, g);
    auto x =  M(Tree::multi_insert_sorted(m.get_root(),
    				       B.first, B.second, f));
    return x;
  }

  // basic search routines
  maybe_V find(const K& key) const {
    return maybe_entry_to_val(Tree::find(ptr(root, true), key));}

  // returns default value if not found
  V find(const K& key, V defaultv) const {
    auto et_opt = Tree::find(ptr(root, true), key);
    return et_opt.has_value() ? Entry::get_val(*et_opt) : defaultv;
  }

  bool contains(const K& key) const {
    return Tree::find(ptr(root, true), key).has_value();}

// TODO
//  maybe_E next(const K& key) const {
//    return node_to_entry(Tree::next(root, key));}
//
//  maybe_E previous(const K& key) const {
//      return node_to_entry(Tree::previous(root, key));}
//
//  // rank and select
//  size_t rank(const K& key) { return Tree::rank(root, key);}

  maybe_E select(const size_t rank) const {
    auto ret = Tree::select(ptr(root, true), rank);
    auto unpack = *ret;
    // std::cout << Entry::get_key(unpack) << std::endl;
    // std::cout << Entry::get_val(unpack).size() << std::endl;
    return ret;
  }

// TODO
//  maybe_E last() const { return select(size()-1); }

// TODO
//  M take(const size_t k) {
//    return M(Tree::take(root, k));}

  // extra indicates if there is an extra pointer to b2 not included in the
  // ref-cnt.
  template<class F>
  static M map_union(M a, M b, const F& op, bool extra = false) {
    return M(Tree::uniont(a.get_root(), ptr(b.get_root(), extra), op));
  }

  static M map_union(M a, M b, bool extra = false) {
    auto get_right = [] (V a, V b) {return b;};
    auto x = M(Tree::uniont(a.get_root(), ptr(b.get_root(), extra), get_right));
    return x;
  }

  template<class M1, class M2, class F>
  static M map_intersect(M1 a, M2 b, const F& op) {
    using T1 = typename M1::Tree;
    using T2 = typename M2::Tree;
    return M(Tree::template intersect<T1,T2>(a.get_root(),
					     b.get_root(), op));
  }

  static M map_intersect(M a, M b) {
    auto get_right = [] (V a, V b) {return b;};
    return M(Tree::template intersect<Tree,Tree>(a.get_root(),
						 b.get_root(), get_right));
  }

  static M map_difference(M a, M b) {
    return M(Tree::difference(a.get_root(), b.get_root()));
  }

// TODO
//  static M range(M& a, K kl, K kr) {
//    return M(Tree::range(a.root, kl, kr));
//  }
//
//  static M range_number(M& a, K kl, size_t r) {
//    return M(Tree::range_num(a.root, kl, r));
//  }
//
//  template<class Map, class Reduce>
//  static typename Reduce::T range_number_mr(M& a, K kl, size_t r, const Map& mp, const Reduce& rdc) {
//    auto x = Tree::range_num_mr(a.root, kl, r, mp, rdc);
//    return x.first;
//  }
//
//  static M upTo(M& a, K kr) {
//    return M(Tree::left(a.root, kr));
//  }
//
//  static M downTo(M& a, K kr) {
//    return M(Tree::right(a.root, kr));
//  }
//
//  template<class Ma, class F>
//  static M map(Ma a, const F& f) {
//    GC::init();
//    return M(Tree::template map<typename Ma::Tree>(a.root, f));
//  }
//
//  template<class Ma, class F>
//  static M map_set(Ma a, const F& f) {
//    GC::init();
//    return M(Tree::template map_set<typename Ma::Tree>(a.root, f));
//  }
//
//  template<class F>
//  static void map_index(M& m, const F& f,
//			size_t granularity=utils::node_limit,
//			size_t start=0) {
//    Tree::foreach_index(m.root, start, f, granularity, true);
//  }
//
//  template<class R, class F>
//  static typename R::T semi_map_reduce(const M& m, const F& f, const R& r,
//				       size_t grain=utils::node_limit) {
//    return Tree::template semi_map_reduce<R>(m.root, f, r, grain);}
//
//  template<class R, class F>
//  static typename R::T map_reduce(const M& m, const F& f, const R& r,
//				   size_t grain=utils::node_limit) {
//    GC::init();
//    return Tree::template map_reduce<R>(m.root, f, r, grain);
//  }
//
//  template<class F>
//  static void map_void(M& m, const F& f,
//		       size_t granularity=utils::node_limit) {
//    struct do_nothing {
//      using T = bool;
//      static T identity() {return false;}
//      static T add(T a, T b) {return false;}
//    };
//    auto g = [&] (E v) {f(v); return false;};
//    Tree::map_reduce(m.root, g, do_nothing(), granularity);
//  }
//
//  template<class Ma, class F>
//  static  M map_filter(const Ma& a, const F& f, size_t granularity=utils::node_limit) {
//    return M(Tree::template map_filter<typename Ma::Tree>(a.root, f, granularity));
//  }

  // grabs root by "moving" it.  Important for reuse
  node* get_root() {node* t = root; root = NULL; return t;};

  node* root;

  // construct from a node (perhaps should be private)
  map_(node* n) : root(n) { GC::init(); }

  // TODO: remove?
  //maybe_V node_to_val(node* a) const {
  //  if (a != NULL) return maybe_V(Entry::get_val(Tree::get_entry(a)));
  //  else return maybe_V();
  //}

  maybe_V maybe_entry_to_val(maybe_E e) const {
    if (e.has_value()) return maybe_V(Entry::get_val(*e));
    else return maybe_V();
  }

  // TODO: remove?
  //maybe_E node_to_entry(node* a) const {
  //  if (a != NULL) return maybe_E(Tree::get_entry(a));
  //  else return maybe_E();
  //}


  // Some useful debugging utilities
  bool check_balance() const {
	  return Seq_Tree::check_balance(root);
  }

  template <class F>
  size_t size_in_bytes(const F& f) const {
    return Tree::size_in_bytes(root, f);
  }

  size_t ref_cnt() const { return Tree::ref_cnt(root); }
  size_t root_is_compressed() const { return Tree::is_compressed(root); }
  void print_root_info() const {
    Tree::print_node_info(root, "root");
  }
  void print_tree() const {
    //Tree::print_inorder(root);
    Tree::print_rec(root);
  }
  bool check_structure() {
    return Tree::check_structure(root);
  }
  void print_stats() {
    GC::print_stats();
  }


};

// creates a key-value pair for the entry
template <class entry>
struct map_full_entry : entry {
  using val_t = typename entry::val_t;
  using key_t = typename entry::key_t;
  using entry_t = std::tuple<key_t, val_t>;
  static inline key_t get_key(const entry_t& e) { return std::get<0>(e); }
  static inline val_t get_val(const entry_t& e) { return std::get<1>(e); }
  static inline void set_val(entry_t& e, const val_t& v) { std::get<1>(e) = v; }
  static inline entry_t to_entry(const key_t& k, const val_t& v) {
    return std::make_tuple(k, v);
  };
};

template <class _Entry, class Balance = weight_balanced_tree, class Encoder=default_entry_encoder>
using pam_map = map_<
    map_full_entry<_Entry>, /* Entry type */
    typename Balance::template balance<basic_node<
        typename Balance::data,
        typename map_full_entry<_Entry>::entry_t,
        typename Encoder::template encoder<map_full_entry<_Entry>>>>>;

template <class _Entry, class Balance=weight_balanced_tree>
using diff_encoded_map = pam_map<_Entry, Balance, diffencoded_entry_encoder>;

// entry is just the key (no value), for use in sets
template <class entry>
struct set_full_entry : entry {
  using key_t = typename entry::key_t;
  using val_t = bool;  // not used
  using entry_t = key_t;
  static inline key_t get_key(const entry_t& e) { return e; }
  static inline val_t get_val(const entry_t& e) { return 0; }
  static inline void set_val(entry_t& e, const val_t& v) {}
};

template <class _Entry, class Balance = weight_balanced_tree, class Encoder=default_entry_encoder>
using pam_set =
    map_<set_full_entry<_Entry>,
         typename Balance::template balance<basic_node<
             typename Balance::data,
             typename _Entry::key_t,
             typename Encoder::template encoder<_Entry>>>>;

// entry is just a value (no key), for use in sequences
template <class data>
struct seq_full_entry {
  using key_t = data;
  using val_t = bool;  // not used
  using entry_t = data;
  static bool comp(const key_t& a, const key_t& b) { return true; }
  static inline key_t get_key(const entry_t& e) { return e; }
  static inline val_t get_val(const entry_t& e) { return 0; }
  static inline void set_val(entry_t& e, const val_t& v) {}
};

// data can be any type
template <typename data, class Balance = weight_balanced_tree, class Encoder=default_entry_encoder>
using pam_seq =
    map_<seq_full_entry<data>,
          typename Balance::template balance<basic_node<
                typename Balance::data,
                data,
                typename Encoder::template encoder<seq_full_entry<data>>>>>;



//template <class _Entry, class Balance=weight_balanced_tree>
//using pam_map =
//  map_<map_full_entry<_Entry>,
//       typename Balance::template
//       balance<basic_node<typename Balance::data,
//			  typename map_full_entry<_Entry>::entry_t>>>;
//
//// entry is just the key (no value), for use in sets
//template <class entry>
//struct set_full_entry : entry {
//  using key_t = typename entry::key_t;
//  using val_t = bool;  // not used
//  using entry_t = key_t;
//  static inline key_t get_key(const entry_t& e) {return e;}
//  static inline val_t get_val(const entry_t& e) {return 0;}
//  static inline void set_val(entry_t& e, const val_t& v) {}
//};
//
//template <class _Entry, class Balance=weight_balanced_tree>
//using pam_set =
//  map_<set_full_entry<_Entry>,
//       typename Balance::template
//       balance<basic_node<typename Balance::data,
//			  typename _Entry::key_t>>>;
//
//// entry is just a value (no key), for use in sequences
//template <class data>
//struct seq_full_entry {
//  using key_t = data;
//  using val_t = bool;  // not used
//  using entry_t = data;
//  static bool comp(const key_t& a, const key_t& b) {return true;}
//  static inline key_t get_key(const entry_t& e) {return e;}
//  static inline val_t get_val(const entry_t& e) {return 0;}
//  static inline void set_val(entry_t& e, const val_t& v) {}
//};
//
//// data can be any type
//template <typename data, class Balance=weight_balanced_tree>
//using pam_seq =
//  map_<seq_full_entry<data>,
//       typename Balance::template
//       balance<basic_node<typename Balance::data, data>>>;


}  // namespace cpam
