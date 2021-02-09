#pragma once
#include "basic_node.h"

// *******************************************
//   AUGMENTED NODE
// *******************************************

namespace cpam {

// Creates an augmented node from a basic node.
// The new augmented entry is a pair of the original entry and the aumented
// value.   The Entry struct must have the static inteface:
//   entry_t;
//   aug_t;
//   get_empty() -> aug_t;
//   from_entry(entry_t) -> aug_t;
//   combine(aut_t, aug_t) -> aug_t;
template<class balance, class Entry, class AugEntryEncoder>
struct aug_node : private basic_node<balance,
    std::pair<typename Entry::entry_t, typename Entry::aug_t>,
    null_encoder> {
  using AT = typename Entry::aug_t;
  using ET = typename Entry::entry_t;
  using basic = basic_node<balance, std::pair<ET,AT>, null_encoder>;
  using node = typename basic::node;
  using regular_node = typename basic::regular_node;
  using compressed_node = typename basic::compressed_node;
  using aug = aug_node<balance, Entry, AugEntryEncoder>;
  using allocator = typename basic::allocator;


  using basic::increment_count;
  using basic::empty;
  using basic::ref_cnt;
  using basic::is_regular;
  using basic::is_compressed;
  using basic::cast_to_regular;
  using basic::cast_to_compressed;
  using basic::size;
  using basic::will_be_compressed;
  using basic::check_compressed_node;


  static ET& get_entry(node *a) {return basic::cast_to_regular(a)->entry.first;}
  static ET* get_entry_p(node *a) {return &basic::cast_to_regular(a)->entry.first;}
  static void set_entry(node *a, ET e) {basic::cast_to_regular(a)->entry.first = e;}

  struct aug_compressed_node {
    node_size_t r;  // reference count, top-bit is always "0"
    node_size_t s;  // number of entries used (size)
    node_size_t size_in_bytes;
    AT aug_val;
  };

  // TODO: include case for compressed
  static AT aug_val(node* a) {
    if (a == NULL) return Entry::get_empty();
    if (basic::is_regular(a)) return (basic::cast_to_regular(a)->entry).second;
    auto c = (aug_compressed_node*)a;
    return c->aug_val;
  }

  static void update(node* ov) {
    regular_node* a = basic::cast_to_regular(ov);
    basic::update(a);
    AT av = Entry::from_entry(get_entry(a));
    // TODO: handle extracting aug_val from children using call to aug_val
    // (above)
    if (a->lc) av = Entry::combine(aug_val(a->lc), av);
    if (a->rc) av = Entry::combine(av, aug_val(a->rc));
    (a->entry).second = av;
  }

  // Handles both regular and compressed nodes.
  static void free_node(node* va) {
    if (basic::is_regular(va)) {
      auto a = basic::cast_to_regular(va);
      std::get<0>((a->entry)).~ET();
      allocator::free(a);
    } else {
      auto c = basic::cast_to_compressed(va);
      uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t) + sizeof(AT));
      AugEntryEncoder::destroy(data_start, c->s);
      auto array_size = c->size_in_bytes;
      utils::free_array<uint8_t>((uint8_t*)va, array_size);
    }
  }


  // precondition: a is non-null, and caller has a reference count on a.
  static bool decrement_count(node* a) {
    if ((utils::fetch_and_add(&basic::generic_node(a)->r, -1) & basic::kLowBitMask) == 1) {
      free_node(a);
      return true;
    }
    return false;
  }

  // atomically decrement ref count and deletes node if zero
  static bool decrement(node* t) {
    if (t) {
      return decrement_count(t);
    }
    return false;
  }

  // atomically decrement ref count and if zero:
  //   delete node and recursively decrement the two children
  static void decrement_recursive(node* t) {
    if (!t) return;
    if (basic::is_regular(t)) {
      node* lsub = basic::cast_to_regular(t)->lc;
      node* rsub = basic::cast_to_regular(t)->rc;
      if (decrement(t)) {
        utils::fork_no_result(basic::size(lsub) >= utils::node_limit,
                              [&]() { decrement_recursive(lsub); },
                              [&]() { decrement_recursive(rsub); });
      }
    } else {
      decrement(t);
    }
  }


// TODO
//  // updates augmented value using f, instead of recomputing
//  template<class F>
//  static void lazy_update(node* a, F f) {
//    basic::update(a);
//    (a->entry).second = f((a->entry).second);
//  }

  static regular_node* make_regular_node(ET e) {
    std::pair<ET,AT> ea;
    ea.first = e;
    return basic::make_regular_node(ea);
  }

  static regular_node* single(ET e) {
    AT av = Entry::from_entry(e);
    return basic::single(std::make_pair(e,av));
  }


  template<typename F>
  static void iterate_seq(node* a, const F& f) {
    if (!a) return;
    if (basic::is_regular(a)) {
      auto r = basic::cast_to_regular(a);
      iterate_seq(r->lc, f);
      f(get_entry(r));
      iterate_seq(r->rc, f);
    } else {
      auto c = basic::cast_to_compressed(a);
      uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t) + sizeof(AT));
      AugEntryEncoder::decode(data_start, c->s, f);
    }
  }

  template<typename F>
  static bool iterate_cond(node* a, const F& f) {
    if (!a) return true;
    if (is_regular(a)) {
      auto r = cast_to_regular(a);
      bool ret = iterate_cond(r->lc, f);
      if (!ret) return ret;
      ret = f(get_entry(r));
      if (!ret) return ret;
      ret = iterate_cond(r->rc, f);
      return ret;
    } else {
      auto c = cast_to_compressed(a);
      uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t) + sizeof(AT));
      return AugEntryEncoder::decode_cond(data_start, c->s, f);
    }
  }


  static node* finalize(node* root) {
    auto sz = basic::size(root);
    assert(sz > 0);
    if (sz < B) {
      auto ret = make_compressed_node(root);
      decrement_recursive(root);
      return ret;
    }
    return root;
  }

  // Used by GC to copy a compressed node. TODO: update to work correctly with
  // diff-encoding.
  static node* make_compressed_node(node* b) {
    return basic_node_helpers::make_compressed_node<aug>(b);
  }

  // takes a pointer to an array of ETs, and a length of the number of ETs to
  // construct, and returns a compressed node.
  static compressed_node* make_single_compressed_node(ET* e, size_t s) {
//    assert(s >= B);
    assert(s <= 2*B);

    size_t encoded_size = AugEntryEncoder::encoded_size(e, s);
    size_t node_size = 3*sizeof(node_size_t) + sizeof(AT) + encoded_size;
    aug_compressed_node* c_node = (aug_compressed_node*)utils::new_array_no_init<uint8_t>(node_size);

    c_node->r = 1;
    c_node->s = s;
    c_node->size_in_bytes = node_size;

    uint8_t* encoded_data = (((uint8_t*)c_node) + 3*sizeof(node_size_t) + sizeof(AT));
    c_node->aug_val = AugEntryEncoder::encode(e, s, encoded_data);

    // basic::check_compressed_node(c_node);
    return (compressed_node*)c_node;
  }

  static node* make_compressed(node* l, node* r, regular_node* e) {
    return basic_node_helpers::make_compressed<aug>(l, r, e);
  }

  static node* make_compressed(ET* stack, size_t tot) {
    return basic_node_helpers::make_compressed<aug>(stack, tot);
  }

  static node* make_compressed(node* l, node* r) {
    return make_compressed(l, r, nullptr);
  }

  static ET* compressed_node_elms(compressed_node* c, ET* tmp_arr) {
    uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t) + sizeof(AT));
    size_t i = 0;
    auto f = [&] (const ET& et) {
      parlay::assign_uninitialized(tmp_arr[i++], et);
    };
    AugEntryEncoder::decode(data_start, c->s, f);
    return tmp_arr;
  }

  // F applied to entries.
  template <class F>
  static size_t size_in_bytes(node* a, const F& f) {
    size_t total = 0;
    if (!a) return total;
    if (basic::is_compressed(a)) {
      auto c = basic::cast_to_compressed(a);
      total += c->size_in_bytes;
      auto fn = [&] (const auto& et) {
        total += f(et);
      };
      iterate_seq(c, fn);
    } else {
      auto r = basic::cast_to_regular(a);
      total += size_in_bytes(r->lc, f);
      total += size_in_bytes(r->rc, f);
      total += sizeof(regular_node);
      total += f(get_entry(r));
    }
    return total;
  }

// Unused now.
//  static node* make_node(ET e) {
//    std::pair<ET,AT> ea;
//    ea.first = e;
//    return basic::make_node(ea);
//  }
//
//  static node* single(ET e) {
//    AT av = Entry::from_entry(e);
//    return basic::single(std::make_pair(e,av));
//  }
};

}  // namespace cpam
