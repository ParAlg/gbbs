#pragma once
#include "parlay/alloc.h"

#include "utils.h"
#include "basic_node_helpers.h"
#include "byte_encode.h"
#include "compression.h"

namespace cpam {

using node_size_t = unsigned int;
//using node_size_t = size_t;
static constexpr size_t B = utils::compression_block_size;

// *******************************************
//   BASIC NODE
// *******************************************

template<class balance_data, class _Entry, class EntryEncoder>
struct basic_node {
  using ET = _Entry;
  using node = void;
  using basic = basic_node<balance_data, _Entry, EntryEncoder>;

  struct regular_node : balance_data {
    node_size_t r; // reference count, top-bit is always "1"
    node_size_t s;
    node* lc;
    node* rc;
    ET entry;
  };
  using allocator = parlay::type_allocator<regular_node>;

  struct compressed_node {
    node_size_t r;              // reference count, top-bit is always "0"
    node_size_t s;              // number of entries used (size)
    node_size_t size_in_bytes;  // space allocated in bytes.
  };

  static bool is_complex(node* a) {
    assert(is_regular(a));
    return size(a) > 2*B; }
  static bool is_simplex(node* a) {
    return size(a) <= 2*B; }

  static bool is_regular(node* a) {
    return !a || ((regular_node*)a)->r & kTopBit; }
  static bool is_compressed(node* a) { return !is_regular(a); }

  // Useful if running in debug mode for catching bugs.
  static regular_node* cast_to_regular(node* a) {
    assert(is_regular(a));
    return (regular_node*)a;
  }
  static compressed_node* cast_to_compressed(node* a) {
    assert(is_compressed(a));
    return (compressed_node*)a;
  }
  static regular_node* generic_node(node* a) {
    return (regular_node*)a;
  }
  static node_size_t size(node* a) { return (a == NULL) ? 0 : generic_node(a)->s; }
  static node_size_t ref_cnt(node* a) { return generic_node(a)->r & kLowBitMask; }



  // Used by balance_utils
  static void update(node* na) {
    auto a = cast_to_regular(na);
    a->s = size(a->lc) + size(a->rc) + 1; }

  static regular_node* make_regular_node(const ET& e) {
    regular_node* o = allocator::alloc();
    o->r = 1;
    o->r |= kTopBit;
    parlay::assign_uninitialized(o->entry, e);
    return o;
  }

  static regular_node* single(const ET& e) {
    regular_node* r = make_regular_node(e);
    r->lc = r->rc = NULL;
    r->s = 1;
    return r;
  }


  static node* empty() {return NULL;}

  // TODO: seems that get_entry should handle the case where a is a leaf!
  inline static ET& get_entry(node* a) { return cast_to_regular(a)->entry; }
  inline static ET* get_entry_p(node* a) { return &(cast_to_regular(a)->entry); }
  static void set_entry(node* a, ET e) { cast_to_regular(a)->entry = e; }

  static constexpr node_size_t kTopBit = ((node_size_t)1) << (sizeof(node_size_t)*8-1);
  static constexpr node_size_t kLowBitMask = kTopBit - ((node_size_t)1);

  // TODO
  // static node* left(node a) {return a.lc;}
  // static node* right(node* a) {return a.rc;}


  /* ========================== Compression ============================= */

  template<typename F>
  static void iterate_seq(node* a, const F& f) {
    if (!a) return;
    if (is_regular(a)) {
      auto r = cast_to_regular(a);
      iterate_seq(r->lc, f);
      f(get_entry(r));
      iterate_seq(r->rc, f);
    } else {
      auto c = cast_to_compressed(a);
      uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t));
      EntryEncoder::decode(data_start, c->s, f);
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
      uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t));
      return EntryEncoder::decode_cond(data_start, c->s, f);
    }
  }

  // Used by GC to copy a compressed node. TODO: update to work correctly with
  // diff-encoding.
  static node* make_compressed_node(node* b) {
    return basic_node_helpers::make_compressed_node<basic>(b);
  }

  static bool check_compressed_node(node* a) {
    auto c = cast_to_compressed(a);
    assert(c);
    assert(!is_regular(c));
    //assert(B <= size(c) && size(c) <= 2*B);
    assert(size(c) <= 2*B);
    //return (!is_regular(c)) && (B <= size(c) && size(c) <= 2*B);
    return (!is_regular(c)) && (size(c) <= 2*B);
  }

  static bool will_be_compressed(size_t sz) {
    return (sz >= B) && (sz <= 2*B);
  }

  static bool will_be_compressed(node* l, node* r, node* join) {
    assert(size(join) == 1);
    size_t sz = size(l) + size(r) + size(join);
    return will_be_compressed(sz);
  }

  static node* finalize(node* root) {
    auto sz = size(root);
    assert(sz > 0);
    if (sz < B) {
      auto ret = make_compressed_node(root);
      decrement_recursive(root);
      return ret;
    }
    return root;
  }

  // takes a pointer to an array of ETs, and a length of the number of ETs to
  // construct, and returns a compressed node.
  static compressed_node* make_single_compressed_node(ET* e, size_t s) {
//    assert(s >= B);
    assert(s <= 2*B);

    size_t encoded_size = EntryEncoder::encoded_size(e, s);
    size_t node_size = 3*sizeof(node_size_t) + encoded_size;
    compressed_node* c_node = (compressed_node*)utils::new_array_no_init<uint8_t>(node_size);

    c_node->r = 1;
    c_node->s = s;
    c_node->size_in_bytes = node_size;

    uint8_t* encoded_data = (((uint8_t*)c_node) + 3*sizeof(node_size_t));
    EntryEncoder::encode(e, s, encoded_data);

    check_compressed_node(c_node);
    return c_node;
  }

  static node* make_compressed(node* l, node* r, regular_node* e) {
    return basic_node_helpers::make_compressed<basic>(l, r, e);
  }

  static node* make_compressed(ET* stack, size_t tot) {
    return basic_node_helpers::make_compressed<basic>(stack, tot);
  }

  static node* make_compressed(node* l, node* r) {
    return make_compressed(l, r, nullptr);
  }

  // TODO: change return type to void.
  static ET* compressed_node_elms(compressed_node* c, ET* tmp_arr) {
    uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t));
    size_t i = 0;
    auto f = [&] (const ET& et) {
      parlay::assign_uninitialized(tmp_arr[i++], et);
      // tmp_arr[i++] = et;
    };
    EntryEncoder::decode(data_start, c->s, f);
    // TODO: can optimize in the case of a simple node (return c->arr)
    return tmp_arr;
  }

  /* ================================ GC ================================== */

  // Handles both regular and compressed nodes.
  static void free_node(node* va) {
    if (is_regular(va)) {
      auto a = cast_to_regular(va);
      (a->entry).~ET();
      allocator::free(a);
    } else {
      auto c = cast_to_compressed(va);
      uint8_t* data_start = (((uint8_t*)c) + 3*sizeof(node_size_t));
      EntryEncoder::destroy(data_start, c->s);
      auto array_size = c->size_in_bytes;
      utils::free_array<uint8_t>((uint8_t*)va, array_size);
    }
  }

  // precondition: a is non-null, and caller has a reference count on a.
  static bool decrement_count(node* a) {
    if ((utils::fetch_and_add(&generic_node(a)->r, -1) & kLowBitMask) == 1) {
      free_node(a);
      return true;
    }
    return false;
  }

  // precondition: a is non-null
  static void increment_count(node* a, node_size_t i) {
    utils::write_add(&(generic_node(a)->r), i);
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
    if (is_regular(t)) {
      node* lsub = cast_to_regular(t)->lc;
      node* rsub = cast_to_regular(t)->rc;
      if (decrement(t)) {
        utils::fork_no_result(size(lsub) >= utils::node_limit,
                              [&]() { decrement_recursive(lsub); },
                              [&]() { decrement_recursive(rsub); });
      }
    } else {
      decrement(t);
    }
  }


  /* ================================  Debug routines ================================ */

  // F applied to entries.
  template <class F>
  static size_t size_in_bytes(node* a, const F& f) {
    size_t total = 0;
    if (!a) return total;
    if (is_compressed(a)) {
      auto c = cast_to_compressed(a);
      total += c->size_in_bytes;
      auto fn = [&] (const auto& et) {
        total += f(et);
      };
      iterate_seq(c, fn);
    } else {
      auto r = cast_to_regular(a);
      total += size_in_bytes(r->lc, f);
      total += size_in_bytes(r->rc, f);
      total += sizeof(regular_node);
      total += f(get_entry(r));
    }
    return total;
  }

  static void print_node_info(node* a, std::string node_name) {
//    if (!a) { std::cout << "Empty node!" << std::endl; return; }
//    if (is_regular(a)) {
//      auto r = cast_to_regular(a);
//      std::cout << "Node " << node_name << " is regular. Key = " << get_entry(a).first << " Size = " << size(a) << " left child = " << size(r->lc) << " right_child = " << size(r->rc) <<  std::endl;
//    } else {
//      auto c = cast_to_compressed(a);
//      std::cout << "Node " << node_name << " is compressed. Size = " << c->s << std::endl;
//    }
  }

  static void print_inorder_rec(node* a) {
//    if (!a) return;
//    if (is_regular(a)) {
//      auto r = cast_to_regular(a);
//      print_inorder_rec(r->lc);
//      std::cout << get_entry(r).first << ". Size(" << size(a) << ") ";
//      print_inorder_rec(r->rc);
//    } else {
//      auto c = cast_to_compressed(a);
//      std::cout << "[";
//      for (size_t i=0; i<c->s; i++) {
//        std::cout << c->arr[i].first;
//        if (i < c->s - 1) {
//          std::cout << " ";
//        }
//      }
//      std::cout << "] ";
//    }
  }

  static void print_inorder(node* a) {
    if (a) {
      print_inorder_rec(a);
    } else {
      std::cout << "Empty tree.";
    }
    std::cout << std::endl;
  }

  // TODO: fix when needed.
  static void print_rec_impl(node* a) {
//    if (!a) return;
//    if (is_regular(a)) {
//      auto r = cast_to_regular(a);
////      std::cout << get_entry(r).first << ". Size(" << size(a) << ") " << " RefCnt(" << ref_cnt(a) << ") " << std::endl;
//      print_rec_impl(r->lc);
//      print_rec_impl(r->rc);
//    } else {
//      auto c = cast_to_compressed(a);
////      std::cout << ((size_t)a) << " RefCnt(" << ref_cnt(a) << ") [ size = " << size(a) << "]" << std::endl;
////      uint8_t* data_start = (((uint8_t*)c) + 2*sizeof(node_size_t));
////      auto f = [&] (const ET& e) {
////        EntryEncoder::print_info(e);
////        std::cout << " ";
////      };
////      EntryEncoder::decode(data_start, c->s, f);
//////      for (size_t i=0; i<c->s; i++) {
//////        std::cout << c->arr[i].first;
//////        if (i < c->s - 1) {
//////          std::cout << " ";
//////        }
//////      }
////      std::cout << "] " << std::endl;;
//    }
  }

  static void print_rec(node* a) {
    if (a) {
      print_rec_impl(a);
    } else {
      std::cout << "Empty tree.";
    }
    std::cout << std::endl;
  }

};

}  // namespace cpam
