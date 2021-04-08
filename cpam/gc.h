#pragma once
#include "basic_node.h"
#include "utils.h"

// *******************************************
//   GC
// *******************************************

namespace cpam {

template <class Node>
struct gc {
  using node = typename Node::node;
  using regular_node = typename Node::regular_node;
  using compressed_node = typename Node::compressed_node;
  using alloc = typename Node::allocator;
  using ET = typename Node::ET;

  static void init() { alloc::init();}
  static bool initialized() { return alloc::initialized;}
  static void reserve(size_t n, bool randomize = false) {
    alloc::reserve(n);}
  static void finish() { alloc::finish();}
  static size_t num_used_nodes() {return alloc::num_used_blocks();}

  static void decrement(node* t) {
    Node::decrement(t);
  }

  static void decrement_recursive(node* t) {
    Node::decrement_recursive(t);
  }

//  // atomically decrement ref count and delete node if zero
//  static bool decrement(node* t) {
//    if (t) {
//      if (utils::fetch_and_add(&t->ref_cnt, -1) == 1) {
//	Node::free_node(t);
//	return true;
//      }
//    }
//    return false;
//  }
//
//  // atomically decrement ref count and if zero:
//  //   delete node and recursively decrement the two children
//  static void decrement_recursive(node* t) {
//    if (!t) return;
//    node* lsub = t->lc;
//    node* rsub = t->rc;
//    if (decrement(t)) {
//      utils::fork_no_result(Node::size(lsub) >= utils::node_limit,
//         [&]() {decrement_recursive(lsub);},
//         [&]() {decrement_recursive(rsub);});
//    }
//  }

  // atomically increment the reference count
  static void increment(node* t) {
    if (t) {
      Node::increment_count(t, 1);
    }
  }

  static void increment(node* t, int i) {
    if (t) {
      Node::increment_count(t, i);
    }
  }

  static inline node* inc(node* t) {
    increment(t);
    return t;
  }

//  // atomically increment the reference count
//  static void increment(node* t) {
//    if (t) { utils::write_add(&t->ref_cnt, 1);}
//  }
//
//  static inline node* inc(node* t) {
//    increment(t); return t;
//  }

  // copies node with entry and children, incrementing children's ref counts
  // does not update
  static inline node* copy(node* t) {
    if (Node::is_regular(t)) {
      auto o = Node::make_regular_node(Node::get_entry(t));
      o->lc = inc(Node::cast_to_regular(t)->lc);
      o->rc = inc(Node::cast_to_regular(t)->rc);
      return o;
    } else {
      auto o = Node::make_compressed_node(t);
      return o;
    }
  }

  // copy node if reference count is > 1, incrementing the children's ref counts
  // then decrements the copied node's count
  // does not update
  static inline node* copy_if_needed(node* t) {
    if (!t) std::cout << "copy if needed fail" << std::endl;
    node* res = t;
    if (Node::ref_cnt(t) > 1) {
      res = copy(t);
      decrement_recursive(t);
    }
    return res;
  }

  // increments ref count if copy flag is true
  static inline node* inc_if(node* t, bool copy) {
    if (copy) increment(t);
    return t;
  }

  // decrements if there is no extra pointer to t
  // if t is a copy then t needs to be recursively decremented
  // if t is not a copy then t's children should already be decremented
  static inline void dec_if(node* t, bool copy, bool extra_ptr) {
    if (!extra_ptr) {
      if (copy) decrement_recursive(t);
      else decrement(t);
    }
  }

  static size_t used_node() {
    if (!alloc::initialized) {
      std::cout << "allocator not initialized" << std::endl;
      return 0;
    } else {
      return alloc::num_used_blocks();
    }
  }

  static void print_stats() {
    if (!alloc::initialized)
      std::cout << "allocator not initialized" << std::endl;
    else {
      size_t used = alloc::num_used_blocks();
      size_t allocated = alloc::num_allocated_blocks();
      size_t size = alloc::block_size();
      std::cout << "used: " << used << ", allocated: " << allocated
		<< ", node size: " << size
		<< ", bytes: " << size*allocated << std::endl;
    }
  }

  // TODO:used?
  struct read_ptr {
    read_ptr(node* p) : p(p) {}
    node* p;

    bool empty() { return !p; }
  };

  struct ptr {
    // Constructs a ptr object with the given node pointer.
    // extra_ptr = true if there is an implicit extra pointer to this node, but
    // we don't want to bump the ref-cnt to encode this explicitly. Instead, we
    // store a bi to remember that there's an extra pointer.
    ptr(node* a, bool extra_ptr = false) : p(set_flag(a, extra_ptr)) {
      // cout << "construct: " << a << ", " << extra_ptr << endl;
    }

    ptr(const ptr& a) {
      std::cout  << "Extra ptr copy " << a.p << std::endl;
      assert(false);
      p = strip_flag(a.p);
      increment(p, 1 + has_flag(a.p));
    }

    ptr(ptr&& a) {
      // cout << "move: " << a.p << endl;
      p = a.p;
      a.p = nullptr;
    }

    ptr& operator=(ptr a) {
      // cout << "assign" << endl;
      swap(p, a.p);
      return *this;
    }

    node_size_t ref_cnt() {
      return Node::ref_cnt(strip_flag(p));
    }

    bool empty() { return (strip_flag(p) == nullptr); }

    // Equivalent to GC::inc_if(r, extra_ptr);
    // This destructs the underlying pointer.
    node* node_ptr() {
      if (has_flag(p)) increment(strip_flag(p));
      node* r = strip_flag(p);
      p = nullptr;
      return r;
    }

    node* unsafe_ptr() {
      return strip_flag(p);
    }

    bool requires_copy() {
      return has_flag(p) || Node::ref_cnt(strip_flag(p)) > 1;
    }

    size_t size() const { return Node::size(strip_flag(p)); }

    ET entry() {
      return Node::get_entry(strip_flag(p));
    }

    ~ptr() {
      // cout << "destruct: " << p << endl;
      if (p && !has_flag(p)) {
//        cout << "decrement recursive from : " << ((size_t)p) << endl;
        decrement_recursive(p);
      }
    }

    /* Should be private */
    node* p;

    static node* set_flag(node* x, bool f) { return (node*)(((size_t)x) + f); }
    static bool has_flag(node* x) { return ((size_t)x) & 1; }
    static node* strip_flag(node* x) { return (node*)(((size_t)x) & ~((size_t)1)); }
    bool extra() { return has_flag(p); }

    bool is_compressed() { return Node::is_compressed(strip_flag(p)); }
    bool is_regular() { return Node::is_regular(strip_flag(p)); }
  };

};

}  // namespace cpam
