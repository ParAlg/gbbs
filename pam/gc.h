#pragma once
#include "basic_node.h"
#include "utils.h"

// *******************************************
//   GC
// *******************************************

template <class Node>
struct gc {
  using node = typename Node::node;
  using alloc = typename Node::allocator;

  static void init() { alloc::init();}
  static bool initialized() { return alloc::initialized;}
  static void reserve(size_t n, bool randomize = false) {
    alloc::reserve(n);}
  static void finish() { alloc::finish();}
  static size_t num_used_nodes() {return alloc::num_used_blocks();}

  // atomically decrement ref count and delete node if zero
  static bool decrement(node* t) {
    if (t) {
      if (utils::fetch_and_add(&t->ref_cnt, -1) == 1) {
	Node::free_node(t);
	return true;
      }
    }
    return false;
  }

  // atomically decrement ref count and if zero:
  //   delete node and recursively decrement the two children
  static void decrement_recursive(node* t) {
    if (!t) return;
    node* lsub = t->lc;
    node* rsub = t->rc;
    if (decrement(t)) {
      utils::fork_no_result(Node::size(lsub) >= utils::node_limit,
         [&]() {decrement_recursive(lsub);},
         [&]() {decrement_recursive(rsub);});
    }
  }

  // atomically increment the reference count
  static void increment(node* t) {
    if (t) { utils::write_add(&t->ref_cnt, 1);}
  }

  static inline node* inc(node* t) {
    increment(t); return t;
  }

  // copies node with entry and children, incrementing children's ref counts
  // does not update
  static inline node* copy(node* t) {
    node* o = Node::make_node(Node::get_entry(t));
    o->lc = inc(t -> lc);
    o->rc = inc(t -> rc);
    return o;
  }

  // copy node if reference count is > 1, incrementing the children's ref counts
  // then decrements the copied node's count
  // does not update
  static inline node* copy_if_needed(node* t) {
	if (!t) std::cout << "copy if needed fail" << std::endl;
    node* res = t;
    if (t->ref_cnt > 1) {
      res = copy(t);
      decrement_recursive(t);
    }
    return res;
  }

  // This is a "lazy" version of copy_if_needed
  // It does not increment the children (or even copy pointers to them)
  // It also only decrements the old version if extra_ptr==false
  // Significantly improves performance since it avoids ref_cnt updates
  // Be very careful of placement:
  // needs to go after last use of t or t's children, since it can
  // delete t, and its descendants
  // typical usage:
  //      bool copy = (t->ref_cnt > 1) || extra_ptr
  //      recursive calls on left and right children passing copy as extra_ptr
  //    or inc_if(t->xc, copy) if no recursive call on child xc (x in {l,r})
  //      node* o = copy_if(t, copy, extra_ptr)
  //      return node_join(l, r, o) -- needs to set children of o, and update
  //    or dec_if(t, copy extra_ptr), if t not needed
  // important to only read ref_cnt once, so it is atomic.
  static inline node* copy_if(node* t, bool copy, bool extra_ptr) {
    if (copy) {
      node* r = Node::make_node(Node::get_entry(t));
      if (!extra_ptr) decrement_recursive(t);
      return r;
    } else return t;
  }

template <typename E>
  static inline node* copy_if(node* t, E e, bool copy, bool extra_ptr) {
    if (copy) {
      node* r = Node::make_node(e);
      if (!extra_ptr) decrement_recursive(t);
      return r;
    } else return t;
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
};
