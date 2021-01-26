#pragma once
#include "parlay/alloc.h"

using node_size_t = unsigned int;
//using node_size_t = size_t;

// *******************************************
//   BASIC NODE
// *******************************************

template<class balance, class _ET>
struct basic_node {
  using ET = _ET;

  struct node : balance {
    node* lc;  node* rc;
    ET entry;
    node_size_t s;
    node_size_t ref_cnt;
  };

  using allocator = parlay::type_allocator<node>;

  static node_size_t size(node* a) {
    return (a == NULL) ? 0 : a->s;
  }

  static void update(node* a) {
    a->s = size(a->lc) + size(a->rc) + 1;
  }

  static node* make_node(ET e) {
    node *o = allocator::alloc();
    o->ref_cnt = 1;
    //o->entry = e;
    parlay::assign_uninitialized(o->entry,e);
    return o;
  }

  static node* single(ET e) {
    node* r = make_node(e);
    r->lc = r->rc = NULL; r->s = 1;
    return r;
  }

  static void free_node(node* a) {
	(a->entry).~ET();
    allocator::free(a);
  }

  static node* empty() {return NULL;}
  inline static ET& get_entry(node *a) {return a->entry;}
  inline static ET* get_entry_p(node *a) {return &(a->entry);}
  static void set_entry(node *a, ET e) {a->entry = e;}
  static node* left(node a) {return a.lc;}
  static node* right(node* a) {return a.rc;}
};

//template <class E>
//using node_type = basic_node<weight_balanced_tree, E>;
