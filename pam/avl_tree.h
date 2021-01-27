#pragma once
#include "balance_utils.h"

// *******************************************
//   AVL TREES
//   From the paper:
//   Just Join for Parallel Ordered Sets
//   G. Blelloch, D. Ferizovic, and Y. Sun
//   SPAA 2016
// *******************************************

struct avl_tree {

  struct data { int height; };

  // defines: node_join, is_balanced
  // redefines: update, single
  // inherits: make_node
  template<class Node>
  struct balance : Node {
    using node = typename Node::node;
    using t_utils = balance_utils<balance>; // bulk of work done here
    friend t_utils;
    using ET = typename Node::ET;

    static node* node_join(node* t1, node* t2, node* k) {
      return t_utils::node_join(t1, t2, k);
    }

    static inline bool is_balanced(node* t) {
      //return !t || !(is_left_heavy(t->lc,t->rc) || is_left_heavy(t->rc,t->lc));
	  return !t || !((is_left_heavy(t->lc,t->rc) || is_left_heavy(t->rc,t->lc)));
    }

    static void update(node* t) {
      Node::update(t);
      t->height = std::max(height(t->lc),height(t->rc)) + 1;
    }

    static node* single(ET e) {
      node *o = Node::single(e);
      o->height = 1;
      return o;
    }

  private:
    static int height(node* a) {
      if (a == NULL) return 0;
      else return a->height;
    }

    // following two needed by balance_utils
    static inline bool is_left_heavy(node* t1, node* t2) {
      return height(t1) > height(t2) + 1;
    }

    static inline bool is_single_rotation(const node* t, const bool dir) {
      bool heavier = height(t->lc) > height(t->rc);
      return dir ? heavier : !heavier;
    }

	// static inline bool is_single_rotation(const node* t, const bool dir) {
      // return dir ? (height(t->lc) > height(t->rc)) : (height(t->lc) < height(t->rc));
    // }


  };
};
