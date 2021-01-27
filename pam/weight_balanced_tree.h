#pragma once
#include "balance_utils.h"

// *******************************************
//   WEIGHT BALANCED TREES
//   From the paper:
//   Just Join for Parallel Ordered Sets
//   G. Blelloch, D. Ferizovic, and Y. Sun
//   SPAA 2016
// *******************************************

struct weight_balanced_tree {

  // no additional balance criterial needed for weight balanced trees
  struct data { };

  // defines: node_join, is_balanced
  // inherits: update, make_node, single
  template<class Node>
  struct balance : Node {
    using node = typename Node::node;
    using t_utils = balance_utils<balance>;
    using ET = typename Node::ET;
    friend t_utils;

    static node* node_join(node* t1, node* t2, node* k) {
      return t_utils::node_join(t1, t2, k);
    }

    static inline bool is_balanced(node* t) {
      return !t || !(is_left_heavy(t->lc,t->rc) || is_left_heavy(t->rc,t->lc));
    }

  private:
    static constexpr double alpha = 0.29;
    static constexpr double ratio = alpha / (1 - alpha);
    static constexpr double beta  = (1-2*alpha) / (1-alpha);

    // following two needed by balance_utils
    static inline bool is_left_heavy(node* t1, node* t2) {
      return ratio * weight(t1) > weight(t2);
    }

    static inline bool is_single_rotation(node* t, bool dir) {
      double balance = dir ? weight(t->rc) : weight(t->lc);
      return balance <= beta * weight(t);
    }

    static double weight(node* t) {
      return (double) (Node::size(t) + 1);
    }
  };
};
