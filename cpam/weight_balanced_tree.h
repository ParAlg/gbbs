#pragma once
#include "balance_utils.h"

// *******************************************
//   WEIGHT BALANCED TREES
//   From the paper:
//   Just Join for Parallel Ordered Sets
//   G. Blelloch, D. Ferizovic, and Y. Sun
//   SPAA 2016
// *******************************************

namespace cpam {

struct weight_balanced_tree {

  // no additional balance criterial needed for weight balanced trees
  struct data { };

  // defines: node_join, is_balanced
  // inherits: update, make_node, single
  template<class Node>
  struct balance : Node {
    using node = typename Node::node;
    using regular_node = typename Node::regular_node;
    using t_utils = balance_utils<balance>;
    using ET = typename Node::ET;
    friend t_utils;

    static node* node_join(node* t1, node* t2, regular_node* k) {
      return t_utils::node_join(t1, t2, k);
    }

    static regular_node* regular_node_join(regular_node* t1, regular_node* t2, regular_node* k) {
      return t_utils::regular_node_join(t1, t2, k);
    }

    static inline bool is_balanced(regular_node* t) {
      bool lc_heavy = is_left_heavy(t->lc, t->rc);
      bool rc_heavy = is_left_heavy(t->rc, t->lc);
      // std::cout << "lc_heavy = " << lc_heavy << " rc_heavy = " << rc_heavy << std::endl;
      bool ok =  !t || !(lc_heavy || rc_heavy);
      // check tight balance constraint on nodes with compressed children
//      if (t) {
//        if (Node::is_compressed(t->lc) || Node::is_compressed(t->rc)) {
//          bool l_ok = is_single_rotation(t, 0);
//          bool r_ok = is_single_rotation(t, 1);
//          if (!r_ok || !l_ok) {
//            std::cout << "t weight =" << weight(t) << std::endl;
//            std::cout << "l weight =" << weight(t->lc) << std::endl;
//            std::cout << "r weight =" << weight(t->rc) << std::endl;
//          }
//          assert(l_ok);
//          assert(r_ok);
//          return ok && l_ok && r_ok;
//        }
//      }
      return ok;
    }

  static bool check_structure(node* b, size_t orig_tree_size = std::numeric_limits<size_t>::max()) {
  //return true;
    if (!b) return true;
    // assert(Node::ref_cnt(b) == 1);
    if (orig_tree_size == std::numeric_limits<size_t>::max()) orig_tree_size = Node::size(b);
    if (Node::is_regular(b)) {
      size_t sz = Node::size(b);
      if (!(sz >= 2*B || sz < B)) {
        assert(false);
        exit(-1);
      }
      if (orig_tree_size > 2*B && sz < B) {
        assert(false);
        exit(-1);
      }
      assert(sz >= 2*B || sz < B);
      assert(!(orig_tree_size > 2*B && sz < B));
      auto r = Node::cast_to_regular(b);
      bool ret = is_balanced(r);
      if (!ret) {
        std::cout << "Check balance failed" << std::endl;
        assert(false);
        exit(0);
      }
      ret &= check_structure(r->lc, orig_tree_size);
      ret &= check_structure(r->rc, orig_tree_size);
      return ret;
    } else {
      return Node::check_compressed_node(b);
    }
  }

  private:
    static constexpr double alpha = 0.29;
    static constexpr double ratio = alpha / (1 - alpha);
    static constexpr double beta  = (1-2*alpha) / (1-alpha);

    // following two needed by balance_utils
    static inline bool is_left_heavy(node* t1, node* t2) {
      return ratio * weight(t1) > weight(t2);
    }

    static inline bool is_left_heavy(size_t s_t1, size_t s_t2) {
      return ratio * (s_t1 + 1) > (s_t2 + 1);
    }

    static inline bool is_single_rotation(regular_node* t, bool dir) {
      double balance = dir ? weight(t->rc) : weight(t->lc);
      return balance <= beta * weight(t);
    }

    // used?
    static inline bool could_single_rotate(node* lc, node* rc) {
      size_t l_s = weight(lc);
      size_t r_s = weight(rc);
      size_t tot = (Node::size(lc) + Node::size(rc) + 1) + 1;
      std::cout << "beta = " << beta << " ls = " << l_s << " rs = " << r_s << " tot = " << tot << std::endl;
      return (l_s <= beta * tot) && (r_s <= beta * tot);
    }

    static double weight(node* t) {
      return (double) (Node::size(t) + 1);
    }
  };
};

}  // namespace cpam
