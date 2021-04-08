#pragma once

// *******************************************
//   BALANCE UTILS
//   Mostly generic across balancing schemes
//   although not all schemes use these functions
// *******************************************

namespace cpam {

template <class Node>
struct balance_utils {
  using node = typename Node::node;
  using regular_node = typename Node::regular_node;
  using GC = gc<Node>;

  // requires t and t->lc have ref_cnt == 1
  static regular_node* rotate_right(regular_node* t) {
    regular_node* root = Node::cast_to_regular(t->lc);
    node* rsub = root->rc;
    root->rc = t, t->lc = rsub;
    Node::update(root->rc);
    Node::update(root);
    return root;
  }

  // requires t and t->rc have ref_cnt == 1
  static regular_node* rotate_left(regular_node* t) {
    regular_node* root = Node::cast_to_regular(t->rc);
    node* lsub = root->lc;
    root->lc = t, t->rc = lsub;
    Node::update(root->lc);
    Node::update(root);
    return root;
  }

  // requires t and t->lc have ref_cnt == 1
  static regular_node* double_rotate_right(regular_node* t) {
    regular_node* l = Node::cast_to_regular(t->lc);
    regular_node* root = Node::cast_to_regular(GC::copy_if_needed(l->rc));
    l->rc = root->lc;
    Node::update(l);
    root->lc = l;
    t->lc = root;
    return rotate_right(t);
  }

  // requires t and t->rc have ref_cnt == 1
  static regular_node* double_rotate_left(regular_node* t) {
    regular_node* r = Node::cast_to_regular(t->rc);
    regular_node* root = Node::cast_to_regular(GC::copy_if_needed(r->lc));
    r->lc = root->rc;
    Node::update(r);
    root->rc = r;
    t->rc = root;
    return rotate_left(t);
  }

  static regular_node* regular_right_join(regular_node* t1, regular_node* t2,
                                          regular_node* k) {
    if (!Node::is_left_heavy(t1, t2)) return regular_balanced_join(t1, t2, k);
    regular_node* t = Node::cast_to_regular(GC::copy_if_needed(t1));
    t->rc = regular_right_join(Node::cast_to_regular(t->rc), t2, k);

    // rebalance if needed
    if (Node::is_left_heavy(t->rc, t->lc)) {
      if (Node::is_single_rotation(Node::cast_to_regular(t->rc), 0))
        t = rotate_left(t);
      else
        t = double_rotate_left(t);
    } else
      Node::update(t);
    return t;
  }

  static regular_node* regular_left_join(regular_node* t1, regular_node* t2,
                                         regular_node* k) {
    if (!Node::is_left_heavy(t2, t1)) return regular_balanced_join(t1, t2, k);
    regular_node* t = Node::cast_to_regular(GC::copy_if_needed(t2));
    t->lc = regular_left_join(t1, Node::cast_to_regular(t->lc), k);

    // rebalance if needed
    if (Node::is_left_heavy(t->lc, t->rc)) {
      if (Node::is_single_rotation(Node::cast_to_regular(t->lc), 1))
        t = rotate_right(t);
      else
        t = double_rotate_right(t);
    } else
      Node::update(t);
    return t;
  }

  // a version without rotation
  static regular_node* regular_balanced_join(regular_node* l, regular_node* r,
                                             regular_node* e) {
    e->lc = l;
    e->rc = r;
    Node::update(e);
    return e;
  }

  // main function
  static regular_node* regular_node_join(regular_node* t1, regular_node* t2,
                                         regular_node* k) {
    if (Node::is_left_heavy(t1, t2)) return regular_right_join(t1, t2, k);
    if (Node::is_left_heavy(t2, t1)) return regular_left_join(t1, t2, k);
    return regular_balanced_join(t1, t2, k);
  }

/* ============================ Compressed Join =============================== */

  static node* check_ret(node* ret) {
#ifdef DEBUG
      if (!Node::check_structure(ret)) {
        Node::print_rec(ret);
        assert(false);
      }
#endif
    return ret;
  }

  static node* right_join(node* t1_o, node* t2, regular_node* k) {
    if (Node::is_compressed(t1_o)) {
      return Node::make_compressed(t1_o, t2, k);
    }
    regular_node* t1 = Node::cast_to_regular(t1_o);
    if (!Node::is_left_heavy(t1, t2)) {
      auto ret = balanced_join(t1, t2, k);
#ifdef DEBUG
      if (!Node::check_structure(ret)) { assert(false); }
#endif
      return ret;
    }
    // If t1 will become imbalanced after updating its right subtree with the
    // recursive right join, and the right subtree is going to be a compressed
    // node, then just do a full rebalance.
    // TODO(check): the preemptive rebalancing to eliminate single rotations
    // actually
    // prohibits this, right? So I believe it can be removed. Ditto in
    // left_join.
    // size_t new_right_size = Node::size(t1->rc) + Node::size(t2) + 1;
    // std::cout << "new_right_size = " << new_right_size << std::endl;
    ////std::cout << "Right join: t1.size = " << Node::size(t1) << " t2.size = "
    ///<< Node::size(t2) << std::endl;
    // if (Node::is_left_heavy(new_right_size, Node::size(t1->lc)) &&
    //    Node::will_be_compressed(new_right_size)) {
    //  return Node::full_rebalance(t1, t2, k);
    //}

    // New logic: if |t1| + |t2| + 1 fits inside of a compressed node, just
    // create it.
    size_t tot = Node::size(t1) + Node::size(t2) + 1;
    if (Node::will_be_compressed(tot)) {
      return Node::make_compressed(t1, t2, k);
    }

    // Otherwise even if t1 is imbalanced after the join, we know that this call
    // can successfully perform a single rotation.
    regular_node* t = Node::cast_to_regular(GC::copy_if_needed(t1));
    t->rc = right_join(t->rc, t2, k);

    // rebalance if needed
    if (Node::is_left_heavy(t->rc, t->lc)) {
      assert(Node::is_regular(t->rc));
      if (Node::is_single_rotation(Node::cast_to_regular(t->rc), 0)) {
        t = rotate_left(t);
      } else {
        t = double_rotate_left(t);
      }
    } else {
      Node::update(t);
    }
#ifdef DEBUG
      if (!Node::check_structure(t)) { assert(false); }
#endif
    return t;
  }

  static node* left_join(node* t1, node* t2_o, regular_node* k) {
    // std::cout << "in left_join" << std::endl;
    if (Node::is_compressed(t2_o)) {
      return Node::make_compressed(t1, t2_o, k);
    }
    regular_node* t2 = Node::cast_to_regular(t2_o);
    if (!Node::is_left_heavy(t2, t1)) {
      // std::cout << "left no longer heavy; balanced join" << std::endl;
      auto ret = balanced_join(t1, t2, k);
#ifdef DEBUG
      if (!Node::check_structure(ret)) { assert(false); }
#endif
      return ret;
    }

    size_t tot = Node::size(t1) + Node::size(t2) + 1;
    if (Node::will_be_compressed(tot)) {
      return Node::make_compressed(t1, t2, k);
    }

    regular_node* t = Node::cast_to_regular(GC::copy_if_needed(t2));
    t->lc = left_join(t1, t->lc, k);

    // rebalance if needed
    if (Node::is_left_heavy(t->lc, t->rc)) {
      assert(Node::is_regular(t->lc));
      if (Node::is_single_rotation(Node::cast_to_regular(t->lc), 1)) {
        t = rotate_right(t);
      } else {
        t = double_rotate_right(t);
      }
    } else {
      Node::update(t);
    }
#ifdef DEBUG
      if (!Node::check_structure(t)) { assert(false); }
#endif
    return t;
  }

  static constexpr double wb_alpha = 0.29;
  static constexpr double wb_ratio = wb_alpha / (1 - wb_alpha);
  static constexpr double wb_beta = (1 - 2 * wb_alpha) / (1 - wb_alpha);

  // Assumes that e will be used. need to decrement if it is not (using gc)
  static node* balanced_join(node* l, node* r, regular_node* e) {
    size_t l_s = Node::size(l), r_s = Node::size(r);
    assert(!Node::is_left_heavy(l, r)); assert(!Node::is_left_heavy(r, l));

    size_t tot = l_s + r_s + 1;
    // upper_threshold \approx 3.45B
    // size_t upper_threshold = B*(1 + (static_cast<double>(1) / wb_ratio));
    // 2.45B + 1
    //size_t upper_threshold = ceil(B * (static_cast<double>(1) / wb_ratio)) + 1;
    // size_t upper_threshold = ceil(B*(1 + (static_cast<double>(1) / wb_ratio)));

    // First check conditions for compressing.
    // If either left or right side is regular, and the overall tree will be at
    // least B, re-make the whole thing (the left  or right side needs to be
    // compressed anyway).
    if ((l_s < B || r_s < B) && (tot >= B)) {
      auto ret = Node::make_compressed(l, r, e);
      if (Node::is_regular(ret)) {
        assert(Node::is_balanced(Node::cast_to_regular(ret)));
      }
      return ret;
    }

    e->lc = l;
    e->rc = r;
    Node::update(e);

    assert(Node::is_balanced(e));
    return e;
  }

  static node* regular_join(regular_node* t1, regular_node* t2,
                            regular_node* k) {
    size_t tot = Node::size(t1) + Node::size(t2) + 1;
    size_t upper_threshold = ceil(B*(1 + (static_cast<double>(1) / wb_ratio)));
    if (tot >= B && tot <= upper_threshold) {
      return Node::make_compressed(t1, t2, k);
    }
    // check sizes of children
    assert(Node::size(t1) >= 2 * B || Node::size(t1) < B);
    assert(Node::size(t2) >= 2 * B || Node::size(t2) < B);

    if (Node::is_left_heavy(t1, t2)) return right_join(t1, t2, k);
    if (Node::is_left_heavy(t2, t1)) return left_join(t1, t2, k);
    return balanced_join(t1, t2, k);
  }


  // main function
  // note that if k is not used (e.g., it gets put into a compressed node, we
  // need to free it)
  //
  // Similarly, calls pass ownership of t1 and t2 to node_join here. In
  // particular, if node_join does not use t1 or t2 (e.g., if it decides to
  // return a compressed node), it will call decrement(t1) and decrement(t2).
  static node* node_join(node* t1, node* t2, regular_node* k) {
    bool t1_reg = Node::is_regular(t1), t2_reg = Node::is_regular(t2);
#ifdef DEBUG
      // std::cout << "t1 is: " << std::endl; Node::print_rec(t1);
      // std::cout << "t2 is: " << std::endl; Node::print_rec(t2);
      if (!Node::check_structure(t1)) { assert(false); }
      if (!Node::check_structure(t2)) { assert(false); }
#endif

    if (t1_reg && t2_reg) {
      auto ret =
          regular_join(Node::cast_to_regular(t1), Node::cast_to_regular(t2),
                       Node::cast_to_regular(k));
#ifdef DEBUG
      if (!Node::check_structure(ret)) {
        Node::print_rec(ret);
//        size_t ret_size = Node::size(ret);
        assert(false);
      }
#endif
      return ret;
    } else if (t1_reg) {  // t2 compressed
      if (Node::is_left_heavy(t1, t2)) {
        auto ret = right_join(Node::cast_to_regular(t1), t2, k);
        return check_ret(ret);
      } else {
        if (Node::is_left_heavy(t2, t1)) {
          auto ret = Node::make_compressed(t1, t2, k);
          return check_ret(ret);
        } else {
          auto ret = balanced_join(t1, t2, k);
          return check_ret(ret);
        }
      }
    } else if (t2_reg) {  // t1 compressed
      if (Node::is_left_heavy(t2, t1)) {
        auto ret = left_join(t1, Node::cast_to_regular(t2), k);
        return check_ret(ret);
      } else {
        if (Node::is_left_heavy(t1, t2)) {
          auto ret = Node::make_compressed(t1, t2, k);
          return check_ret(ret);
        } else {
          auto ret = balanced_join(t1, t2, k);
          return check_ret(ret);
        }
      }
    }
    // both compressed
    bool balanced =
        !(Node::is_left_heavy(t1, t2) || Node::is_left_heavy(t2, t1));
    if (balanced || !Node::could_single_rotate(t1, t2)) {
      // std::cout << "Performing a balanced join from node_join" << std::endl;
      auto ret = balanced_join(t1, t2, k);
      return check_ret(ret);
    }
    // std::cout << "Performing a compressed join from node_join" << std::endl;
    auto ret = Node::make_compressed(t1, t2, k);
    return check_ret(ret);
  }

};

}  // namespace cpam
