#pragma once

// *******************************************
//   BALANCE UTILS
//   Mostly generic across balancing schemes
//   although not all schemes use these functions
// *******************************************

template <class Node>
struct balance_utils {
  using node = typename Node::node;
  using GC = gc<Node>;

  // requires t and t->lc have ref_cnt == 1
  static node* rotate_right(node* t) {
    node* root = t->lc;
    node* rsub = root->rc;
    root->rc = t, t->lc = rsub;
    Node::update(root->rc);
    Node::update(root);
    return root;
  }

  // requires t and t->rc have ref_cnt == 1
  static node* rotate_left(node* t) {
    node* root = t->rc;
    node* lsub = root->lc;
    root->lc = t, t->rc = lsub;
    Node::update(root->lc);
    Node::update(root);
    return root;
  }

  // requires t and t->lc have ref_cnt == 1
  static node* double_rotate_right(node* t) {
    node* l = t->lc;
    node* root = GC::copy_if_needed(l->rc);
    l->rc = root->lc;
    Node::update(l);
    root->lc = l;
    t->lc = root;
    return rotate_right(t);
  }

  // requires t and t->rc have ref_cnt == 1
  static node* double_rotate_left(node* t) {
    node* r = t->rc;
    node* root = GC::copy_if_needed(r->lc);
    r->lc = root->rc;
    Node::update(r);
    root->rc = r;
    t->rc = root;
    return rotate_left(t);
  }

  static node* right_join(node* t1, node* t2, node* k) {
    if (!Node::is_left_heavy(t1, t2)) return balanced_join(t1, t2, k);
    node* t = GC::copy_if_needed(t1);
    t->rc = right_join(t->rc, t2, k);

    // rebalance if needed
    if (Node::is_left_heavy(t->rc, t->lc)) {
      if (Node::is_single_rotation(t->rc, 0)) t = rotate_left(t);
      else t = double_rotate_left(t);
    } else Node::update(t);
    return t;
  }

  static node* left_join(node* t1, node* t2, node* k) {
    if (!Node::is_left_heavy(t2, t1)) return balanced_join(t1, t2, k);
    node* t = GC::copy_if_needed(t2);
    t->lc = left_join(t1, t->lc, k);

    // rebalance if needed
    if (Node::is_left_heavy(t->lc, t->rc)) {
      if (Node::is_single_rotation(t->lc, 1)) t = rotate_right(t);
      else t = double_rotate_right(t);
    } else Node::update(t);
    return t;
  }

  // a version without rotation
  static node* balanced_join(node* l, node* r, node* e) {
    e->lc = l;
    e->rc = r;
    Node::update(e);
    return e;
  }

  // main function
  static node* node_join(node* t1, node* t2, node* k) {
    if (Node::is_left_heavy(t1, t2)) return right_join(t1, t2, k);
    if (Node::is_left_heavy(t2, t1)) return left_join(t1, t2, k);
    return balanced_join(t1,t2,k);}

};
