#pragma once
#include "balance_utils.h"

// *******************************************
//   TREAPS
// *******************************************

template<class Entry>
struct treap {

  struct data { };

  // defines: node_join, is_balanced
  // inherits: update, make_node, single
  template<class Node>
  struct balance : Node {
    using node = typename Node::node;
    using t_utils = balance_utils<balance>;
    using ET = typename Node::ET;
    using GC = gc<balance>;

    static inline bool is_balanced(node* t) {
      return !t || (priority(t) >= priority(t->lc) && priority(t) >= priority(t->rc));
    }

    static node* node_join(node* t1, node* t2, node* k) {
      if (priority(k) >= priority(t1) && priority(k) >= priority(t2)) {
	k->lc = t1; k->rc = t2;
	Node::update(k);
	return k;
      } else if (priority(t1) > priority(t2)) {
	node* t = GC::copy_if_needed(t1);
	t->rc = node_join(t->rc, t2, k);
	Node::update(t);
	return t;
      } else {
	node* t = GC::copy_if_needed(t2);
	t->lc = node_join(t1, t->lc, k);
	Node::update(t);
	return t;
      }
    }

  private:
    static size_t priority(node* a) {
      return (a == NULL) ? 0 :Entry::hash(Node::get_entry(a));
    }
  };
};
