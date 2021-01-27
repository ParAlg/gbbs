#pragma once
#include "balance_utils.h"

// *******************************************
//   RED BLACK TREES
//   From the paper:
//   Just Join for Parallel Ordered Sets
//   G. Blelloch, D. Ferizovic, and Y. Sun
//   SPAA 2016
// *******************************************


struct red_black_tree {

  typedef enum { RED, BLACK } Color;
  //static const unsigned char RED = 0;
  //static const unsigned char BLACK = 1;


  struct data { unsigned char height; Color color;};
  //using Color = unsigned char;
  //struct data { unsigned char height; unsigned char color = 1;};

  // defines: node_join, is_balanced
  // redefines: update, single
  // inherits: make_node
  template<class Node>
  struct balance : Node {
    using node = typename Node::node;
    using t_utils = balance_utils<balance>;
    using ET = typename Node::ET;
    using GC = gc<balance>;

    static bool is_balanced(node* t) {
      // return (!t ||
	      // (((t->color == BLACK) ||
		// (color(t->lc) == BLACK && color(t->rc) == BLACK))
	       // && (height(t->lc) == height(t->rc))));
	  if (!t) return true;
	  bool f = true;
	  if (!((t->color == BLACK) ||
		(color(t->lc) == BLACK && color(t->rc) == BLACK))) {
			std::cout << "bad color" << std::endl;
			std::cout << color(t->lc) << " " << color(t) << " " << color(t->rc) << std::endl;
			f = false;
		}
	  if  (!(height(t->lc) == height(t->rc))) {
		  std::cout << "bad height" << std::endl;
		  f = false;
	  }
	  return f;
    }

    static void update(node* t) {
      Node::update(t);
      t->height = height(t->lc) + (t->color == BLACK);
    }

	static node* make_node(ET e) {
		node *o = Node::make_node(e);
		o->height = 2;
		o->color = BLACK;
		return o;
	}

    static node* single(ET e) {
      node *o = Node::single(e);
      o->height = 2;
      o->color = BLACK;
      return o;
    }

    // main function
    static node* node_join(node* t1, node* t2, node* k) {
      //right join
	  if (height(t1) > height(t2)) {
		node* t = right_join(t1, t2, k);

	    if (t->color == BLACK && color(t->rc) == RED && color(t->rc->rc) == RED) {
			t->rc->rc->color = BLACK;
			t->rc->rc->height++;
			t = t_utils::rotate_left(t);
        } else update(t);

		if (t->color == RED && color(t->rc) == RED) {
		  t->color = BLACK; t->height++;}
		return t;
      }

	  //left join
      if (height(t2) > height(t1)) {
		node* t = left_join(t1, t2, k);

		if (t->color == BLACK && color(t->lc) == RED && color(t->lc->lc) == RED) {
			t->lc->lc->color = BLACK;
			t->lc->lc->height++;
			t = t_utils::rotate_right(t);
        } else update(t);

		if (t->color == RED && color(t->lc) == RED) {
		  t->color = BLACK; t->height++;}
		return t;
      }

	  //balanced join
      if (color(t1) == BLACK && color(t2) == BLACK)
		return balanced_join(t1,t2,k,RED);
      return balanced_join(t1,t2,k,BLACK);
    }

  private:

    static int height(node* a) {
      return (a == NULL) ? 1 : a->height;
    }

    static int color(node* a) {
      return (a == NULL) ? BLACK : a->color;
    }

    // a version without rotation
    static node* balanced_join(node* l, node* r, node* e, Color color) {
      e->color = color;
      e->lc = l;
      e->rc = r;
      update(e);
      return e;
    }

    static node* right_join(node* t1, node* t2, node* k) {
		// if (!t1) {
			// if (!t2) std::cout << "null and null" << std::endl; else {
				// std::cout << "null and non-null" << std::endl;
				// std::cout << height(t1) << " " << height(t2) << std::endl;
			// }
		// }
      if (height(t1) == height(t2) && color(t1) == BLACK)
		return balanced_join(t1, t2, k, RED);
	  //if (!t1) std::std::cout << "right_join fail" << std::std::endl;
      node* t = GC::copy_if_needed(t1);
      t->rc = right_join(t->rc, t2, k);

      // rebalance if needed
      if (t->color == BLACK && color(t->rc) == RED && color(t->rc->rc) == RED) {
		t->rc->rc->color = BLACK;
		t->rc->rc->height++;
		t = t_utils::rotate_left(t);
      } else update(t);
      return t;
    }

    static node* left_join(node* t1, node* t2, node* k) {
		// if (!t2) {
			// if (!t1) std::cout << "null and null" << std::endl; else {
				// std::cout << "null and non-null" << std::endl;
				// std::cout << height(t1) << " " << height(t2) << std::endl;
			// }
		// }
      if (height(t1) == height(t2) && color(t2) == BLACK)
		return balanced_join(t1, t2, k, RED);
	//if (!t2) std::std::cout << "left_join fail" << std::std::endl;
      node* t = GC::copy_if_needed(t2);
      t->lc = left_join(t1, t->lc, k);

      // rebalance if needed
      if (t->color == BLACK && color(t->lc) == RED && color(t->lc->lc) == RED) {
	t->lc->lc->color = BLACK;
	t->lc->lc->height++;
	t = t_utils::rotate_right(t);
      } else update(t);
      return t;
    }

  };
};
