#pragma once

#include <assert.h>

#include "utils.h"

namespace cpam {
namespace basic_node_helpers {

static constexpr size_t B = utils::compression_block_size;

template <class Node>
static typename Node::node* single_compressed_node(typename Node::ET* stack, size_t tot) {
  auto c = Node::make_single_compressed_node(stack, tot);
  Node::check_compressed_node(c);
  return c;
}

template <class Node>
static typename Node::node* two_compressed_nodes(typename Node::ET* stack, size_t tot,
                                  typename Node::regular_node* e = nullptr) {
  assert(tot >= (2 * B + 1));
  size_t left_size = tot / 2, right_size = tot / 2 - (!(tot & 1));
  if (e == nullptr) e = Node::single(stack[left_size]);

  auto c_l = Node::make_single_compressed_node(stack, left_size);
  Node::set_entry(e, stack[left_size]);
  auto c_r = Node::make_single_compressed_node(stack + left_size + 1, right_size);

  e->lc = c_l;
  e->rc = c_r;
  Node::update(e);
  return e;
}

  template <class Node>
  static typename Node::node* four_compressed_nodes(typename Node::ET* stack, size_t tot,
                                     typename Node::regular_node* e = nullptr) {
    using ET = typename Node::ET;
    using node = typename Node::node;
    auto make_two_nodes = [](ET* start, size_t tot) -> node* {
      // return two_compressed_nodes(start, tot);
      //  TODO: refactor into above
      assert(tot >= 2*B);
      if (tot == 2*B) {
        return Node::make_single_compressed_node(start, tot);
      } else {
        size_t lc_size = tot/2;
        size_t rc_size = tot/2;
        if (tot % 2 == 0) lc_size--;
        auto c_l = Node::make_single_compressed_node(start, lc_size);
        auto c_r = Node::make_single_compressed_node(start + lc_size + 1,
        rc_size);
        auto e = Node::make_regular_node(start[lc_size]);
        e->lc = c_l;
        e->rc = c_r;
        Node::update(e);
        return e;
      }
    };

  auto lc_size = tot / 2, rc_size = tot / 2;
  if (tot % 2 == 0) lc_size--;
  node* lc = make_two_nodes(stack, lc_size);
  node* rc = make_two_nodes(stack + lc_size + 1, rc_size);

  if (!e) e = Node::make_regular_node(stack[lc_size]);
  Node::set_entry(e, stack[lc_size]);  // double set?

  e->lc = lc;
  e->rc = rc;
  Node::update(e);
  return e;
}

template <class Node>
static typename Node::node* make_compressed(typename Node::ET* stack, size_t tot) {
  assert(tot >= B);
  assert(tot <= 7 * B);
  using ET = typename Node::ET;

  typename Node::node* ret = nullptr;
  if (tot <= 2 * B) {
    ret = single_compressed_node<Node>((ET*)stack, tot);
  } else if (tot <= 4 * B + 1) {
    ret = two_compressed_nodes<Node>(stack, tot);
  } else {
    // 2B on lhs, 4.9B on rhs
    assert(tot <= 7 * B);
    LOG("Returning four nodes, tot = " << tot);
    ret = four_compressed_nodes<Node>(stack, tot);
  }

  return ret;
}

// owns e
// will decrement l and r.
template <class Node>
static typename Node::node* make_compressed(typename Node::node* l, typename Node::node* r, typename Node::regular_node* e) {
  size_t l_s = Node::size(l);
  size_t r_s = Node::size(r);
  size_t tot = l_s + r_s + (e != nullptr);
  assert(tot >= B);

  using ET = typename Node::ET;
  ET stack[7 * B];
  size_t offset = 0;
  auto copy_f = [&](const ET& a) {
    parlay::assign_uninitialized(stack[offset++], a);
  };
  Node::iterate_seq(l, copy_f);
  if (e) {
    copy_f(Node::get_entry(e));
  }
  Node::iterate_seq(r, copy_f);
  assert(offset == tot);

  void* ret;
  if (tot <= 2*B) {
    ret = single_compressed_node<Node>((ET*)stack, offset);
    if (e) Node::decrement(e);
  } else if (tot <= 4*B + 1) {
    ret = two_compressed_nodes<Node>(stack, offset, e);
  } else {
    // 2B on lhs, 4.9B on rhs
    assert(tot <= 7*B);
    LOG("Returning four nodes, tot = " << tot);
    ret = four_compressed_nodes<Node>(stack, offset, e);
  }
  Node::decrement_recursive(l);
  Node::decrement_recursive(r);
  return ret;
}

// Used by GC to copy a compressed node. TODO: update to work correctly with
// diff-encoding.
template <class Node>
static typename Node::node* make_compressed_node(typename Node::node* b) {
  using ET = typename Node::ET;
  ET stack[2*B];
  size_t offset = 0;
  auto copy_f = [&] (const ET& e) {
    parlay::assign_uninitialized(stack[offset++], e);
    // stack[offset++] = e;
  };
  Node::iterate_seq(b, copy_f);
  return Node::make_single_compressed_node(stack, offset);
}

}  // namespace basic_node_helpers
}  // namespace cpam
