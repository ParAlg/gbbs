// This code is part of the project "A Simple Parallel Cartesian Tree
// Algorithm and its Application to Parallel Suffix Tree
// Construction", ACM Transactions on Parallel Computing, 2014.
// (earlier version appears in ALENEX 2011).
// Copyright (c) 2014-2019 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Parallel algorithm for Catesian trees
// Takes a sequence of integers and returns a parent sequence of the same length.
// Each location points to its parent location in the Cartesian tree.
// The tree is binary: ties are broken arbitrarily
// The root points to itself
// If a parent points to the left, then it is a right child
//    and if it points to the right it is a left child
#pragma once

namespace pbbs {

  template <class Index>
  struct ct_node { Index value; Index parent;};

  template <class Index>
  void ct_merge(ct_node<Index>* N, Index left, Index right) {
    Index head;
    if (N[left].value > N[right].value) {
      head = left; left = N[left].parent;}
    else {head = right; right= N[right].parent;}

    while(1) {
      if (left == 0) {N[head].parent = right; break;}
      if (right == 0) {N[head].parent = left; break;}
      if (N[left].value > N[right].value) {
	N[head].parent = left; left = N[left].parent;}
      else {N[head].parent = right; right = N[right].parent;}
      head = N[head].parent;
    }
  }

  template <class Index>
  void cartesian_tree_r(ct_node<Index>* Nodes, Index s, Index e) {
    if (e-s < 2) {
    } else if (e-s == 2) {
      if (Nodes[s].value > Nodes[s+1].value)
	Nodes[s].parent=s+1;
      else Nodes[s+1].parent=s;
    } else {
      Index mid = (s+e)/2;
      par_do_if((e-s) > 1000,
		[&] () {cartesian_tree_r(Nodes, s, mid);},
		[&] () {cartesian_tree_r(Nodes, mid, e);});
      ct_merge(Nodes, mid-1, mid);
    }
  }

  template <class Index>
  sequence<Index> cartesian_tree(sequence<Index> const &S) {
    size_t n = S.size();
    sequence<ct_node<Index>> Nodes(n, [&] (Index i) {
	ct_node<Index> r = {S[i], 0};
	return r; });
    cartesian_tree_r(Nodes.begin(), (Index) 0, (Index) n);
    return sequence<Index>(n, [&] (size_t i) {
	return Nodes[i].parent;});
  }

}  // namespace pbbs
