// This code is part of the project "A Simple Parallel Cartesian Tree
// Algorithm and its Application to Parallel Suffix Tree
// Construction", ACM Transactions on Parallel Computing, 2014
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
#pragma once

#include "sequence.h"
#include "get_time.h"
#include "suffix_array.h"
#include "integer_sort.h"
#include "lcp.h"
#include "cartesian_tree.h"
#include "histogram.h"

namespace pbbs {

  template <class Uint>
  struct suffix_tree {
    static constexpr bool verbose = false;
    using uchar = unsigned char;
    enum etype {leaf, internal, empty};
    timer t;

    struct edge {
      Uint child;
      uchar c;
      etype type;
    };

    edge empty_edge = {0, 0, empty};

    using Pair = std::pair<Uint,Uint>;

    struct node {
      Uint lcp;
      Uint location;
      Uint offset;
    };

    sequence<edge> Edges;
    sequence<node> Nodes;
    sequence<Uint> LCP; // can be dropped after construction
    sequence<Uint> SA;
    sequence<uchar> S;
    Uint n;
    Uint num_edges;

    suffix_tree(sequence<uchar> const &Str) {
      t = timer("suffix tree", verbose);
      S = Str;
      n = S.size();
      SA = suffix_array<Uint>(S);
      t.next("suffix array");
      sequence<Uint> LCPs = lcp(S, SA);
      LCP = sequence<Uint>(n, [&] (size_t i) {
	  return (i==0) ? 0 : LCPs[i-1];}); // shift by 1
      LCPs.clear();
      t.next("LCP");
      suffix_tree_from_SA_LCP();
    }

    void suffix_tree_from_SA_LCP() {
      n = SA.size();
      sequence<Uint> Parents = cartesian_tree(LCP);
      t.next("Cartesian Tree");

      // A cluster is a set of connected nodes in the CT with the same LCP
      // Each cluster needs to be converted into a single node.
      // The following marks the roots of each cluster and has Roots point to self
      sequence<Uint> Roots(n);
      sequence<bool> is_root(n, [&] (Uint i) {
	  bool is_root_i = (i == 0) || (LCP[Parents[i]] != LCP[i]);
	  Roots[i] = is_root_i ? i : Parents[i];
	  return is_root_i;
	});

      // this shortcuts each node to point to its root
      parallel_for(0, n, [&] (Uint i) {
	  if (!is_root[i]) {
	    Uint root = Roots[i];
	    while (!is_root[root]) root = Roots[root];
	    Roots[i] = root;
	  }});
      t.next("Cluster roots");

      sequence<Uint> new_labels;
      Uint num_internal;
      std::tie(new_labels, num_internal) = enumerate<Uint>(is_root);
      t.next("Relabel roots");

      // interleave potential internal nodes and leaves
      auto edges = delayed_seq<Pair>(2*n, [&] (size_t j) {
	  if (j & 1) { // is leaf (odd)
	    Uint i = j/2;
	    // for each element of SA find larger LCP on either side
	    // the root of it will be the parent in the suffix tree
	    Uint parent = (i==n-1) ? i : (LCP[i] > LCP[i+1] ? i : i+1);
	    return Pair(new_labels[Roots[parent]], j);
	  } else { // is internal node (even)
	    Uint i = j/2;
	    Uint parent = Parents[i];
	    return Pair(new_labels[Roots[parent]], 2 * new_labels[i]);
	  }
	});

      // keep if leaf or root internal node (except the overall root).
      sequence<bool> flags(2*n, [&] (size_t j) {
	  return (j & 1)  || (is_root[j/2] && (j != 0));});

      sequence<Pair> all_edges = pack(edges, flags);
      t.next("generate edges");
      flags.clear();
      Roots.clear();
      new_labels.clear();

      auto get_first = [&] (Pair p) {return p.first;};

      sequence<Pair> sorted_edges =
	integer_sort(all_edges, get_first, log2_up(num_internal));
      t.next("Sort edges");

      sequence<Uint> offsets = get_counts<Uint>(sorted_edges, get_first, num_internal);
      //sequence<size_t> h = histogram<size_t>(offsets, 200);
      //int i = 2;
      //int total = 0;
      //while (h[i] > 0) std::cout << i  << " : " << h[i] << ", " << (total += h[i++]) << std::endl;

      scan_inplace(offsets.slice(), addm<Uint>());
      t.next("Get Counts");

      if (verbose)
	cout << "leaves = " << n << " internal nodes = " << num_internal << std::endl;

      // tag edges with character from s
      sequence<Uint> root_indices = pack_index<Uint>(is_root);
      Edges = map<edge>(sorted_edges, [&] (Pair p) {
	  Uint child = p.second;
	  Uint i = child/2;
	  bool is_leaf = child & 1;
	  int depth = (is_leaf ?
		       ((i==n-1) ? LCP[i] : std::max(LCP[i], LCP[i+1])) :
		       LCP[Parents[root_indices[i]]]);
	  int start = is_leaf ? SA[i] : SA[root_indices[i]];
	  //cout << i << ", " << depth << ", " << is_leaf << ", " << SA[i] << std::endl;
	  edge e = {i, S[start+depth], is_leaf ? leaf : internal};
	  return e;
	});
      t.next("project edges");
      sorted_edges.clear();
      Parents.clear();

      Nodes = sequence<node>(num_internal, [&] (size_t i) {
	  Uint lcp = LCP[root_indices[i]];
	  Uint offset = offsets[i];
	  Uint location = SA[root_indices[i]];
	  node r = {lcp, location, offset};
	  //cout << i << ", " << offset << ", " << location << ", " << lcp << std::endl;
	  return r;
	});
      t.next("Make nodes");
    }

    range<edge*> get_children(size_t i) {
      if (i == Nodes.size()-1)
	return Edges.slice(Nodes[i].offset, Edges.size());
      else return Edges.slice(Nodes[i].offset, Nodes[i+1].offset);
    }

    edge find_child(Uint i, uchar c) {
      for (edge e : get_children(i))
    	if (e.c == c) return e;
      return empty_edge;
    }

    std::optional<Uint> find(char const *s) {
      Uint node = 0;
      Uint j = 0;
      std::optional<Uint> None;
      while (true) {
	if (s[j] == 0) return std::optional<Uint>(Nodes[node].location);
	//cout << "j = " << j << " node = " << node << std::endl;
	edge e = find_child(node, s[j++]);
	switch (e.type) {
	case empty :
	  return None;
	case leaf :
	  //cout << "leaf" << std::endl;
	  while (true) {
	    if (s[j] == 0) return std::optional<Uint>(SA[e.child]);
	    if (s[j] != S[SA[e.child] + j]) return None;
	    j++;
	  }
	case internal :
	  //cout << "internal: " << Nodes[node].location << std::endl;
	  node = e.child;
	  size_t l = Nodes[node].lcp;
	  while (j < l) {
	    if (s[j] == 0) return std::optional<Uint>(Nodes[node].location);
	    if (s[j] != S[Nodes[node].location + j]) return None;
	    j++;
	  }
	}
      }
    }
  };

}  // namespace pbbs
