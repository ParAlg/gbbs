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

#include "sequence.h"
#include "get_time.h"
#include "suffix_array.h"
#include "integer_sort.h"
#include "lcp.h"
#include "cartesian_tree.h"

namespace pbbs {


  template <class Uint>
  struct suffix_tree {
    using edge = Uint;
    using Pair = std::pair<Uint,Uint>;

    struct node {
      Uint lcp;
      Uint location;
      Uint offset;
    };

    struct ct_node {
      Uint parent;
      Uint value;};
    
    sequence<edge> Edges;
    sequence<node> Nodes;
    sequence<ct_node> CT_Nodes;
    Uint n;
    Uint num_edges;
    
    void set_cluster_root(size_t i) {
      auto root = CT_Nodes[i].parent;
      while (root != 0 &&
	     CT_Nodes[CT_Nodes[root].parent].value == CT_Nodes[root].value)
	root = CT_Nodes[root].parent;
      CT_Nodes[i].parent = root;
    }

    suffix_tree(sequence<unsigned char> const &S) {
      timer t("suffix tree", true);
      n = S.size();
      cout << n << endl;
      sequence<Uint> SA = suffix_array<Uint>(S);
      t.next("suffix array");
      sequence<Uint> LCPs_ = LCP(S, SA);
      sequence<Uint> LCPs(n, [&] (size_t i) {
	  return (i==0) ? 0 : LCPs_[i-1];}); // shift by 1
      t.next("LCP");
      LCPs_.clear();

      sequence<Uint> Parents = cartesian_tree(LCPs);
      t.next("Cartesian Tree");
      
      //for (size_t i=0; i < n; i++)
      //cout << i << ", " << LCPs[i] << ", " << Parents[i] << endl;

      // A cluster is a set of connected nodes in the CT with the same LCP
      // Each cluster needs to be converted into a single node.
      // The following points each node to the root of its cluster
      sequence<Uint> Roots = Parents;
      sequence<bool> is_root(n);
      parallel_for(0, n, [&] (Uint i) {
	  Uint root = i;
	  while (root != 0 && root != Roots[root] && LCPs[Roots[root]] == LCPs[root]) 
	    root = Roots[root];
	  is_root[i] = (root == i);
	  Roots[i] = root;
	}, 1000);
      t.next("Cluster roots");
      is_root[0] = false;
      
      // cout << "next " << endl;
      //for (size_t i=0; i < n; i++)
      //cout << i << ", " << Parents[i] << ", " << Roots[i] << ", " << is_root[i] << ", " << LCPs[i] << endl;

      auto ie = delayed_seq<Pair>(n, [&] (size_t i) {
	  return Pair(Roots[Parents[i]], 2*i);});

      // edges from internal nodes to their parents
      auto internal_edges = pack(ie, is_root);
      //for (size_t i=0; i < idxs.size(); i++)
      //cout << internal_edges[i].first << ", " << internal_edges[i].second << endl;
      t.next("Internal edges");
      Parents.clear();
      is_root.clear();

      // edges from leaf nodes to their parents
      sequence<Pair> leaf_edges(n, [&] (size_t i) {
	  // for each element of SA find larger LCP on either side
	  Uint parent = (i==n-1) ? i : (LCPs[i] > LCPs[i+1] ? i : i+1);
	  // get cluster root
	  parent = Roots[parent];
	  return Pair(parent, 2*i + 1);
	});
      t.next("Leaf edges");
      Roots.clear();
      
      //for (size_t i=0; i < n; i++)
      //cout << leaf_edges[i].first << ", " << leaf_edges[i].second << endl;

      // sort all edges by parent
      auto sc = integer_sort_with_counts<Uint>(append(internal_edges, leaf_edges),
					       [&] (Pair p) {return p.first;},
					       n);
      internal_edges.clear();
      leaf_edges.clear();
      t.next("Sort edges");
      sequence<Pair> sorted_edges = std::move(sc.first);
      num_edges = sorted_edges.size();
      sequence<Uint> offsets = std::move(sc.second);
      scan_inplace(offsets.slice(), addm<Uint>());
      t.next("Scan edges");
      cout << "n = " << n << " m = " << num_edges << endl;
      
      //for (size_t i=0; i < num_edges; i++)
      //  cout << sorted_edges[i].first << ", " << sorted_edges[i].second << endl;
      
      Edges = sequence<edge>(sorted_edges.size(), [&] (size_t i) {
	  return sorted_edges[i].second;});
      t.next("project edges");
      sorted_edges.clear();
      
      Nodes = sequence<node>(n, [&] (size_t i) {
	  Uint lcp = LCPs[i];
	  Uint offset = offsets[i];
	  Uint location = SA[i];
	  node r = {lcp, location, offset};
	  //cout << i << ", " << offset << ", " << location << ", " << lcp << endl;
	  return r;
	});
      t.next("Make nodes");
    }
  };
}
