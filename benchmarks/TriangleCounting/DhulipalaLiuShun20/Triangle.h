// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

#include <algorithm>
#include <cmath>
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "gbbs/gbbs.h"
// #include "gbbs/pbbslib/sparse_table.h"
// #include "two_level_tables.h"


namespace gbbs {
using namespace std;

#define OLD_EDGE 0
#define NEW_INSERTION 1
#define NEW_DELETION 2

// // inline std::tuple<std::pair<uintE, uintE>, int8_t> newKV(uintE s, uintE d, int8_t v){
// //   return make_tuple(make_pair(s,d), v);
// // }

// // inline std::tuple<std::pair<uintE, uintE>, size_t> newKV(uintE s, uintE d, size_t v){
// //   return make_tuple(make_pair(s,d), v);
// // }

// // struct hash_pair {
// //   inline size_t operator () (const std::tuple<uintE, uintE>& t) {
// //     size_t l = std::min(std::get<0>(t), std::get<1>(t));
// //     size_t r = std::max(std::get<0>(t), std::get<1>(t));
// //     size_t key = (l << 32) + r;
// //     return pbbslib::hash64_2(key);
// //   }
// // };

// struct hash_vertex {
//   inline size_t operator () (const uintE t) {
//     return pbbslib::hash64_2(t);
//   }
// };


// template <class Graph, class E, class F>
// struct BPDTriangleCountState {
//   using K = uintE;
//   using BV = int;
//   using BT = std::tuple<K, BV>;
//   using V = sparse_table<K, BV, hash_vertex>; // top value
//   // using T = std::tuple<K, V>; // top key value pair

//   // initialize tables assuming state.D is already initialized
//   // use D to determine the low/high of vertices
//   // update HH, HL, LH, LL
//   struct updateTablesF {
//     updateTablesF() {}

//     inline bool update(uintE s, uintE d) {
//     // can add condition s < d and add both edges in one call
//       add_to_tables(s,d,OLD_EDGE);
//       return 1;
//     }// when to return false?

//     inline bool updateAtomic(uintE s, uintE d) {
//       update(s,d);
//       return 1;
//     }

//     inline bool cond(uintE d) { return cond_true(d); }
//   };

//   // update T
//   struct updateTF {
//     void operator ()(size_t* v0, std::tuple<K, size_t>& kv){
//       pbbslib::write_add(v0, std::get<1>(kv));
//     }

//   };

//   Graph& G;
//   size_t M;
//   size_t t1;
//   size_t t2;
//   pbbs::sequence<uintE> D;
//   E& HH; 
//   E& HL; 
//   E& LH; 
//   E& LL; 
//   F& T;

//   BPDTriangleCountState(Graph& _G, E _HH, E _HL, E _LH, E _LL, F _T, pbbs::sequence<uintE> _D): 
//     G(_G),
//     // D(_D),
//     HH(_HH),
//     HL(_HL),
//     LH(_LH),
//     LL(_LL),
//     T(_T),
//     D(_D){

//     size_t n = G.n;
//     M  = 2 * G.m + 1;
//     t1 = sqrt(M) / 2;
//     t2 = 3 * t1;


//     // auto init_D_f = [&](const uintE& v) {
//     //   D[v] = G.get_vertex(v).getOutDegree();
//     // };
//     // auto frontier = pbbs::sequence<bool>(n, true);
//     // vertexSubset Frontier(n,n,frontier.to_array());
//     // vertexMap(activeAndCts, init_D_f);

    
//     // edgeMap(G, Frontier, updateTablesF());
//     for(T kv : HL.entries()){ // loop over keys?
//       K vhigh = std::get<0>(kv);
//       BV table = std::get<1>(kv);
//       for(BT kv2 : table.entries()){
//         K vlow  = std::get<0>(kv2);
//         V table2 = LH.find(vlow, LH.empty_val);

//         for(BT kv3 : table2.entries()){
//           K vhigh2 = std::get<0>(kv);
//           add_to_T(vhigh, vhigh2, 1);
//           add_to_T(vhigh2, vhigh, 1);
//         }

//       }
//     }
//   }


//   ~BPDTriangleCountState(){
//     // free(D); ??
//     // HH.del();
//     // HL.del();
//     // LH.del();
//     // LL.del();
//     // T.del();
//   }

//   inline bool is_high(uintE s){
//     return D[s] >= t2;
//   }

//   inline void add_to_tables(uintE s, uintE d, int8_t V){
//     bool highS = is_high(s);
//     bool highD = is_high(d);
//     if( highS && highD){ //HH
//       HH.insert(s,d,V);
//     }else if(highS){ //HL
//       HL.insert(s,d,V);
//     }else if(highD){ //LH
//       LH.insert(s,d,V);
//     }else{ //LL
//       LL.insert(s,d,V);
//     }
//   }

//   inline void add_to_T(uintE s, uintE d, size_t V){
//     T.insert_f(s,d,V, updateTF());
//   }



// };




// template <class Graph>
// inline auto Initialize(Graph& G){
//   using K = uintE;
//   using V = int;
//   using BT = std::tuple<K, V>;
//   // using W = typename Graph::weight_type;
//   size_t n = G.n;
//   BT empty = std::make_tuple(UINT_E_MAX, -1);

//   pbbs::sequence<uintE> D = pbbs::sequence<uintE>(n, [&] (size_t i) { return G.get_vertex(i).getOutDegree(); });
//   auto HH = NestHash::array_table<K, V, hash_vertex>(n, empty, hash_vertex(), D.begin());
//   auto HL = NestHash::array_table<K, V, hash_vertex>(n, empty, hash_vertex(), D.begin());
//   auto LH = NestHash::array_table<K, V, hash_vertex>(n, empty, hash_vertex(), D.begin());
//   auto LL = NestHash::array_table<K, V, hash_vertex>(n, empty, hash_vertex(), D.begin());
//   auto T  = NestHash::array_table<K, size_t, hash_vertex>(2*n, empty, hash_vertex(), D.begin());
//   // what to initialize?

//   using E = decltype(HH);
//   using F = decltype(T);
  
//   auto state = BPDTriangleCountState<Graph, E, F>(G, HH, HL, LH, LL, T, D);

//   return state;
// }

// template <class Graph>
// inline auto test(Graph& G){
//   using K = uintE;
//   using V = int;
//   using BT = std::tuple<K, V>;
//   // using W = typename Graph::weight_type;
//   size_t n = G.n;
//   BT empty = std::make_tuple(UINT_E_MAX, -1);

//   auto HH = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   cout << "constructed" << endl;
//   HH.insert(1, 2, 3);
//   cout << "inserted" << endl;

//   auto low = HH.find(1, HH.empty_val);
//   cout << "found table" << endl;

//   low.insert(make_tuple(3, 4));
//   V v = HH.find(1,3,0);
//   cout << v << endl;


// }


template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, commandLine& P) {
  // auto C0 = P.getOptionIntValue("-c", 0);
  // test(G);

  // auto state = Initialize<Graph>(G);

  return 0;
}


// template <class Graph>
// inline auto Initialize(Graph& G){
//   using K = uintE;
//   using V = int;
//   using BT = std::tuple<K, V>;
//   // using W = typename Graph::weight_type;
//   size_t n = G.n;
//   BT empty = std::make_tuple(UINT_E_MAX, -1);

//   auto HH = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto HL = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto LH = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto LL = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto T  = NestHash::nested_table<K, size_t, hash_vertex>(2*n, empty, hash_vertex());
//   // auto D  = sequence<size_t>(n);

//   using E = decltype(HH);
//   using F = decltype(T);
  
//   auto state = BPDTriangleCountState<Graph, E, F>(G, HH, HL, LH, LL, T);

//   return state;
// }

}  // namespace gbbs
