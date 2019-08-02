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

#include "ligra.h"

template <class W, class GW>
struct BF_F {
  intE* SP;
  intE* Visited;
  GW& get_weight;
  BF_F(intE* _SP, intE* _Visited, GW& get_weight) : SP(_SP), Visited(_Visited), get_weight(get_weight) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    intE edgeLen = get_weight(s, d, wgh);
    intE newDist = SP[s] + edgeLen;
    if (SP[d] > newDist) {
      SP[d] = newDist;
      if (Visited[d] == 0) {
        Visited[d] = 1;
        return 1;
      }
    }
    return 0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d,
                           const W& wgh) {
    intE edgeLen = get_weight(s, d, wgh);
    intE newDist = SP[s] + edgeLen;
    return (pbbslib::write_min(&SP[d], newDist) && pbbslib::atomic_compare_and_swap(&Visited[d], 0, 1));
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

template <class W, class GW>
BF_F<W, GW> make_bf_f(intE* SP, intE* Visited, GW& get_weight) {
  return BF_F<W, GW>(SP, Visited, get_weight);
}

// reset visited vertices
struct BF_Vertex_F {
  intE* Visited;
  BF_Vertex_F(intE* _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 0;
    return 1;
  }
};

template <class G>
inline sequence<intE> BellmanFord(G& GA, const uintE& start) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  auto Visited = sequence<int>(n, 0);
  auto SP = sequence<intE>(n, INT_MAX / 2);
  SP[start] = 0;

  auto degree_f = [&](size_t i) {
    auto vtx = GA.get_vertex(i);
    return std::max(vtx.getInDegree(), vtx.getOutDegree());
  };
  auto degree_im = pbbslib::make_sequence<size_t>(n, degree_f);
  size_t max_degree = pbbslib::reduce_max(degree_im);
  size_t normalize = 2*max_degree+1;

  // edge -> custom_weight
  auto get_weight = [&] (const uintE& u, const uintE& v, const W& wgh) -> intE {
    return wgh;
//    uintE deg_u = GA.get_vertex(u).getOutDegree();
//    uintE deg_v = GA.get_vertex(v).getOutDegree();
//    return pbbs::log2_up((size_t)((1/static_cast<double>(deg_u + deg_v + 1))*normalize));
  };

  vertexSubset Frontier(n, start);
  size_t round = 0;
  while (!Frontier.isEmpty()) {
    // Check for a negative weight cycle
    if (round == n) {
      par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { SP[i] = -(INT_E_MAX / 2); });
      break;
    }
    auto em_f = make_bf_f<W>(SP.begin(), Visited.begin(), get_weight);
    auto output =
        edgeMap(GA, Frontier, em_f, GA.m / 10, no_dense);
    vertexMap(output, BF_Vertex_F(Visited.begin()));
    std::cout << output.size() << "\n";
    Frontier.del();
    Frontier = output;
    round++;
  }
  auto dist_im_f = [&](size_t i) { return (SP[i] == (INT_MAX / 2)) ? 0 : SP[i]; };
  auto dist_im = pbbslib::make_sequence<size_t>(n, dist_im_f);
  std::cout << "max dist = " << pbbslib::reduce_max(dist_im) << "\n";
  std::cout << "n rounds = " << round << "\n";
  return SP;
}
