#include "ligra.h"

struct BF_F {
  intE* SP;
  intE* Visited;
  BF_F(intE* _SP, intE* _Visited) : SP(_SP), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d, const intE& edgeLen) {
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
                           const intE& edgeLen) {
    intE newDist = SP[s] + edgeLen;
    return (writeMin(&SP[d], newDist) && CAS(&Visited[d], 0, 1));
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

// reset visited vertices
struct BF_Vertex_F {
  intE* Visited;
  BF_Vertex_F(intE* _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 0;
    return 1;
  }
};

template <template <class W> class vertex, class W>
auto BellmanFord(graph<vertex<W>>& GA, const uintE& start) {
  size_t n = GA.n;
  auto Visited = array_imap<int>(n, 0);
  auto SP = array_imap<intE>(n, INT_MAX / 2);
  SP[start] = 0;

  vertexSubset Frontier(n, start);
  size_t round = 0;
  while (!Frontier.isEmpty()) {
    // Check for a negative weight cycle
    if (round == n) {
      parallel_for(long i = 0; i < n; i++) SP[i] = -(INT_E_MAX / 2);
      break;
    }
    auto em_f =
        wrap_with_default<W, intE>(BF_F(SP.start(), Visited.start()), (intE)1);
    auto output =
        edgeMap(GA, Frontier, em_f, GA.m / 10, sparse_blocked | dense_forward);
    vertexMap(output, BF_Vertex_F(Visited.start()));
    cout << output.size() << endl;
    Frontier.del();
    Frontier = output;
    round++;
  }
  auto dist_im = make_in_imap<size_t>(
      n, [&](size_t i) { return (SP[i] == (INT_MAX / 2)) ? 0 : SP[i]; });
  cout << "max dist = " << pbbs::reduce_max(dist_im) << endl;
  cout << "n rounds = " << round << endl;
  return SP;
}
