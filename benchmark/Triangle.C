// Usage:
// numactl -i all ./Triangle -rounds 2 -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

#include "Triangle.h"

template <class vertex>
void Triangle_runner(graph<vertex>& GA, commandLine P) {
  assert(P.getOption("-s"));
  size_t count = 0;
  count = Triangle(GA);
  cout << "triangle count = " << count << endl;
  if (P.getOption("-stats")) {
    auto wedge_im = make_in_imap<size_t>(GA.n, [&] (size_t i) {
      size_t deg = GA.V[i].getOutDegree();
      return (deg*deg-1)/2;
    });
    size_t n_wedges = pbbs::reduce_add(wedge_im);
    cout << "n_wedges = " << n_wedges << endl;
    cout << "triangle density = " << ((3.0*count)/n_wedges) << endl;
  }
}

generate_main(Triangle_runner, false);
