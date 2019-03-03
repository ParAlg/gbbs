#include "sample_sort.h"

int main (int argc, char *argv[]) {
  //using T = size_t;
  //using T = __int128;
  //using T = double;
  using E = double;
  using T = std::pair<E,E>;
  timer t;
  size_t n = 25000;
  size_t blocks = 100000000/n;
  size_t rounds = 10;
  size_t total = n * blocks;
  auto less = [&] (T a, T b) {return a.first <  b.first;};
  //auto less = std::less<T>();
  //sequence<T> in(total, [&] (size_t i) {return hash64(i)%total;});
  sequence<T> in(total, [&] (size_t i) {return T(hash64(i)%total,0.0);});
  sequence<T> out = in;
  for (size_t j=0; j<rounds; j++) {
    out = in;
    t.start();
    parallel_for (0, blocks, [&] (size_t i) {
	range<T*> a = out.slice(i*n,(i+1)*n);
	//sample_sort_inplace(a, less);
	bucket_sort(a, less, false);
	//quicksort(a.begin(), a.size(), less);
      }, 1);
    double tm = t.stop();
    cout << "time = " << tm << endl;
    // parallel_for (0, blocks, [&] (size_t i) {
    // 	parallel_for (0, n-1, [&] (size_t j) {
    // 	    if (less(out[i*n+j+1],out[i*n+j])) abort();
    // 	  });});
  }
  //for (size_t i = 0; i < 10; i++)
  //  cout << in[i] << ", " << out[i] << endl;
}
  
