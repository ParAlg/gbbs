#include "sequence.h"
#include "random.h"
#include "get_time.h"

using namespace std;
using namespace pbbs;

// Maximum Contiguous Subsequence Sum

// 10x improvement by using delayed sequence
template <class Seq>
typename Seq::T mcss(Seq const &A) {
  using T = typename Seq::T;
  using TT = tuple<T,T,T,T>;
  auto f = [&] (TT a, TT b) {
    T aa, pa, sa, ta, ab, pb, sb, tb;
    tie(aa, pa, sa, ta) = a;
    tie(ab, pb, sb, tb) = b;
    return TT(max(aa,max(ab,sa+pb)),
	      max(pa, ta+pb),
	      max(sa + ab, sb),
	      ta + tb);
  };
  //auto S = sequence<TT>(A.size(), [&] (size_t i) {
  auto S = delayed_seq<TT>(A.size(), [&] (size_t i) {
      return TT(A[i],A[i],A[i],A[i]);});
  TT r = reduce(S, make_monoid(f, TT(0,0,0,0)));
  return get<0>(r);
}

int main() {
  using T = double;
  size_t n = 100000000;
  pbbs::random r(0);
  sequence<T> A(n, [&] (size_t i) {return (T) (r[i]%n - n/2);});
  T result;
  for (int i=0; i < 5; i++) {
    result = mcss(A);
    t.next("Total");
  }
  cout << result << endl;
}
