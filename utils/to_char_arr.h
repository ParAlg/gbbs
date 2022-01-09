#pragma once

namespace gbbs {
inline int xToStringLen(long a) { return 21; }
inline void xToString(char* s, long a) { sprintf(s, "%ld", a); }

inline int xToStringLen(unsigned long a) { return 21; }
inline void xToString(char* s, unsigned long a) { sprintf(s, "%lu", a); }

inline uint xToStringLen(uint a) { return 12; }
inline void xToString(char* s, uint a) { sprintf(s, "%u", a); }

inline int xToStringLen(int a) { return 12; }
inline void xToString(char* s, int a) { sprintf(s, "%d", a); }

inline int xToStringLen(double a) { return 18; }
inline void xToString(char* s, double a) { sprintf(s, "%.11le", a); }

inline int xToStringLen(char* a) { return strlen(a) + 1; }
inline void xToString(char* s, char* a) { sprintf(s, "%s", a); }

template <class A, class B>
inline int xToStringLen(std::pair<A, B> a) {
  return xToStringLen(a.first) + xToStringLen(a.second) + 1;
}
template <class A, class B>
inline void xToString(char* s, std::pair<A, B> a) {
  int l = xToStringLen(a.first);
  xToString(s, a.first);
  s[l] = ' ';
  xToString(s + l + 1, a.second);
}

inline int xToStringLen(gbbs::empty a) { return 0; }
inline void xToString(char* s, gbbs::empty a) {}

template <class A, class B, class C>
inline int xToStringLen(std::tuple<A, B, C> a) {
  return xToStringLen(std::get<0>(a)) + xToStringLen(std::get<1>(a)) +
         xToStringLen(std::get<2>(a)) + 1;
}
template <class A, class B, class C>
inline void xToString(char* s, std::tuple<A, B, C> a) {
  int l = xToStringLen(std::get<0>(a));
  xToString(s, std::get<0>(a));
  s[l] = ' ';
  int r = xToStringLen(std::get<1>(a));
  xToString(s + l + 1, std::get<1>(a));
  int e = xToStringLen(std::get<2>(a));
  if (e > 0) {  // empty weight
    s[l + r + 1] = ' ';
    xToString(s + l + r + 2, std::get<2>(a));
  }
}

template <class S>
sequence<char> arrayToString(S& A) {
  auto AA = A.begin();
  size_t n = A.size();
  auto L = sequence<size_t>::from_function(
      n, [&](size_t i) { return xToStringLen(A[i]) + 1; });
  size_t m = parlay::scan_inplace(make_slice(L));
  auto out = sequence<char>(m, (char)0);
  parallel_for(0, n - 1, [&](size_t i) {
    xToString(out.begin() + L[i], AA[i]);
    out[L[i + 1] - 1] = '\n';
  });
  xToString(out.begin() + L[n - 1], A[n - 1]);
  out[m - 1] = '\n';

  return filter(out, [&](char c) { return c != (char)0; });
}

template <class S>
void writeArrayToStream(std::ofstream& os, S& arr) {
  size_t BSIZE = 1000000;
  size_t offset = 0;
  size_t n = arr.size();
  while (offset < n) {
    // Generates a string for a sequence of size at most BSIZE
    // and then wrties it to the output stream
    auto arr_slice = arr.cut(offset, offset + std::min(BSIZE, n - offset));
    auto Q = arrayToString(arr_slice);
    os.write(Q.begin(), Q.size());
    offset += BSIZE;
  }
}
}  // namespace gbbs
