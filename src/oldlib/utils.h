// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
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

#include <stdlib.h>
#include <fstream>
#include <iostream>

#include "pbbslib/utilities.h"

// scan/filter macros; used by sequence implementations
#define _SCAN_LOG_BSIZE 10
#define _SCAN_BSIZE (1 << _SCAN_LOG_BSIZE)
#define _F_BSIZE (2 * _SCAN_BSIZE)

//template <class ET>
//inline bool CAS(ET* ptr, ET oldv, ET newv) {
//  if (sizeof(ET) == 1) {
//    return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv),
//                                        *((bool*)&newv));
//  } else if (sizeof(ET) == 4) {
//    return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv),
//                                        *((int*)&newv));
//  } else if (sizeof(ET) == 8) {
//    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv),
//                                        *((long*)&newv));
//  } else {
//    std::cout << "CAS bad length : " << sizeof(ET) << "\n";
//    abort();
//  }
//}

template <typename ET>
inline bool CAS(ET* ptr, const ET oldv, const ET newv) {
  return __sync_bool_compare_and_swap(ptr, oldv, newv);
}

//template <class ET>
//inline bool writeMin(ET* a, ET b) {
//  ET c;
//  bool r = 0;
//  do
//    c = *a;
//  while (c > b && !(r = CAS(a, c, b)));
//  return r;
//}

//template <class ET, class F>
//inline bool writeMin(ET* a, ET b, F less) {
//  ET c;
//  bool r = 0;
//  do
//    c = *a;
//  while (less(b, c) && !(r = CAS(a, c, b)));
//  return r;
//}
//
//template <class ET>
//inline bool writeMax(ET* a, ET b) {
//  ET c;
//  bool r = 0;
//  do
//    c = *a;
//  while (c < b && !(r = CAS(a, c, b)));
//  return r;
//}

// template <class ET>
// inline ET writeAdd(ET *a, ET b) {
//  return b + pbbslib::xadd<ET>(a, b);
//}

//template <class ET>
//inline ET writeAdd(ET* a, ET b) {
//  volatile ET newV, oldV;
//  do {
//    oldV = *a;
//    newV = oldV + b;
//  } while (!pbbslib::atomic_compare_and_swap(a, oldV, newV));
//  return newV;
//}

namespace ligra_utils {
template <class E>
struct identityF {
  E operator()(const E& x) { return x; }
};

template <class E>
struct addF {
  E operator()(const E& a, const E& b) const { return a + b; }
};

template <class E>
struct minF {
  E operator()(const E& a, const E& b) const { return (a < b) ? a : b; }
};

template <class E>
struct maxF {
  E operator()(const E& a, const E& b) const { return (a > b) ? a : b; }
};

struct nonMaxF {
  bool operator()(uintE& a) { return (a != UINT_E_MAX); }
};

template <class E, class I, class P>
struct filter_iter {
  I& iter;
  P& pred;
  E cur_val;

  filter_iter(I& _it, P& _pr) : iter(_it), pred(_pr) {
    cur_val = iter.cur();
    while (!pred(cur_val) && iter.has_next()) {
      cur_val = iter.next();
    }
  }

  E cur() { return cur_val; }

  E next() {
    while (iter.has_next()) {
      cur_val = iter.next();
      if (pred(cur_val)) {
        break;
      }
    }
    return cur_val;
  }

  // has_next
};

template <class E, class I, class P>
inline filter_iter<E, I, P> make_filter_iter(I& _it, P& _pr) {
  return filter_iter<E, I, P>(_it, _pr);
}

// Sugar to pass in a single f and get a struct suitable for edgeMap.
template <class W, class F>
struct EdgeMap_F {
  F f;
  EdgeMap_F(F _f) : f(_f) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    return f(s, d, wgh);
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    return f(s, d, wgh);
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class W, class F>
inline EdgeMap_F<W, F> make_em_f(F f) {
  return EdgeMap_F<W, F>(f);
}

template <class T>
struct _seq {
  T* A;
  long n;
  _seq() {
    A = NULL;
    n = 0;
  }
  _seq(T* _A, long _n) : A(_A), n(_n) {}
  void clear() { pbbslib::free_array(A); }
};

namespace seq {
template <class intT>
struct boolGetA {
  bool* A;
  boolGetA(bool* AA) : A(AA) {}
  intT operator()(intT i) { return (intT)A[i]; }
};

template <class ET, class intT>
struct getA {
  ET* A;
  getA(ET* AA) : A(AA) {}
  ET operator()(intT i) { return A[i]; }
};

template <class IT, class OT, class intT, class F>
struct getAF {
  IT* A;
  F f;
  getAF(IT* AA, F ff) : A(AA), f(ff) {}
  OT operator()(intT i) { return f(A[i]); }
};

#define nblocks(_n, _bsize) (1 + ((_n)-1) / (_bsize))

#define blocked_for(_i, _s, _e, _bsize, _body) \
  {                                            \
    intT _ss = _s;                             \
    intT _ee = _e;                             \
    intT _n = _ee - _ss;                       \
    intT _l = nblocks(_n, _bsize);             \
    par_for(0, _l, 1, [&] (size_t _i) {     \
      intT _s = _ss + _i * (_bsize);           \
      intT _e = std::min(_s + (_bsize), _ee);  \
      _body                                    \
    });                                        \
  }

template <class OT, class intT, class F, class G>
inline OT reduceSerial(intT s, intT e, F f, G g) {
  OT r = g(s);
  for (intT j = s + 1; j < e; j++) r = f(r, g(j));
  return r;
}

template <class OT, class intT, class F, class G>
inline OT reduce(intT s, intT e, F f, G g) {
  intT l = nblocks(e - s, _SCAN_BSIZE);
  if (l <= 1) return reduceSerial<OT>(s, e, f, g);
  OT* Sums = pbbslib::new_array_no_init<OT>(l);
  blocked_for(i, s, e, _SCAN_BSIZE, Sums[i] = reduceSerial<OT>(s, e, f, g););
  OT r = reduce<OT>((intT)0, l, f, getA<OT, intT>(Sums));
  pbbslib::free_array(Sums);
  return r;
}

template <class OT, class intT, class F>
inline OT reduce(OT* A, intT n, F f) {
  return reduce<OT>((intT)0, n, f, getA<OT, intT>(A));
}

template <class OT, class intT>
inline OT plusReduce(OT* A, intT n) {
  return reduce<OT>((intT)0, n, addF<OT>(), getA<OT, intT>(A));
}

// g is the map function (applied to each element)
// f is the reduce function
// need to specify OT since it is not an argument
template <class OT, class IT, class intT, class F, class G>
inline OT mapReduce(IT* A, intT n, F f, G g) {
  return reduce<OT>((intT)0, n, f, getAF<IT, OT, intT, G>(A, g));
}

template <class intT>
inline intT sum(bool* In, intT n) {
  return reduce<intT>((intT)0, n, addF<intT>(), boolGetA<intT>(In));
}

template <class ET, class intT, class F, class G>
inline ET scanSerial(ET* Out, intT s, intT e, F f, G g, ET zero, bool inclusive,
                     bool back) {
  ET r = zero;
  if (inclusive) {
    if (back)
      for (intT i = e - 1; i >= s; i--) Out[i] = r = f(r, g(i));
    else
      for (intT i = s; i < e; i++) Out[i] = r = f(r, g(i));
  } else {
    if (back)
      for (intT i = e - 1; i >= s; i--) {
        ET t = g(i);
        Out[i] = r;
        r = f(r, t);
      }
    else
      for (intT i = s; i < e; i++) {
        ET t = g(i);
        Out[i] = r;
        r = f(r, t);
      }
  }
  return r;
}

template <class ET, class intT, class F>
inline ET scanSerial(ET* In, ET* Out, intT n, F f, ET zero) {
  return scanSerial(Out, (intT)0, n, f, getA<ET, intT>(In), zero, false, false);
}

// back indicates it runs in reverse direction
template <class ET, class intT, class F, class G>
inline ET scan(ET* Out, intT s, intT e, F f, G g, ET zero, bool inclusive,
               bool back) {
  intT n = e - s;
  intT l = nblocks(n, _SCAN_BSIZE);
  if (l <= 2) return scanSerial(Out, s, e, f, g, zero, inclusive, back);
  ET* Sums = pbbslib::new_array_no_init<ET>(nblocks(n, _SCAN_BSIZE));
  blocked_for(i, s, e, _SCAN_BSIZE, Sums[i] = reduceSerial<ET>(s, e, f, g););
  ET total = scan(Sums, (intT)0, l, f, getA<ET, intT>(Sums), zero, false, back);
  blocked_for(i, s, e, _SCAN_BSIZE,
              scanSerial(Out, s, e, f, g, Sums[i], inclusive, back););
  pbbslib::free_array(Sums);
  return total;
}

template <class ET, class intT, class F>
inline ET scan(ET* In, ET* Out, intT n, F f, ET zero) {
  return scan(Out, (intT)0, n, f, getA<ET, intT>(In), zero, false, false);
}

template <class ET, class intT, class F>
inline ET scanI(ET* In, ET* Out, intT n, F f, ET zero) {
  return scan(Out, (intT)0, n, f, getA<ET, intT>(In), zero, true, false);
}

template <class ET, class intT, class F>
inline ET scanBack(ET* In, ET* Out, intT n, F f, ET zero) {
  return scan(Out, (intT)0, n, f, getA<ET, intT>(In), zero, false, true);
}

template <class ET, class intT, class F>
inline ET scanIBack(ET* In, ET* Out, intT n, F f, ET zero) {
  return scan(Out, (intT)0, n, f, getA<ET, intT>(In), zero, true, true);
}

template <class ET, class intT>
inline ET plusScan(ET* In, ET* Out, intT n) {
  return scan(Out, (intT)0, n, addF<ET>(), getA<ET, intT>(In), (ET)0, false,
              false);
}

// sums a seq of n boolean flags
// an optimized version that sums blocks of 4 booleans by treating
// them as an integer
// Only optimized when n is a multiple of 512 and Fl is 4byte aligned
template <class intT>
inline intT sumFlagsSerial(bool* Fl, intT n) {
  intT r = 0;
  if (n >= 128 && (n & 511) == 0 && ((long)Fl & 3) == 0) {
    int* IFl = (int*)Fl;
    for (size_t k = 0; k < (n >> 9); k++) {
      size_t rr = 0;
      for (size_t j = 0; j < 128; j++) rr += IFl[j];
      r += (rr & 255) + ((rr >> 8) & 255) + ((rr >> 16) & 255) +
           ((rr >> 24) & 255);
      IFl += 128;
    }
  } else
    for (intT j = 0; j < n; j++) r += Fl[j];
  return r;
}

template <class ET, class intT, class F>
inline _seq<ET> packSerial(ET* Out, bool* Fl, intT s, intT e, F f) {
  if (Out == NULL) {
    intT m = sumFlagsSerial(Fl + s, e - s);
    Out = pbbslib::new_array_no_init<ET>(m);
  }
  intT k = 0;
  for (intT i = s; i < e; i++)
    if (Fl[i]) Out[k++] = f(i);
  return _seq<ET>(Out, k);
}

template <class ET, class intT, class F>
inline _seq<ET> pack(ET* Out, bool* Fl, intT s, intT e, F f) {
  intT l = nblocks(e - s, _F_BSIZE);
  if (l <= 1) return packSerial(Out, Fl, s, e, f);
  intT* Sums = pbbslib::new_array_no_init<intT>(l);
  blocked_for(i, s, e, _F_BSIZE, Sums[i] = sumFlagsSerial(Fl + s, e - s););
  intT m = plusScan(Sums, Sums, l);
  if (Out == NULL) Out = pbbslib::new_array_no_init<ET>(m);
  blocked_for(i, s, e, _F_BSIZE, packSerial(Out + Sums[i], Fl, s, e, f););
  pbbslib::free_array(Sums);
  return _seq<ET>(Out, m);
}

template <class ET, class intT>
inline intT pack(ET* In, ET* Out, bool* Fl, intT n) {
  return pack(Out, Fl, (intT)0, n, getA<ET, intT>(In)).n;
}

template <class intT>
inline _seq<intT> packIndex(bool* Fl, intT n) {
  return pack((intT*)NULL, Fl, (intT)0, n, identityF<intT>());
}

template <class intT>
inline intT packIndex(intT* Out, bool* Fl, intT n) {
  return pack(Out, Fl, (intT)0, n, identityF<intT>()).n;
}

template <class ET, class intT, class PRED>
inline intT filter(ET* In, ET* Out, bool* Fl, intT n, PRED p) {
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { Fl[i] = (bool)p(In[i]); });
  intT m = pack(In, Out, Fl, n);
  return m;
}

template <class ET, class intT, class PRED>
inline intT filter(ET* In, ET* Out, intT n, PRED p) {
  bool* Fl = pbbslib::new_array_no_init<bool>(n);
  intT m = filter(In, Out, Fl, n, p);
  pbbslib::free_array(Fl);
  return m;
}

}  // namespace seq

inline uint hashInt(uint a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

inline ulong hashInt(ulong a) {
  a = (a + 0x7ed55d166bef7a1d) + (a << 12);
  a = (a ^ 0xc761c23c510fa2dd) ^ (a >> 9);
  a = (a + 0x165667b183a9c0e1) + (a << 59);
  a = (a + 0xd3a2646cab3487e3) ^ (a << 49);
  a = (a + 0xfd7046c5ef9ab54c) + (a << 3);
  a = (a ^ 0xb55a4f090dd4a67b) ^ (a >> 32);
  return a;
}

// Remove duplicate integers in [0,...,n-1].
// Assumes that flags is already allocated and cleared to UINT_E_MAX.
// Sets all duplicate values in the array to UINT_E_MAX and resets flags to
// UINT_E_MAX.
template <class G>
inline void remDuplicates(G& get_key, uintE* flags, long m, long n) {
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintE key = get_key(i);
    if (key != UINT_E_MAX && flags[key] == UINT_E_MAX) {
      CAS(&flags[key], (uintE)UINT_E_MAX, static_cast<uintE>(i));
    }
  });
  // reset flags
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintE key = get_key(i);
    if (key != UINT_E_MAX) {
      if (flags[key] == i) {      // win
        flags[key] = UINT_E_MAX;  // reset
      } else {
        get_key(i) = UINT_E_MAX;  // lost
      }
    }
  });
}
}  // namespace ligra_utils
