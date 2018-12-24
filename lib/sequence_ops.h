// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2016 Guy Blelloch and the PBBS team
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

#include <iostream>

#include "macros.h"
#include "seq.h"
#include "utilities.h"

namespace pbbs {

constexpr const size_t _log_block_size = 12;
constexpr const size_t _block_size = (1 << _log_block_size);

inline size_t num_blocks(size_t n, size_t block_size) {
  return (1 + ((n)-1) / (block_size));
}

template <class F>
inline void sliced_for(size_t n, size_t block_size, const F& f) {
  size_t l = num_blocks(n, block_size);
  par_for(0, l, 1, [&] (size_t i) {
    size_t s = i * block_size;
    size_t e = std::min(s + block_size, n);
    f(i, s, e);
  });
}

template <class Seq, class F>
inline auto reduce_serial(Seq A, const F& f) -> typename Seq::T {
  using T = typename Seq::T;
  T r = A[0];
  for (size_t j = 1; j < A.size(); j++) r = f(r, A[j]);
  return r;
}

template <class Seq, class F>
inline auto reduce(Seq A, const F& f, flags fl = no_flag) -> typename Seq::T {
  using T = typename Seq::T;
  size_t n = A.size();
  size_t l = num_blocks(n, _block_size);
  if (l <= 1 || (fl & fl_sequential)) return reduce_serial(A, f);
  sequence<T> Sums(l);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    Sums[i] = reduce_serial(A.slice(s, e), f);
  });
  T r = reduce_serial(Sums, f);
  return r;
}

// sums I with + (can be any type with + defined)
template <class Seq>
inline auto reduce_add(Seq I, flags fl = no_flag) -> typename Seq::T {
  using T = typename Seq::T;
  auto add = [](T x, T y) { return x + y; };
  return reduce(I, add, fl);
}

// sums I with max
template <class Seq>
inline auto reduce_max(Seq I, flags fl = no_flag) -> typename Seq::T {
  using T = typename Seq::T;
  auto add = [](T x, T y) { return std::max(x, y); };
  return reduce(I, add, fl);
}

// sums I with max
template <class Seq>
inline auto reduce_min(Seq I, flags fl = no_flag) -> typename Seq::T {
  using T = typename Seq::T;
  auto add = [](T x, T y) { return std::min(x, y); };
  return reduce(I, add, fl);
}

// sums I with max
template <class Seq>
inline auto reduce_xor(Seq I, flags fl = no_flag) -> typename Seq::T {
  using T = typename Seq::T;
  auto add = [](T x, T y) { return x ^ y; };
  return reduce(I, add, fl);
}

const flags fl_scan_inclusive = (1 << 4);

// serial scan with combining function f and start value zero
// fl_scan_inclusive indicates it includes the current location
template <class In_Seq, class Out_Seq, class F>
inline auto scan_serial(In_Seq In, Out_Seq Out, const F& f,
                        typename In_Seq::T zero, flags fl = no_flag) ->
    typename In_Seq::T {
  using T = typename In_Seq::T;
  T r = zero;
  size_t n = In.size();
  bool inclusive = fl & fl_scan_inclusive;
  if (inclusive) {
    for (size_t i = 0; i < n; i++) {
      r = f(r, In[i]);
      Out.update(i, r);
    }
  } else {
    for (size_t i = 0; i < n; i++) {
      T t = In[i];
      Out.update(i, r);
      r = f(r, t);
    }
  }
  return r;
}

// parallel version of scan_serial -- see comments above
template <class In_Seq, class Out_Seq, class F>
inline auto scan(In_Seq In, Out_Seq Out, const F& f, typename In_Seq::T zero,
                 flags fl = no_flag) -> typename In_Seq::T {
  using T = typename In_Seq::T;
  size_t n = In.size();
  size_t l = num_blocks(n, _block_size);
  if (l <= 2 || fl & fl_sequential) return scan_serial(In, Out, f, zero, fl);
  sequence<T> Sums(l);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    Sums[i] = reduce_serial(In.slice(s, e), f);
  });
  T total = scan_serial(Sums, Sums, f, zero, 0);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    scan_serial(In.slice(s, e), Out.slice(s, e), f, Sums[i], fl);
  });
  return total;
}

template <class In_Seq, class Out_Seq>
inline auto scan_add(In_Seq In, Out_Seq Out, flags fl = no_flag) ->
    typename In_Seq::T {
  using T = typename In_Seq::T;
  auto add = [](T x, T y) { return x + y; };
  return scan(In, Out, add, (T)0, fl);
}

template <class Seq>
inline size_t sum_flags_serial(Seq I) {
  size_t r = 0;
  for (size_t j = 0; j < I.size(); j++) r += I[j];
  return r;
}

template <class In_Seq, class Bool_Seq>
inline auto pack_serial(In_Seq In, Bool_Seq Fl,
                        typename In_Seq::T* _Out = nullptr)
    -> sequence<typename In_Seq::T> {
  using T = typename In_Seq::T;
  size_t n = In.size();
  size_t m = sum_flags_serial(Fl);
  T* Out = (_Out) ? _Out : new_array_no_init<T>(m);
  size_t k = 0;
  for (size_t i = 0; i < n; i++)
    if (Fl[i]) assign_uninitialized(Out[k++], In[i]);
  return std::move(sequence<T>(Out, m, (_Out) ? false : true));
}

template <class In_Seq, class Bool_Seq>
inline void pack_serial_at(In_Seq In, typename In_Seq::T* Out, Bool_Seq Fl) {
  size_t k = 0;
  for (size_t i = 0; i < In.size(); i++)
    if (Fl[i]) assign_uninitialized(Out[k++], In[i]);
}

template <class In_Seq, class Bool_Seq>
inline auto pack(In_Seq In, Bool_Seq Fl, flags fl = no_flag,
                 typename In_Seq::T* _Out = nullptr)
    -> sequence<typename In_Seq::T> {
  using T = typename In_Seq::T;
  size_t n = In.size();
  size_t l = num_blocks(n, _block_size);
  if (l <= 1 || fl & fl_sequential) {
    return std::move(pack_serial(In, Fl.slice(0, In.size()), _Out));
  }
  sequence<size_t> Sums(l);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    Sums[i] = sum_flags_serial(Fl.slice(s, e));
  });
  size_t m = scan_add(Sums, Sums);
  T* Out = (_Out) ? _Out : new_array_no_init<T>(m);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    pack_serial_at(In.slice(s, e), Out + Sums[i], Fl.slice(s, e));
  });
  return std::move(sequence<T>(Out, m, (_Out) ? false : true));
}

template <class Idx_Type, class Bool_Seq>
inline sequence<Idx_Type> pack_index(Bool_Seq Fl, flags fl = no_flag) {
  auto identity = [](size_t i) { return (Idx_Type)i; };
  return pack(make_sequence<Idx_Type>(Fl.size(), identity), Fl, fl);
}

template <class Idx_Type, class D, class F>
inline sequence<std::tuple<Idx_Type, D> > pack_index_and_data(
    F& f, size_t size, flags fl = no_flag) {
  auto identity = [&](size_t i) {
    return std::make_tuple((Idx_Type)i, std::get<1>(f(i)));
  };
  auto flgs_in =
      make_sequence<bool>(size, [&](size_t i) { return std::get<0>(f(i)); });
  return pack(make_sequence<std::tuple<Idx_Type, D> >(size, identity), flgs_in,
              fl);
}

template <class In_Seq, class Pred>
inline auto filter_serial(In_Seq In, Pred p, flags fl = no_flag,
                          typename In_Seq::T* _Out = nullptr)
    -> sequence<typename In_Seq::T> {
  using T = typename In_Seq::T;
  size_t n = In.size();
  auto Fl = make_sequence<bool>(n, [&](size_t i) { return p(In[i]); });
  size_t m = sum_flags_serial(Fl);
  T* Out = (_Out) ? _Out : new_array_no_init<T>(m);
  size_t k = 0;
  for (size_t i = 0; i < n; i++) {
    if (p(In[i])) {
      assign_uninitialized(Out[k++], In[i]);
    }
  }
  auto ret = sequence<T>(Out, m, (_Out) ? false : true);
  return std::move(ret);
}

template <class In_Seq, class Pred>
inline sequence<typename In_Seq::T> filter(In_Seq In, Pred p,
                                           flags fl = no_flag,
                                           typename In_Seq::T* _Out = nullptr) {
  using T = typename In_Seq::T;
  size_t n = In.size();
  size_t l = num_blocks(n, _block_size);
  if (l <= 1 || fl & fl_sequential) {
    return std::move(filter_serial(In, p, fl, _Out));
  }
  auto Flags = sequence<bool>(n);
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                  { Flags[i] = (bool)p(In[i]); });

  auto ret = pack(In, Flags, fl, _Out);
  return std::move(ret);
}

template <class In_Seq, class Out_Seq, class Char_Seq>
inline std::pair<size_t, size_t> split_three(In_Seq In, Out_Seq Out,
                                             Char_Seq Fl, flags fl = no_flag) {
  size_t n = In.size();
  size_t l = num_blocks(n, _block_size);
  sequence<size_t> Sums0(l);
  sequence<size_t> Sums1(l);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    Sums0[i] = Sums1[i] = 0;
    for (size_t j = s; j < e; j++) {
      if (Fl[j] == 0)
        Sums0[i]++;
      else if (Fl[j] == 1)
        Sums1[i]++;
    }
  });
  size_t m0 = scan_add(Sums0, Sums0);
  size_t m1 = scan_add(Sums1, Sums1);
  sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
    size_t c0 = Sums0[i];
    size_t c1 = m0 + Sums1[i];
    size_t c2 = m0 + m1 + (s - Sums0[i] - Sums1[i]);
    for (size_t j = s; j < e; j++) {
      if (Fl[j] == 0)
        Out[c0++] = In[j];
      else if (Fl[j] == 1)
        Out[c1++] = In[j];
      else
        Out[c2++] = In[j];
    }
  });
  return std::make_pair(m0, m1);
}

template <class T, class Pred>
inline size_t filter_seq(T* in, T* out, size_t n, Pred p) {
  size_t k = 0;
  for (size_t i = 0; i < n; i++)
    if (p(in[i])) out[k++] = in[i];
  return k;
}

// Faster for a small number in output (about 40% or less)
// Destroys the input.   Does not need a bool array.
template <class T, class PRED>
inline size_t filterf(T* In, T* Out, size_t n, PRED p) {
  size_t b = _F_BSIZE;
  if (n < b) return filter_seq(In, Out, n, p);
  size_t l = num_blocks(n, b);
  size_t* Sums = new_array_no_init<size_t>(l + 1);
  par_for(0, l, 1, [&] (size_t i) {
    size_t s = i * b;
    size_t e = std::min(s + b, n);
    size_t k = s;
    for (size_t j = s; j < e; j++) {
      if (p(In[j])) In[k++] = In[j];
    }
    Sums[i] = k - s;
  });
  auto isums = sequence<size_t>(Sums, l);
  size_t m = scan_add(isums, isums);
  Sums[l] = m;
  par_for(0, l, 1, [&] (size_t i) {
    T* I = In + i * b;
    T* O = Out + Sums[i];
    for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
      O[j] = I[j];
    }
  });
  pbbs::free_array(Sums);
  return m;
}

// Faster for a small number in output (about 40% or less)
// Destroys the input.   Does not need a bool array.
template <class T, class PRED, class OUT>
inline size_t filterf(T* In, size_t n, PRED p, OUT out, size_t out_off) {
  size_t b = _F_BSIZE;
  if (n < b) {
    size_t k = out_off;
    for (size_t i = 0; i < n; i++) {
      if (p(In[i])) out(k++, In[i]);
    }
    return k - out_off;
  }
  size_t l = num_blocks(n, b);
  size_t* Sums = new_array_no_init<size_t>(l + 1);
  par_for(0, l, 1, [&] (size_t i) {
    size_t s = i * b;
    size_t e = std::min(s + b, n);
    size_t k = s;
    for (size_t j = s; j < e; j++) {
      if (p(In[j])) In[k++] = In[j];
    }
    Sums[i] = k - s;
  });
  auto isums = sequence<size_t>(Sums, l);
  size_t m = scan_add(isums, isums);
  Sums[l] = m;
  par_for(0, l, 1, [&] (size_t i) {
    T* I = In + i * b;
    size_t si = out_off + Sums[i];
    for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
      out(si + j, I[j]);
    }
  });
  pbbs::free_array(Sums);
  return m;
}

template <class T, class PRED>
inline size_t filterf_and_clear(T* In, T* Out, size_t n, PRED p, T& empty,
                                size_t* Sums) {
  size_t b = _F_BSIZE;
  if (n < b) return filter_seq(In, Out, n, p);
  size_t l = num_blocks(n, b);
  b = num_blocks(n, l);
  par_for(0, l, 1, [&] (size_t i) {
    size_t s = i * b;
    size_t e = std::min(s + b, n);
    size_t k = s;
    for (size_t j = s; j < e; j++) {
      if (p(In[j])) {
        In[k] = In[j];
        if (j > k) {
          In[j] = empty;
        }
        k++;
      }
    }
    Sums[i] = k - s;
  });
  auto isums = sequence<size_t>(Sums, l);
  size_t m = scan_add(isums, isums);
  Sums[l] = m;
  par_for(0, l, 1, [&] (size_t i) {
    T* I = In + i * b;
    T* O = Out + Sums[i];
    for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
      O[j] = I[j];
      I[j] = empty;
    }
  });
  return m;
}

// template <class T, class F>
// void split_serial(T* OutFalse, T* OutTrue, bool* Fl,
// 		    size_t s, size_t e, const F& f) {
//   size_t kT = 0;
//   size_t kF = 0;
//   for (size_t i=s; i < e; i++)
//     if (Fl[i]) { assign_uninitialized(OutTrue[kT], f(i)); ++kT; }
//     else { assign_uninitialized(OutFalse[kF], f(i)); ++kF; }
// }

// // Given a boolean array, splits so false (0) elements are at the bottom
// // and true (1) elements are at the top of the output (of lenght e-s).
// // As usual s is a start index, e is an end index and
// // f is a function of type [s,e-1) -> T
// template <class T, class F>
//   int split(T* Out, bool*  Fl, size_t s, size_t e, const F& f) {
//   size_t l = nblocks(e-s, _F_BSIZE);
//   size_t *sums_true = new_array_no_init<size_t>(l);
//   blocked_for (i, s, e, _F_BSIZE,
// 		 sums_true[i] = sum_flags_serial(Fl+s, e-s););
//   auto sums = sequence<size_t>(sums_true, l);
//   size_t num_true = scan_add(sums, sums);
//   size_t num_false = (e - s) - num_true;
//   T* out_true = Out + num_false;
//   blocked_for(i, s, e, _F_BSIZE,
// 		split_serial(Out + _F_BSIZE*i - sums_true[i],
// 			    out_true + sums_true[i],
// 			    Fl, s, e, f););
//   pbbs::free_array(sums_true);
//   return num_false;
// }

// template <class T>
// size_t split(T* In, T* Out, bool* Fl, size_t n) {
//   auto get = [In] (size_t i) {return In[i];};
//   return split(Out, Fl, (size_t) 0, n, get);
// }
}  // namespace pbbs
