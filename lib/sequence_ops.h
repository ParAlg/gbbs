// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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
#include "utilities.h"
#include "seq.h"
#include "monoid.h"

namespace pbbs {

  constexpr const size_t _log_block_size = 10;
  constexpr const size_t _block_size = (1 << _log_block_size);

  inline size_t num_blocks(size_t n, size_t block_size) {
    if (n == 0) return 0;
    else return (1 + ((n)-1)/(block_size));}

  template <class F>
  void sliced_for(size_t n, size_t block_size, const F& f, flags fl = no_flag) {
    size_t l = num_blocks(n, block_size);
    auto body = [&] (size_t i) {
      size_t s = i * block_size;
      size_t e = std::min(s + block_size, n);
      f(i, s, e);
    };
    parallel_for(0, l, body, 1, 0 != (fl & fl_conservative));
  }

  template <class Seq, class Monoid>
  auto reduce_serial(Seq const &A, Monoid m) -> typename Seq::T {
    using T = typename Seq::T;
    T r = A[0];
    for (size_t j=1; j < A.size(); j++) r = m.f(r,A[j]);
    return r;
  }

  template <class Seq, class Monoid>
  auto reduce(Seq const &A, Monoid m, flags fl = no_flag)
    -> typename Seq::T
  {
    using T = typename Seq::T;
    size_t n = A.size();
    size_t block_size = std::max(_block_size, 4 * (size_t) ceil(sqrt(n)));
    size_t l = num_blocks(n, block_size);
    if (l == 0) return m.identity;
    if (l == 1 || (fl & fl_sequential))
      return reduce_serial(A, m);
    sequence<T> Sums(l);
    sliced_for (n, block_size,
		[&] (size_t i, size_t s, size_t e)
		{ Sums[i] = reduce_serial(A.slice(s,e), m);});
    T r = reduce(Sums, m);
    return r;
  }

  // sums I with + (can be any type with + defined)
  template <class Seq>
  auto reduce_add(Seq const &I, flags fl = no_flag)
    -> typename Seq::T
  {
    using T = typename Seq::T;
    return reduce(I, addm<T>(), fl);
  }

  const flags fl_scan_inclusive = (1 << 4);

  template <class Seq, class F>
  auto reduce_serialo(Seq A, const F& f) -> typename Seq::T {
    using T = typename Seq::T;
    T r = A[0];
    for (size_t j=1; j < A.size(); j++) r = f(r,A[j]);
    return r;
  }

    // serial scan with combining function f and start value zero
  // fl_scan_inclusive indicates it includes the current location
  template <class In_Seq, class Out_Seq, class F>
  auto scan_serialo(In_Seq In, Out_Seq Out,
		   const F& f, typename In_Seq::T zero,
		   flags fl = no_flag)  -> typename In_Seq::T
  {
    using T = typename In_Seq::T;
    T r = zero;
    size_t n = In.size();
    bool inclusive = fl & fl_scan_inclusive;
    if (inclusive) {
      for (size_t i = 0; i < n; i++) {
	r = f(r,In[i]);
	Out.update(i,r);
      }
    } else {
      for (size_t i = 0; i < n; i++) {
	T t = In[i];
	Out.update(i,r);
	r = f(r,t);
      }
    }
    return r;
  }

  // parallel version of scan_serial -- see comments above
  template <class In_Seq, class Out_Seq, class F>
  auto scano(In_Seq In, Out_Seq Out,
	    const F& f, typename In_Seq::T zero,
	    flags fl = no_flag)  -> typename In_Seq::T
  {
    using T = typename In_Seq::T;
    size_t n = In.size();
    size_t l = num_blocks(n,_block_size);
    if (l <= 2 || fl & fl_sequential)
      return scan_serialo(In, Out, f, zero, fl);
    sequence<T> Sums(l);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{ Sums[i] = reduce_serialo(In.slice(s,e),f);});
    T total = scan_serialo(Sums, Sums, f, zero, 0);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{ scan_serialo(In.slice(s,e), Out.slice(s,e), f, Sums[i], fl);});
    return total;
  }

  template <class In_Seq, class Out_Seq>
  auto scan_add(In_Seq In, Out_Seq Out, flags fl = no_flag)
    -> typename In_Seq::T
  {
    using T = typename In_Seq::T;
    auto add = [] (T x, T y) {return x + y;};
    return scano(In, Out, add, (T) 0, fl);
  }

  template <class In_Seq, class Out_Seq, class Monoid>
  auto scan_serial(In_Seq const &In, Out_Seq &Out,
		   Monoid const &m, typename In_Seq::T offset,
		   flags fl = no_flag)  -> typename In_Seq::T
  {
    using T = typename In_Seq::T;
    T r = offset;
    size_t n = In.size();
    bool inclusive = fl & fl_scan_inclusive;
    if (inclusive) {
      for (size_t i = 0; i < n; i++) {
	r = m.f(r,In[i]);
	Out.update(i,r);
      }
    } else {
      for (size_t i = 0; i < n; i++) {
	T t = In[i];
	Out.update(i,r);
	r = m.f(r,t);
      }
    }
    return r;
  }

  // parallel version of scan_serial -- see comments above
  template <class In_Seq, class Out_Seq, class Monoid>
  auto scan_(In_Seq const &In, Out_Seq &Out, Monoid const &m,
	     flags fl = no_flag) -> typename In_Seq::T
  {
    using T = typename In_Seq::T;
    size_t n = In.size();
    size_t l = num_blocks(n,_block_size);
    if (l <= 2 || fl & fl_sequential)
      return scan_serial(In, Out, m, m.identity, fl);
    sequence<T> Sums(l);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{ Sums[i] = reduce_serial(In.slice(s,e), m);});
    T total = scan_serial(Sums, Sums, m, m.identity, 0);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{ auto O = Out.slice(s,e);
		  scan_serial(In.slice(s,e), O, m, Sums[i], fl);});
    return total;
  }

  template <class In_Seq, class Monoid>
  auto scan_inplace(In_Seq &In, Monoid m, flags fl = no_flag)
    -> typename In_Seq::T
  { return scan_(In, In, m, fl); }

  template <class In_Seq, class Monoid>
  auto scan(In_Seq const &In, Monoid m, flags fl = no_flag)
    ->  std::pair<sequence<typename In_Seq::T>,typename In_Seq::T>
  {
    using T = typename In_Seq::T;
    sequence<T> Out(In.size());
    return std::make_pair(std::move(Out), scan_(In, Out, m, fl));
  }

  // do in place if rvalue reference to a sequence<T>
  template <class T, class Monoid>
  auto scan(sequence<T> &&In, Monoid m, flags fl = no_flag)
    ->  std::pair<sequence<T>, T> {
    sequence<T> Out = std::move(In);
    T total = scan_(Out, Out, m, fl);
    return std::make_pair(std::move(Out), total);
  }

  template <class Seq>
  size_t sum_bools_serial(Seq const &I) {
    size_t r = 0;
    for (size_t j=0; j < I.size(); j++) r += I[j];
    return r;
  }

  template <class In_Seq, class Bool_Seq>
  inline auto pack_serial(In_Seq const &In, Bool_Seq const &Fl,
                          typename In_Seq::T* _Out = nullptr)
      -> sequence<typename In_Seq::T> {
    using T = typename In_Seq::T;
    T* Out;
    size_t n = In.size();
    Out = (_Out) ? _Out : new_array_no_init<T>(sum_bools_serial(Fl));
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
      if (Fl[i]) assign_uninitialized(Out[k++], In[i]);
    return sequence<T>(Out, k, (_Out) ? false : true);
  }

  template <class In_Seq, class Bool_Seq>
  void pack_serial_at(In_Seq const &In, typename In_Seq::T* Out,
		      Bool_Seq const &Fl) {
    size_t k = 0;
    for (size_t i=0; i < In.size(); i++)
      if (Fl[i]) assign_uninitialized(Out[k++], In[i]);
  }

  template <class In_Seq, class Bool_Seq>
  inline auto pack(In_Seq const &In, Bool_Seq const &Fl, flags fl = no_flag,
                   typename In_Seq::T* _Out = nullptr)
      -> sequence<typename In_Seq::T> {
    using T = typename In_Seq::T;
    size_t n = In.size();
    size_t l = num_blocks(n, _block_size);
    if (l <= 1 || fl & fl_sequential) {
      return pack_serial(In, Fl.slice(0, In.size()), _Out);
    }
    sequence<size_t> Sums(l);
    sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
      Sums[i] = sum_bools_serial(Fl.slice(s, e));
    });
    size_t m = scan_inplace(Sums, addm<size_t>());
    T* Out = (_Out) ? _Out : new_array_no_init<T>(m);
    sliced_for(n, _block_size, [&](size_t i, size_t s, size_t e) {
	pack_serial(In.slice(s, e),  Fl.slice(s, e), Out + Sums[i]);
    });
    return sequence<T>(Out, m, (_Out) ? false : true);
  }

  template <class In_Seq, class F>
  auto filter(In_Seq const &In, F const &f, flags fl = no_flag,
	      typename In_Seq::T* _Out = nullptr)
    -> sequence<typename In_Seq::T>
  {
    using T = typename In_Seq::T;
    size_t n = In.size();
    size_t l = num_blocks(n,_block_size);
    sequence<size_t> Sums(l);
    sequence<bool> Fl(n);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{ size_t r = 0;
		  for (size_t j=s; j < e; j++)
		    r += (Fl[j] = f(In[j]));
		  Sums[i] = r;});
    size_t m;
    std::tie(Sums,m) = scan(std::move(Sums), addm<size_t>());
    T* Out = (_Out != nullptr) ? _Out : new_array_no_init<T>(m);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{ pack_serial(In.slice(s,e),
			      Fl.slice(s,e),
			      Out + Sums[i]);});
    return sequence<T>(Out,m,(_Out == nullptr));
  }

  template <class Idx_Type, class Bool_Seq>
  sequence<Idx_Type> pack_index(Bool_Seq const &Fl, flags fl = no_flag) {
    auto identity = [] (size_t i) {return (Idx_Type) i;};
    return pack(make_sequence<Idx_Type>(Fl.size(),identity), Fl, fl);
  }

  template <class In_Seq, class Out_Seq, class Char_Seq>
  std::pair<size_t,size_t> split_three(In_Seq const &In, Out_Seq &Out, Char_Seq const &Fl,
				       flags fl = no_flag) {
    size_t n = In.size();
    if (In.begin() == Out.begin()) {
      std::cout << "In and Out cannot be the same in split_tree" << std::endl;
      return std::pair<size_t,size_t>(0,0);
    }
    size_t l = num_blocks(n,_block_size);
    sequence<size_t> Sums0(l);
    sequence<size_t> Sums1(l);
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{
		  size_t c0 = 0;
		  size_t c1 = 0;
		  Sums0[i] = Sums1[i] = 0;
		  for (size_t j=s; j < e; j++) {
		    if (Fl[j] == 0) c0++;
		    else if (Fl[j] == 1) c1++;
		  }
		  Sums0[i] = c0;
		  Sums1[i] = c1;
		}, fl);
    size_t m0 = scan_inplace(Sums0, addm<size_t>());
    size_t m1 = scan_inplace(Sums1, addm<size_t>());
    sliced_for (n, _block_size,
		[&] (size_t i, size_t s, size_t e)
		{
		  size_t c0 = Sums0[i];
		  size_t c1 = m0 + Sums1[i];
		  size_t c2 = m0 + m1 + (s - Sums0[i] - Sums1[i]);
		  for (size_t j=s; j < e; j++) {
		    if (Fl[j] == 0) Out[c0++] = In[j];
		    else if (Fl[j] == 1) Out[c1++] = In[j];
		    else Out[c2++] = In[j];
		  }
		}, fl);
    return std::make_pair(m0,m1);
  }
}

