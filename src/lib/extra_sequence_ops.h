#pragma once

#include "lib/seq.h"
#include "lib/sequence_ops.h"

namespace pbbs {

  template <class Seq>
  inline auto reduce_max(Seq I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    auto add = [](T x, T y) { return std::max(x, y); };
    return reduce(I, add, fl);
  }

  template <class Seq>
  inline auto reduce_min(Seq I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    auto add = [](T x, T y) { return std::min(x, y); };
    return reduce(I, add, fl);
  }

  template <class Seq>
  inline auto reduce_xor(Seq I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    auto add = [](T x, T y) { return x ^ y; };
    return reduce(I, add, fl);
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

}
