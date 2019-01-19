#pragma once

#include "utilities.h"

template <typename E>
struct sequence {
public:
  using T = E;

  sequence() {}

  // copy constructor
  sequence(const sequence& a) : s(a.s), e(a.e), allocated(false) {}

  // move constructor
  sequence(sequence&& b)
    : s(b.s), e(b.e), allocated(b.allocated) {
    b.s = b.e = NULL; b.allocated = false;}

  // copy assignment
  sequence& operator = (const sequence& b) {
    if (this != &b) {clear(); s = b.s; e = b.e; allocated = false;}
    return *this;
  }

  // move assignment
  sequence& operator = (sequence&& b) {
    if (this != &b) {
      clear(); s = b.s; e = b.e; allocated = b.allocated;
      b.s = b.e = NULL; b.allocated=false;}
    return *this;
  }

  sequence(const size_t n)
    : s(pbbs::new_array<E>(n)), allocated(true) {
    e = s + n;
  };

  static sequence<E> alloc_no_init(const size_t n) {
    sequence<E> r;
    r.s = pbbs::new_array_no_init<E>(n);
    r.e = r.s + n;
    r.allocated = true;
    return r;
  };

  sequence(const size_t n, T v)
    : s(pbbs::new_array_no_init<E>(n, true)), allocated(true) {
    e = s + n;
    auto f = [=] (size_t i) {new ((void*) (s+i)) T(v);};
    //par_for(0, n, pbbs::granularity(n), f);
    parallel_for(0, n, f, 0);
    //parallel_for (size_t i=0; i < n; i++) f(i);
  };

  template <typename Func>
  sequence(const size_t n, const Func& f)
    : s(pbbs::new_array_no_init<E>(n)), allocated(true) {
    e = s + n;
    auto g = [&] (size_t i) {
      T x = f(i);
      new ((void*) (s+i)) T(x);};
    //par_for(0, n, pbbs::granularity(n), g);
    parallel_for(0, n, g, 0);
  }

  sequence(E* s, const size_t n, bool allocated = false)
    : s(s), e(s + n), allocated(allocated) {};

  sequence(E* s, E* e, bool allocated = false)
    : s(s), e(e), allocated(allocated) {};

  ~sequence() { clear();}

  // template <typename X, typename F>
  // static sequence<X> tabulateo(size_t n, F f) {
  //   X* r = pbbs::new_array_no_init<X>(n);
  //   parallel_for (size_t i = 0; i < n; i++)
  //     new ((void*) (r+i)) X(f(i));
  //   sequence<X> y(r,n);
  //   y.allocated = true;
  //   return y;
  // }

  template <typename F>
  static sequence<T> tabulate(size_t n, F f) {
    T* r = pbbs::new_array_no_init<T>(n);
    //parallel_for (size_t i = 0; i < n; i++)
    auto g = [&] (size_t i) {
      new ((void*) (r+i)) T(f(i));};
    par_for(0, n, pbbs::granularity(n), g);
    sequence<T> y(r,n);
    y.allocated = true;
    return y;
  }

  sequence copy(sequence& a) {
    return tabulate(e-s, [&] (size_t i) {return a[i];});
  }

  E& operator[] (const size_t i) const {return s[i];}
  E& operator() (const size_t i) const {return s[i];}

  sequence slice(size_t ss, size_t ee) {
    return sequence(s + ss, s + ee);
  }

  void update(const size_t i, T& v) {
    s[i] = v;
  }

  size_t size() { return e - s;}

  sequence as_sequence() {
    return sequence(s, e);
  }

  T* as_array() {return s;}
  T* start() {return s;}
  T* end() {return e;}
  T* get_array() {
    allocated=false;
    return s;
  }

  //private:
  void clear() {
    if (s && e) {
      if (allocated) pbbs::delete_array<E>(s,e-s);
      s = e = NULL;
    }
  }
  E *s; // = NULL;
  E *e; // = NULL;
  bool allocated = false;
};

template <typename E, typename F>
struct func_sequence {
  using T = E;
  func_sequence(size_t n, F& _f) : f(&_f), s(0), e(n) {};
  func_sequence(size_t s, size_t e, F& _f) : f(&_f), s(s), e(e) {};
  T operator[] (const size_t i) {return (*f)(i+s);}
  T operator() (const size_t i) {return (*f)(i+s);}
  func_sequence<T,F> slice(size_t ss, size_t ee) {
    return func_sequence<T,F>(s+ss,s+ee,*f); }
  size_t size() { return e - s;}
  sequence<T> as_sequence() {
    return sequence<T>::tabulate(e-s, [&] (size_t i) {return (*f)(i+s);});
  }
private:
  F *f;
  size_t s, e;
};

// used so second template argument can be inferred
template <class E, class F>
func_sequence<E,F> make_sequence (size_t n, F& f) {
  return func_sequence<E,F>(n,f);
}

template <class E>
inline sequence<E> make_sequence(E* A, size_t n) {
  return sequence<E>(A, n);
}
