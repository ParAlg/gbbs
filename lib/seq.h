#pragma once

#include "utilities.h"

template <typename E>
struct sequence {
public:
  using T = E;

  sequence() {}

  // copy constructor
  sequence(const sequence& a) : s(a.s), e(a.e), allocated(false) {
  }

  // move constructor
  sequence(sequence&& b)
    : s(b.s), e(b.e), allocated(b.allocated) {
    b.s = b.e = NULL; b.allocated = false;}

  // copy assignment
  sequence& operator = (const sequence& b) {
    if (this != &b) {
      clear(); s = b.s; e = b.e; allocated = false;
    }
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
    : s(NULL), e(NULL), allocated(false) {
    if (n > 0) {
      s = pbbs::new_array<E>(n);
      e = s + n;
      allocated = true;
    }
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
    parallel_for(0, n, f, 0);
  };

  template <typename Func>
  sequence(const size_t n, Func f)
    : s(pbbs::new_array_no_init<E>(n)), allocated(true) {
    e = s + n;
    auto g = [&] (size_t i) {
      T x = f(i);
      new ((void*) (s+i)) T(x);};
    parallel_for(0, n, g, 0);
  }

  sequence(E* s, const size_t n, bool allocated = false)
    : s(s), e(s + n), allocated(allocated) {};

  sequence(E* s, E* e, bool allocated = false)
    : s(s), e(e), allocated(allocated) {};

  ~sequence() { clear();}

  sequence copy(sequence& a) {
    return sequence<E>(e-s, [&] (size_t i) {return a[i];});
  }

  E& operator[] (const size_t i) const {return s[i];}

  sequence slice(size_t ss, size_t ee) const {
    return sequence(s + ss, s + ee);
  }

  sequence slice() const {
    return sequence(s, e);
  }

  void update(const size_t i, T& v) {
    s[i] = v;
  }

  size_t size() const { return e - s;}

  sequence as_sequence() const {
    return sequence(s, e);
  }

  T* begin() const {return s;}
  T* end() const {return e;}
  bool is_allocated() const {return allocated;}
  void set_allocated(bool a) {allocated = a;}

  // TODO(ldhulipa): crutch before porting over code
  T* start() const {return s;}
  E& operator() (const size_t i) const {return s[i];}
  T* get_array() {
    allocated = false;
    return s;
  }

  void clear() {
    if (allocated) pbbs::delete_array<E>(s,e-s);
    allocated = false;
    s = e = NULL;
  }

  //private:
  E *s; // = NULL;
  E *e; // = NULL;
  bool allocated = false;
};

template <typename E, typename F>
struct func_sequence {
  using T = E;
  func_sequence(size_t n, F _f) : f(_f), s(0), e(n) {};
  func_sequence(size_t s, size_t e, F _f) : f(_f), s(s), e(e) {};
  const T operator[] (size_t i) const {return (f)(i+s);}
  func_sequence<T,F> slice(size_t ss, size_t ee) const {
    return func_sequence<T,F>(s+ss,s+ee,f); }
  func_sequence<T,F> slice() const {
    return func_sequence<T,F>(s,e,f); }
  size_t size() const { return e - s;}
  sequence<T> as_sequence() const {
    return sequence<T>::tabulate(e-s, [&] (size_t i) {return (f)(i+s);});
  }
  // TODO(ldhulipa): crutch before porting over code
  const T operator() (size_t i) const {return (f)(i+s);}
private:
  const F f;
  const size_t s, e;
};

// used so second template argument can be inferred
template <class E, class F>
func_sequence<E,F> make_sequence (size_t n, F f) {
  return func_sequence<E,F>(n,f);
}

template <class E>
inline sequence<E> make_sequence(E* A, size_t n) {
  return sequence<E>(A, n);
}
