#pragma once

#include <iostream>
#include <ctype.h>
#include <memory>
#include <stdlib.h>
#include <type_traits>
#include <atomic>
#include <cstring>

#include "parallel.h"

using std::cout;
using std::endl;

void* my_alloc(size_t);
void my_free(void*);

template <typename Lf, typename Rf >
static void par_do_if(bool do_parallel, Lf left, Rf right, bool cons=false) {
  if (do_parallel) par_do(left, right, cons);
  else {left(); right();}
}

template <typename Lf, typename Mf, typename Rf >
inline void par_do3(Lf left, Mf mid, Rf right) {
  auto left_mid = [&] () {par_do(left,mid);};
  par_do(left_mid, right);
}

template <typename Lf, typename Mf, typename Rf >
static void par_do3_if(bool do_parallel, Lf left, Mf mid, Rf right) {
  if (do_parallel) par_do3(left, mid, right);
  else {left(); mid(); right();}
}

#if defined(__APPLE__)
inline void* aligned_alloc(size_t a, size_t n) {return malloc(n);}
#else
#ifdef USEMALLOC
#include <malloc.h>
struct __mallopt {
  __mallopt() {
    mallopt(M_MMAP_MAX,0);
    mallopt(M_TRIM_THRESHOLD,-1);
  }
};

__mallopt __mallopt_var;
#endif
#endif

namespace pbbs {

  struct empty {};

  typedef uint32_t flags;
  const flags no_flag = 0;
  const flags fl_sequential = 1;
  const flags fl_debug = 2;
  const flags fl_time = 4;
  const flags fl_conservative = 8;
  const flags fl_inplace = 16;

  template<typename T>
  inline void assign_uninitialized(T& a, const T& b) {
    using TT = typename std::remove_volatile<T>::type;
    new (static_cast<void*>((TT*)std::addressof(a))) TT(b);
  }

  template<typename T>
  inline void move_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(std::move(b));
  }

  // a 32-bit hash function
  uint32_t hash32(uint32_t a);

  uint32_t hash32_2(uint32_t a);

  uint32_t hash32_3(uint32_t a);

  uint64_t hash64(uint64_t u);

  // a slightly cheaper, but possibly not as good version
  // based on splitmix64
  uint64_t hash64_2(uint64_t x);

  // Does not initialize the array
  template<typename E>
  E* new_array_no_init(size_t n, bool touch_pages=false) { //true) {
    // pads in case user wants to allign with cache lines
    size_t line_size = 64;
    size_t bytes = ((n * sizeof(E))/line_size + 1)*line_size;
    //E* r = (E*) aligned_alloc(line_size, bytes);
    E* r = (E*) my_alloc(bytes);
    if (r == NULL) {fprintf(stderr, "Cannot allocate space: %lu bytes", bytes); exit(1);}
    //parallel_for (size_t i = 0; i < bytes; i = i + (1 << 21)) ((bool*) r)[i] = 0;
    return r;
  }

  // Initializes in parallel
  template<typename E>
  E* new_array(size_t n) {
    E* r = new_array_no_init<E>(n);
    if (!std::is_trivially_default_constructible<E>::value) {
    //if (!std::is_default_constructible<E>::value) {
      if (n > 2048) {
	auto f = [&] (size_t i) { new ((void*) (r+i)) E;};
	parallel_for(0, n, f);
      }
      else
	for (size_t i = 0; i < n; i++) new ((void*) (r+i)) E;
    }
    return r;
  }

  void free_array(void* a);

  // Destructs in parallel
  template<typename E>
  void delete_array(E* A, size_t n) {
    // C++14 -- suppored by gnu C++11
    if (!std::is_trivially_destructible<E>::value) {
      //if (!std::is_destructible<E>::value) {
      if (n > 2048) {
	auto f = [&] (size_t i) {A[i].~E();};
	parallel_for(0, n, f);
      } else for (size_t i = 0; i < n; i++) A[i].~E();
    }
    using NVE = typename std::remove_volatile<E>::type;
    NVE* AA = (NVE*)A;
    my_free(AA);
  }

  template <typename ET>
  inline bool atomic_compare_and_swap(ET* a, ET oldval, ET newval) {
    if (sizeof(ET) == 1) {
      uint8_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 4) {
      uint32_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 8) {
      uint64_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 16) {
      __int128 r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap_16(reinterpret_cast<__int128*>(a), r_oval, r_nval);
    } else {
      std::cout << "Bad CAS Length" << sizeof(ET) << std::endl;
      exit(0);
    }
  }

  template <typename ET>
  inline bool atomic_compare_and_swap(volatile ET* a, ET oldval, ET newval) {
    if (sizeof(ET) == 1) {
      uint8_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint8_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 4) {
      uint32_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint32_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 8) {
      uint64_t r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint64_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 16) {
      __int128 r_oval, r_nval;
      std::memcpy(&r_oval, &oldval, sizeof(ET));
      std::memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap_16(reinterpret_cast<volatile __int128*>(a), r_oval, r_nval);
    } else {
      std::cout << "Bad CAS Length" << sizeof(ET) << std::endl;
      exit(0);
    }
  }

  template <typename E, typename EV>
  inline E fetch_and_add(E *a, EV b) {
    volatile E newV, oldV;
    do {oldV = *a; newV = oldV + b;}
    while (!atomic_compare_and_swap(a, oldV, newV));
    return oldV;
  }

  template <typename E, typename EV>
  inline void write_add(E *a, EV b) {
    //volatile E newV, oldV;
    E newV, oldV;
    do {oldV = *a; newV = oldV + b;}
    while (!atomic_compare_and_swap(a, oldV, newV));
  }

  template <typename E, typename EV>
  inline void write_add(std::atomic<E> *a, EV b) {
    //volatile E newV, oldV;
    E newV, oldV;
    do {oldV = a->load(); newV = oldV + b;}
    while (!std::atomic_compare_exchange_strong(a, &oldV, newV));
  }

  template <typename ET, typename F>
  inline bool write_min(ET *a, ET b, F less) {
    ET c; bool r=0;
    do c = *a;
    while (less(b,c) && !(r=atomic_compare_and_swap(a,c,b)));
    return r;
  }


  template <typename ET, typename F>
  inline bool write_min(volatile ET *a, ET b, F less) {
    ET c; bool r=0;
    do c = *a;
    while (less(b,c) && !(r=atomic_compare_and_swap(a,c,b)));
    return r;
  }

  template <typename ET, typename F>
  inline bool write_min(std::atomic<ET> *a, ET b, F less) {
    ET c; bool r=0;
    do c = a->load();
    while (less(b,c) && !(r=std::atomic_compare_exchange_strong(a, &c, b)));
    return r;
  }

  template <typename ET, typename F>
  inline bool write_max(ET *a, ET b, F less) {
    ET c; bool r=0;
    do c = *a;
    while (less(c,b) && !(r=atomic_compare_and_swap(a,c,b)));
    return r;
  }

  template <typename ET, typename F>
  inline bool write_max(volatile ET *a, ET b, F less) {
    ET c; bool r=0;
    do c = *a;
    while (less(c,b) && !(r=atomic_compare_and_swap(a,c,b)));
    return r;
  }

  template <typename ET, typename F>
  inline bool write_max(std::atomic<ET> *a, ET b, F less) {
    ET c; bool r=0;
    do c = a->load();
    while (less(c,b) && !(r=std::atomic_compare_exchange_strong(a, &c, b)));
    return r;
  }

  // returns the log base 2 rounded up (works on ints or longs or unsigned versions)
  template <class T>
  size_t log2_up(T i) {
    size_t a=0;
    T b=i-1;
    while (b > 0) {b = b >> 1; a++;}
    return a;
  }

  size_t granularity(size_t n);

  void assert_str(int cond, const std::string& s);

}
