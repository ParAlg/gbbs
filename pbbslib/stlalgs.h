#pragma once

#include "kth_smallest.h"
#include "sample_sort.h"
#include "sequence_ops.h"

namespace pbbs {

template <class IntegerPred>
size_t count_if_index(size_t n, IntegerPred p) {
  auto BS =
      pbbs::delayed_seq<size_t>(n, [&](size_t i) -> size_t { return p(i); });
  size_t r = pbbs::reduce(BS, pbbs::addm<size_t>());
  return r;
}

template <class IntegerPred>
size_t find_if_index(size_t n, IntegerPred p, size_t granularity = 1000) {
  {
    size_t i;
    for (i = 0; i < std::min(granularity, n); i++)
      if (p(i)) return i;
    if (i == n) return n;
  }
  size_t start = granularity;
  while (start < n) {
    size_t end = std::min(n, start + granularity);
    auto f = [&](size_t i) -> size_t { return p(i + start) ? i + start : n; };
    size_t i = pbbs::reduce(delayed_seq<size_t>(end - start, f), minm<size_t>());
    if (i < n) return i;
    start += granularity;
    granularity *= 2;
  }
  return n;
}

template <class Seq, class UnaryFunction>
void for_each(Seq const &S, UnaryFunction f) {
  parallel_for(S.size_t(), [&](size_t i) { f(S[i]); });
}

template <class Seq, class UnaryPred>
size_t count_if(Seq const &S, UnaryPred p) {
  return count_if_index(S.size(), [&](size_t i) { return p(S[i]); });
}

template <class Seq, class T>
size_t count(Seq const &S, T const &value) {
  return count_if_index(S.size(), [&](size_t i) { return S[i] == value; });
}

template <class Seq, class UnaryPred>
bool all_of(Seq const &S, UnaryPred p) {
  return count_if(S, p) == S.size();
}

template <class Seq, class UnaryPred>
bool any_of(Seq const &S, UnaryPred p) {
  return count_if(S, p) > 1;
}

template <class Seq, class UnaryPred>
bool none_of(Seq const &S, UnaryPred p) {
  return count_if(S, p) == 0;
}

template <class Seq, class UnaryPred>
size_t find_if(Seq const &S, UnaryPred p) {
  return find_if_index(S.size(), [&](size_t i) { return p(S[i]); });
}

template <class Seq, class T>
size_t find(Seq const &S, T const &value) {
  return find_if(S, [&](T x) { return x == value; });
}

template <class Seq, class UnaryPred>
size_t find_if_not(Seq const &S, UnaryPred p) {
  return find_if_index(S.size(), [&](size_t i) { return !p(S[i]); });
}

template <class Seq, class BinaryPred>
size_t adjacent_find(Seq const &S, BinaryPred pred) {
  return find_if_index(S.size() - 1,
                       [&](size_t i) { return S[i] == S[i + 1]; });
}

template <class Seq, class BinaryPred>
size_t mismatch(Seq const &S1, Seq const &S2, BinaryPred pred) {
  return find_if_index(std::min(S1.size(), S2.size()),
                       [&](size_t i) { return S1[i] != S2[i]; });
}

template <class Seq, class BinaryPred>
size_t search(Seq const &S1, Seq const &S2, BinaryPred pred) {
  return find_if_index(S1.size() - S2.size() + 1, [&](size_t i) {
    size_t j;
    for (j = 0; j < S2.size(); j++)
      if (S1[i + j] != S2[j]) break;
    return (j == S2.size());
  });
}

template <class Seq, class BinaryPred>
size_t find_end(Seq const &S1, Seq const &S2, BinaryPred pred) {
  size_t n1 = S1.size();
  size_t n2 = S2.size();
  size_t idx = find_if_index(S1.size() - S2.size() + 1, [&](size_t i) {
    size_t j;
    for (j = 0; j < n2; j++)
      if (S1[(n1 - i - n2) + j] != S2[j]) break;
    return (j == S2.size());
  });
  return n1 - idx - n2;
}

template <class Seq1, class Seq2, class BinaryPred>
bool equal(Seq1 s1, Seq2 s2, BinaryPred p) {
  return count_if_index(s1.size(), [&](size_t i) { return p(s1[i], s2[i]); });
}

template <class Seq1, class Seq2, class BinaryPred>
bool equal(Seq1 s1, Seq2 s2) {
  return count_if_index(s1.size(), [&](size_t i) { return s1[i] == s2[i]; });
}

template <class Seq1, class Seq2, class Compare>
bool lexicographical_compare(Seq1 s1, Seq2 s2, Compare less) {
  auto s = delayed_seq(s1.size(), [&](size_t i) {
    return less(s1[i], s2[i]) ? -1 : (less(s2[i], s1[i]) ? 1 : 0);
  });
  auto f = [&](char a, char b) { return (a == 0) ? b : a; };
  return reduce(s, make_monoid(f, (char)0)) == -1;
}

template <class Seq, class Eql>
sequence<typename Seq::value_type> unique(Seq const &s, Eql eq) {
  sequence<bool> b(s.size(),
                   [&](size_t i) { return (i == 0) || !eq(s[i], s[i - 1]); });
  return pack(s, b);
}

// needs to return location, and take comparison
template <class Seq, class Compare>
size_t min_element(Seq const &S, Compare comp) {
  auto SS = delayed_seq<size_t>(S.size(), [&](size_t i) { return i; });
  auto f = [&](size_t l, size_t r) { return (!comp(S[r], S[l]) ? l : r); };
  return pbbs::reduce(SS, make_monoid(f, (size_t)S.size()));
}

template <class Seq, class Compare>
size_t max_element(Seq const &S, Compare comp) {
  using T = typename Seq::value_type;
  return min_element(S, [&](T a, T b) { return comp(b, a); });
}

template <class Seq, class Compare>
std::pair<size_t, size_t> minmax_element(Seq const &S, Compare comp) {
  size_t n = S.size();
  using P = std::pair<size_t, size_t>;
  auto SS = delayed_seq<P>(S.size(), [&](size_t i) { return P(i, i); });
  auto f = [&](P l, P r) -> P {
    return (P(!comp(S[r.first], S[l.first]) ? l.first : r.first,
              !comp(S[l.second], S[r.second]) ? l.second : r.second));
  };
  return pbbs::reduce(SS, make_monoid(f, P(n, n)));
}

template <class Seq>
sequence<typename Seq::value_type> reverse(Seq const &S) {
  size_t n = S.size();
  return sequence<typename Seq::value_type>(
      S.size(), [&](size_t i) { return S[n - i - 1]; });
}

template <class Seq>
sequence<typename Seq::value_type> rotate(Seq const &S, size_t r) {
  size_t n = S.size();
  return sequence<typename Seq::value_type>(S.size(), [&](size_t i) {
    size_t j = (i < r) ? n - r + i : i - r;
    return S[j];
  });
}

template <class Seq, class Compare>
bool is_sorted(Seq const &S, Compare comp) {
  auto B = delayed_seq<bool>(
      S.size() - 1, [&](size_t i) -> size_t { return f(S[i + 1], S[i]); });
  return (reduce(B, addm<size_t>()) != 0);
}

template <class Seq, class Compare>
size_t is_sorted_until(Seq const &S, Compare comp) {
  return find_if_index(S.size() - 1, [&](size_t i) { f(S[i + 1], S[i]); }) + 1;
}

template <class Seq, class UnaryPred>
size_t is_partitioned(Seq const &S, UnaryPred f) {
  auto B = delayed_seq<bool>(
      S.size() - 1, [&](size_t i) -> size_t { return !f(S[i + 1]) && S[i]; });
  return (reduce(B, addm<size_t>()) != 0);
}

template <class Seq, class UnaryPred>
size_t remove_if(Seq const &S, UnaryPred f) {
  using T = typename Seq::value_type;
  return filter(S, [&](T a) { return !f(a); });
}

template <class Seq, class Compare>
sequence<typename Seq::value_type> sort(Seq const &S, Compare less) {
  return sample_sort(S, less, false);
}

template <class Iter, class Compare>
void sort_inplace(range<Iter> A, const Compare &f) {
  sample_sort_inplace(A, f);
};

template <class Seq, class Compare>
sequence<typename Seq::value_type> stable_sort(Seq const &S, Compare less) {
  return sample_sort(S, less, true);
}

template <class Seq, class Compare>
sequence<typename Seq::value_type> remove_duplicates_ordered(Seq const &s,
                                                             Compare less) {
  using T = typename Seq::value_type;
  return unique(stable_sort(s, less),
                [&](T a, T b) { return !less(a, b) && !less(b, a); });
}

template <class Seq1, class Seq2>
sequence<typename Seq1::value_type> append(Seq1 const &s1, Seq2 const &s2) {
  using T = typename Seq1::value_type;
  size_t n1 = s1.size();
  return sequence<T>(n1 + s2.size(),
                     [&](size_t i) { return (i < n1) ? s1[i] : s2[i - n1]; });
}

template <class Index, class BoolSeq>
std::pair<sequence<Index>, Index> enumerate(BoolSeq const &s) {
  return scan(
      delayed_seq<Index>(s.size(), [&](size_t i) -> Index { return s[i]; }),
      addm<Index>());
}
}
// template <class Seq, class Compare>
// std::pair<sequence<typename Seq::value_type>, size_t>
// partition(Seq const &S, UnaryPred f) {
//   using T = typename Seq::value_type;
//   auto r = count_sort(S, delayed_seq<size_t>(n, [&] (size_t i) {return
//   f(S[i]);}), 1);
//   return r[0]
// }

/*

Most of these are from the boost libraries, but the boost versions take ranges
(as in slices).
copy
copy_if
copy_n
move
fill
fill_n (takes iter + len)
transform (parallel_for)
generate (applies same function to each location)
generate_n
replace
replace_if
swap_ranges
shift_left
shift_right
partition
partition_copy
stable_partition
sort
partial_sort
partial_sort_copy
stable_sort
nth_element (finds n_th element and paritions on it)
merge
inplace_merge
includes
set_difference
set_intersection
set_symmetric_difference
set_union
is_heap
is_heap_until
adjacent_difference
reduce (requires commutativity)
exclusive_scan
inclusive_scan
transform_reduce (map reduce)
transform_exclusive_scan
transform_inclusive_scan
uninitialized_copy
uninitialized_copy_n
uninitialized_fill
uninitialized_fill_n
uninitialized_move
uninitialized_move_n
uninitialized_default_construct
uninitialized_default_construct_n
uninitialized_value_construct
uninitialized_value_construct_n
destroy
destroy_n

Following do not have parallel versions
random_shuffle
sample
make_heap
sort_heap
is_permutaton
iota
accumulate (does not require associativity
inner_product

*/

/*
Duck Typing

constraints and concepts
-fconcepts in gcc >= 6.1

examples from: www.stroustrup.com/good_concepts.pdf

template<typename T>
concept bool Sequence =
  requires(T t) {
    typename Value_type<T>;
    typename Iterator_of<T>;
    { begin(t) } -> Iterator_of<T>;
    { end(t) } -> Iterator_of<T>;
    requires Input_iterator<Iterator_of<T>>;
    requires Same_type<Value_type<T>,Value_type<Iterator_of<T>>>;
};

template<typename T>
concept bool Equality_comparable =
  requires (T a, T b) {
    { a == b } -> bool;
    { a != b } -> bool;
};

template<typename T>
concept bool Sortable =
  Sequence<T> &&
  Random_access_iterator<Iterator_of<T>> &&
  Less_than_comparable<Value_type<T>>;

template<typename T, typename U>
concept bool Equality_comparable =
  requires (T a, U b) {
    { a == b } -> bool;
    { a != b } -> bool;
    { b == a } -> bool;
    { b != a } -> bool;
};


template<typename T>
concept bool Number = requires(T a, T b) {
  { a+b } -> T;
  { a-b } -> T;
  { a*b } -> T;
  { a/b } -> T;
  { -a } -> T;
  { a+=b } -> T&;
  { a-=b } -> T&;
  { a*=b } -> T&;
  { a/=b } -> T&;
  { T{0} };
  ...
};

template <Sequence S, typename T>
    requires Equality_comparable<Value_type<S>, T>
Iterator_of<S> find(S& seq, const T& value);

Or, equivalently:

template <typename S, typename T>
    requires Sequence<S> && Equality_comparable<Value_type<S>, T>
Iterator_of<S> find(S& seq, const T& value);

If only one constraint:

void sort(Sortable& s);


Problem of accidental matching.  The longer the list the better.
Avoid single property concepts.

class My_number { ...  };
static_assert(Number<My_number>);

*/

/*
ranges are similar to slices
come out of boost, although in there ranges would cover both slices and
sequences, and vectors, ...

"The motivation for the Range concept is that there are many useful
Container-like types that do not meet the full requirements of
Container, and many algorithms that can be written with this reduced
set of requirements. In particular, a Range does not necessarily

    - own the elements that can be accessed through it,
    - have copy semantics,

": From boost  (forward range, bidirectional range, random access range)

following from:
https://www.fluentcpp.com/2017/01/12/ranges-stl-to-the-next-level/
ideas seem to be mostly from boost

transform_iterator -> delayed_seq
filter_iterator -> skips over elements not satisfying a predicate, only makes
sense for linear iterator, not random access iterator.

view adaptors

view::transform adaptor

std::vector numbers = { 1, 2, 3, 4, 5 };
auto range = numbers | view::transform(multiplyBy2);
ranges::accumulate(numbers | view::transform(multiplyBy2), 0);
ranges::accumulate(numbers | view::filter(isEven), 0);
ranges::accumulate(numbers | view::filter(isEven) |
view::transform(multiplyBy2), 0);

This has many similarities to lazy sequences in haskell.  There is a pipelining
from one to the next.   I guess it could not only filter, but actually add new
elelments (e.g. duplicate).

Note that not exactly like delayed seq since they transform an
existing iterator rather than starting from indices.  Although I guess
can start with iota.
*/
