# pbbslib

Contains the PBBS (Problem Based Benchmark Suite) library.  This is a
collection of parallel algorithms and other support for parallelism
written in C++.  It supports a reasonably wide variety of algorithms
that run on shared memory multicores.  In addition to the algorithms
it supports a parallel scheduleer and a memory allocator.  It only
needs the c++ standard threads utilities.

The library encourages a more functional style of programming than the
algorithms available in the standard template library.  It borrows
some ideas from Boost.

******************************************************

## PBBS includes:
  - a library of algorithms, mostly based on sequences and ranges
  - a scheduler (based on arbitrarily nested fork-join parallelism)
  - a memory allocator (optimized for parallelism)
  - parallel random number generation

### Algorithms (all parallel):

  - reduce, scan, scan_inplace
  - pack, filter, pack_index
  - merge
  - random_shuffle
  - histogram
  - integer_sort, counting_sort, sort, sort_inplace, stable_sort
  - collect_reduce
  - kth_smallest

It also includes the following, which are loosely based on the standard
template library.   However none of them mutate
their arguments, but rather are copy based.

  - all_of, any_of, none_of, count, count_if
  - find, find_if, find_if_not, find_first_of, find_end, search
  - mismatch, adjacent_find
  - equal, lexicographical_compare,
  - unique, remove_if
  - min_element, max_element, min_max_element
  - reverse, rotate,
  - is_sorted, is_sorted_until, is_partitioned
  
### Utilities
  - scheduler
  - parallel random number generator
  - memory allocator
  - monoid
  - timer
  - hashing
  
### Concurrency
  - stack
  
******************************************************

## Seqs and Ranges

Much of The library revolves around the concepts of Seq and Range.

###  Seq

A Seq is a type supporting at least a.size(), a.split(),
a.split(start,end) and and a[i], where the last can return an rval.
The split functions return a Seq (possibly also a Range).  It must
also contain a value_type.  No iterator needs to be associated with a
Seq.  A general Seq cannot be mutated (none of its functionality
supports mutation).  Functions that take a Seq as an argument almost
always take it as a const reference.

### Range
A Range can be thought of as a pair of iterators marking the start and
end of the range.  It must support all the operations of Seq (it is
also a Seq), and in addition, it must also support a.begin(), a.end()
and must have a type iterator.  Importantly Range's should not have
any memory associated with them for their underlining container.  They
are meant to reference a range (perhaps the full range) of an
underlining container.  This is different from the concept of Range in
boost.

******************************************************

## sequence, range and delayed_sequence

The library supports three specific types of Seq and Range.

### sequence<T>

A sequence<T> is similar to a vector<T> in that it is associated with
the underlying memory for its contents.  Unlike vector<T> it is
generally meant to be immutable, although it does support mutation by
first extracting a range.  It is designed so all its internal
operations are parallel---including initialization of elements,
destruction of elements, copy assignment and copy construction.  It
supports the following functionality:
    
    typename value_type (shorthand T below)
    sequence<T>(size_t n, IntegerFunc f) :
        applies f to i : [0,n) creating a sequence of length n
    sequence<T>() : creates a zero length sequence
    sequence<T>(size_t n) : creates a sequence of length n
    sequence<T>(size_t n, T v) :
        creates an sequence of length n initializing with v
    sequence<T>(Range r) :
        creates a sequence from the range r
    a.size() -> size_t : size of a
    a.slice() -> range<T*> : a Range encompasing full range of a
    a.slice(size_t s, size_t e) -> range<T*> :
        a Range encompasing [s,e) of a.
    a.begin() : shorthand for a.slice().begin()
    a.end() : shorthand for a.slice().end()
    a[size_t i] -> T& : reference to i-th value of a
    a.swap(sequence<T> b) -> void : swaps contents of a and b
    a.clear() -> clears the contents setting length to zero.    

If passed by rvalue reference it is moved clearing the rvalue.  If
passed by value, it is copied.   

A sequence is a Seq but not a Range.

******************************************************

### range<Iterator>

A range<Iterator> suppors random access into range, and slicing into
subranges.   
    
    typename value_type (shorthand T below)
    typename iterator, where iterator::value_type = value_type
       it must be a random access iterator
    range<I>(I begin, I end) : the range from begin to end
    a.size() -> size_t : size of the range
    a.slice() -> range<I> : returns a 
    a.slice(size_t s, size_t e) -> range<I> :
        the subrange from a.begin()+s to a.begin()+e, 
    a.begin() -> I : an iterator to the begining of the range
    a.end() -> I : an iterator to the end of the range
    a[size_t i] -> value_type& : reference to i-th value of a

A range<Iterator> is both a Seq and a Range.  

******************************************************

### delayed_sequence<T, IntegerFunc>

A delayed_sequence<T, IntegerFunc> is associated with a function that
emits the i-th element.  It is hence immutable, and there is no
notion of a reference to its elements.  It supports:

    typename value_type (shorthand T below)
    delayed_sequence(size_t n, IntegerFunc f) :
      creates a delayed sequence of length n, where the i-th element is f(i)
      the function is applied on demand, not when the delayed sequence is created
    delayed_sequence(size_t n, T v) :
      the function is the constant function returning v for every i
    a[size_t i] -> value_type : f(i)
    a.slice() -> delayed_sequence<T, IntegerFunc> : a
    a.slice(size_t s, size_t e) -> delayed_sequence<T, IntegerFunc> :
        the range from [f(s), f(e))

A delayed_sequence is a Seq but not a Range

 
