#pragma once
#include "utilities.h"
#include "quicksort.h" // needed for insertion_sort

namespace pbbs {

  template <class SeqA, class SeqB, class F> 
  void merge_sort(SeqA Out, SeqB In, const F& f, bool swap=0) {
    size_t n = In.size();
    if (n < 24) {
      pbbs::insertion_sort(In.as_array(), n, f);
      if (!swap)
	for (size_t i=0; i < n; i++) Out[i] = In[i]; //.update(i, In[i]);
      return;
    }
    size_t m = n/2;
    par_do_if(n > 64,
	   [&] () {merge_sort(Out.slice(0,m), In.slice(0,m), f, !swap);},
	   [&] () {merge_sort(Out.slice(m,n), In.slice(m,n), f, !swap);},
	   true);
    if (swap)
      pbbs::merge(Out.slice(0,m), Out.slice(m,n), In.slice(0,n), f, true);
    else
      pbbs::merge(In.slice(0,m), In.slice(m,n), Out.slice(0,n), f, true);
  }
}
