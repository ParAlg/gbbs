#pragma once

#include "sequence.h"
#include "range_min.h"
#include "get_time.h"

// Calculate the Longest Common Prefix (LCP) array from a Suffix array
namespace pbbs {

  //  The suffix array SA are indices into the string s
  template <class Seq1, class Seq2>
  auto lcp(Seq1 const &s, Seq2 const &SA) -> sequence<typename Seq2::value_type> {
    using Uint = typename Seq2::value_type;
    timer t("LCP", false);
    timer t2("LCP total", false);
    size_t len = 111;
    size_t n = SA.size();

    // compare first len characters of adjacent strings from SA.
    sequence<Uint> L(n-1, [&] (size_t i) {
	size_t j = 0;
	while (s[SA[i]+j] == s[SA[i+1]+j] && s[SA[i+1]+j] != 0 && j < len) j++;
	return (j < len) ? j : n;
      });
    t.next("head");

    // keep indices for which we do not yet know their LCP (i.e. LCP >= len)
    sequence<Uint> remain = pack_index<Uint>(map<bool>(L, [&] (Uint l) {
	  return l == n;}));
    t.next("pack");

    if (remain.size() == 0) { t2.next("total"); return  L;}

    // an inverse permutation for SA
    sequence<Uint> ISA(n);
    parallel_for(0, n, [&] (size_t i) {
	ISA[SA[i]] = i;});
    t.next("inverse");

    // repeatedly double len determining LCP by joining next len chars
    // invariant: before i^th round L contains correct LCPs less than len
    //            and n for the rest of them
    //            remain holds indices of the rest of them (i.e., LCP[i] >= len)
    //      after round, len = 2*len and invariant holds for the new len
    do {
      auto rq = make_range_min(L, std::less<Uint>(), 111);
      t.next("make range");

      // see if next len chars resolves LCP
      // set L for those that do, and keep those that do not for next round
      remain = filter(remain, [&] (Uint i) {
	  if (SA[i] + len >= n) {L[i] = len; return false;};
	  Uint i1 = ISA[SA[i]+len];
	  Uint i2 = ISA[SA[i+1]+len];
	  size_t l = L[rq.query(i1, i2-1)];
	  if (l < len) {L[i] = len + l; return false;}
	  else return true;
	});
      t.next("filter");
      len *= 2;
    } while (remain.size() > 0);

    t2.next("total");
    return L;
  }

}  // namespace pbbs
