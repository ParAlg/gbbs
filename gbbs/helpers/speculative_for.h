// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#include <assert.h>
#include <limits.h>

#include "gbbs/bridge.h"

namespace gbbs {

template <class intT>
struct reservation {
  intT r;
  reservation() : r(std::numeric_limits<intT>::max()) {}
  void reserve(intT i) {
    gbbs::write_min(&r, i, [](intT left, intT right) { return left < right; });
  }
  bool reserved() { return (r < std::numeric_limits<intT>::max()); }
  void reset() { r = std::numeric_limits<intT>::max(); }
  bool check(intT i) { return (r == i); }
  bool checkReset(intT i) {
    if (r == i) {
      r = std::numeric_limits<intT>::max();
      return 1;
    } else
      return 0;
  }
};

template <class intT>
inline void reserveLoc(intT* x, intT i) {
  gbbs::write_min<intT>(x, i);
}

// granularity is some constant.
template <class intT, class S>
inline intT eff_for(S step, intT s, intT e, intT granularity, bool hasState = 1,
                    long maxTries = std::numeric_limits<long>::max()) {
  intT maxRoundSize = (e - s) / granularity + 1;
  intT currentRoundSize = maxRoundSize;

  auto I = sequence<intT>(maxRoundSize);
  auto Inext = sequence<intT>(maxRoundSize);
  auto keep = sequence<bool>(maxRoundSize);

  intT round = 0;
  intT numberDone = s;      // number of iterations done
  intT numberKeep = 0;      // number of iterations to carry to next round
  intT totalProcessed = 0;  // number done including wasted tries

  while (numberDone < e) {
    if (round++ > maxTries) {
      std::cout << "speculative_for: too many iterations, increase maxTries"
                << "\n";
      abort();
    }

    intT size = std::min(currentRoundSize, (intT)(e - numberDone));
    totalProcessed += size;

    parallel_for(0, size, [&](size_t i) {
      if (i >= numberKeep) I[i] = numberDone + i;
      keep[i] = step.reserve(I[i]);
    });

    parallel_for(0, size, [&](size_t i) {
      if (keep[i]) keep[i] = !step.commit(I[i]);
    });

    // keep iterations that failed for next round. Written into Inext
    numberKeep = parlay::pack_out(I.cut(0, size), keep, make_slice(Inext));
    //      seq.set_allocated(false);
    //      numberKeep = seq.size();
    numberDone += size - numberKeep;

    std::swap(I, Inext);

    // adjust round size based on number of failed attempts
    if (float(numberKeep) / float(size) < .1) {
      currentRoundSize = std::min(currentRoundSize * 2, maxRoundSize);
    }
  }
  return totalProcessed;
}

template <class intT, class S>
inline intT speculative_for(S step, intT s, intT e, intT granularity,
                            bool hasState = 1, long maxTries = -1) {
  if (maxTries < 0) {
    maxTries = 100 + 200 * granularity;
  }
  intT maxRoundSize = (e - s) / granularity + 1;
  intT currentRoundSize = maxRoundSize;

  auto I = sequence<intT>(maxRoundSize);
  auto Inext = sequence<intT>(maxRoundSize);
  auto keep = sequence<bool>(maxRoundSize);

  intT round = 0;
  intT numberDone = s;      // number of iterations done
  intT numberKeep = 0;      // number of iterations to carry to next round
  intT totalProcessed = 0;  // number done including wasted tries

  while (numberDone < e) {
    if (round++ > maxTries) {
      std::cout << "speculative_for: too many iterations, increase maxTries"
                << "\n";
      abort();
    }

    intT size = std::min(currentRoundSize, (intT)(e - numberDone));
    totalProcessed += size;

    parallel_for(0, size, [&](size_t i) {
      if (i >= numberKeep) I[i] = numberDone + i;
      keep[i] = step.reserve(I[i]);
    });

    parallel_for(0, size, [&](size_t i) {
      if (keep[i]) keep[i] = !step.commit(I[i]);
    });

    // keep iterations that failed for next round. Written into Inext
    numberKeep = parlay::pack_out(I.cut(0, size), keep, make_slice(Inext));
    numberDone += size - numberKeep;

    std::swap(I, Inext);

    // adjust round size based on number of failed attempts
    if (float(numberKeep) / float(size) > .2)
      currentRoundSize =
          std::max(currentRoundSize / 2,
                   std::max(maxRoundSize / 64 + 1, (intT)numberKeep));
    else if (float(numberKeep) / float(size) < .1)
      currentRoundSize = std::min(currentRoundSize * 2, maxRoundSize);
    //    std::cout << size << " : " << numberKeep << " : " << numberDone <<
    //    "\n";
  }
  return totalProcessed;
}

}  // namespace gbbs
