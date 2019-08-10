// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Code based on Julian's implementation of Steriou et al. (see Ligra repo)
// Implementation of the shortcutting label propagation algorithm from
// the paper "Shortcutting Label Propagation for Distributed Connected
// Components", WSDM 2018.


#pragma once

#include "../benchmark/LDD.h"
#include "ligra.h"
#include "utils/contract.h"
#include "utils/stats.h"

namespace stergiou_shortcut {
  struct CC_Shortcut {
    uintE* IDs, *prevIDs;
    CC_Shortcut(uintE* _IDs, uintE* _prevIDs) :
      IDs(_IDs), prevIDs(_prevIDs) {}
    inline bool operator () (uintE i) {
      uintE l = IDs[IDs[i]];
      if(IDs[i] != l) IDs[i] = l;
      if(prevIDs[i] != IDs[i]) {
        prevIDs[i] = IDs[i];
        return 1; }
      else return 0;
    }
  };

  struct CC_F {
    uintE* IDs, *prevIDs;
    CC_F(uintE* _IDs, uintE* _prevIDs) :
      IDs(_IDs), prevIDs(_prevIDs) {}
    inline bool update(uintE s, uintE d){ //Update function writes min ID
      uintE origID = IDs[d];
      if(IDs[s] < origID) {
        IDs[d] = min(origID,IDs[s]);
      } return 1; }
    inline bool updateAtomic (uintE s, uintE d) { //atomic Update
      uintE origID = IDs[d];
      writeMin(&IDs[d],IDs[s]);
      return 1;
    }
    inline bool cond (uintE d) { return cond_true(d); } //does nothing
  };

  template <class G>
  void CC_stergiou_shortcutting(G& GA) {
    long n = GA.n, m = GA.m;
    uintE* IDs = pbbs::new_array_no_init<uintE>(n), *prevIDs = pbbs::new_array_no_init<uintE>(n);
    parallel_for(0, n, [&] (size_t i) {prevIDs[i] = i; IDs[i] = i;}); //initialize unique IDs

    auto all = pbbs::sequence<bool>(n, true);
    vertexSubset All(n,n,all.to_array()); //all vertices
    auto active = pbbs::sequence<bool>(n, true);
    vertexSubset Active(n,n,active.to_array()); //initial frontier contains all vertices

    while(!Active.isEmpty()){ //iterate until IDS converge
      edgeMap(GA, Active, CC_F(IDs,prevIDs),-1,no_output);
      vertexSubset output = vertexFilter(All,CC_Shortcut(IDs,prevIDs));
      Active.del();
      Active = output;
    }
    Active.del(); All.del(); pbbs::free_array(IDs); pbbs::free_array(prevIDs);
  }
}; // namespace stergiou_shortcut
