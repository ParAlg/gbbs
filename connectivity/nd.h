#include <iostream>
#include <limits.h>
#include "sequence.h"
#include "gettime.h"
#include "graph.h"
#include "parallel.h"
#include "unionFind.h"
using namespace std;

struct notMax { bool operator() (uintT i) {return i < UINT_T_MAX;}};

// Assumes root is negative
// Not making parent array volatile improves
// performance and doesn't affect correctness
inline uintT find(uintT i, uintT * parent) {
  //if (parent[i] < 0) return i;
  uintT j = i; //parent[i];
  if (parent[j] == UINT_T_MAX) return j;
  do j = parent[j];
  while (parent[j] < UINT_T_MAX);
  //note: path compression can happen in parallel in the same tree, so
  //only link from smaller to larger to avoid cycles
  uintT tmp;
  while((tmp=parent[i])<j){ parent[i]=j; i=tmp;}
  return j;
}

pair<uintT*, uintT> st(edgeArray<uintT> EA){
  edge<uintT>* E = EA.E;
  long m = EA.nonZeros;
  long n = EA.numRows;
  uintT *parents = newA(uintT,n);
  parallel_for (long i=0; i < n; i++) parents[i] = UINT_T_MAX;
  uintT *hooks = newA(uintT,n);
  parallel_for (long i=0; i < n; i++) hooks[i] = UINT_T_MAX;

  parallel_for (long i = 0; i < m; i++) {
    uintT u = E[i].u, v = E[i].v;
    while(1){
      u = find(u,parents);
      v = find(v,parents);
      if(u == v) break;
      if(u > v) swap(u,v);
      //if successful, store the ID of the edge used in hooks[u]
      if(hooks[u] == UINT_T_MAX && __sync_bool_compare_and_swap(&hooks[u],UINT_T_MAX,i)){
	parents[u]=v;
	break;
      }
    }
  }


  //uncomment this code for connected component labels instead of
  //spanning forest edges
  /*
  parallel_for(long i=0;i<n;i++) {
    hooks[i] = find(i,parents);
  }
  free(parents); free(hooks);
  return pair<uintT*,uintT>(NULL, 0);
  */

  //get the IDs of the edges in the spanning forest
  _seq<uintT> stIdx = sequence::filter((uintT*) hooks, n, notMax());

  free(parents); free(hooks);
  cout<<"nInSt = "<<stIdx.n<<endl;
  return pair<uintT*,uintT>(stIdx.A, stIdx.n);
}
