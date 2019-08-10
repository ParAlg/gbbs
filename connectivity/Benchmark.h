#pragma once

#include "ligra.h"
#include "utils/stats.h"

#include "GBBSCC.h"
#include "GBBSHybridCC.h"
#include "GBBSSF.h"
#include "Async.h"
#include "Rem.h"
#include "NdOpt.h"
#include "StergiouShortcutting.h"

static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

template <class G>
double t_gbbs_cc(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = gbbs_cc::CC(GA, beta, false, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class G>
double t_gbbs_sf(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto E = gbbs_sf::SpanningForest(GA, beta, false, permute));
  E.del();
  return t;
}

template <class G>
double t_async_cc(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto CC = AsyncCC(GA));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class G>
double t_async_sf(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto E = AsyncSF(GA));
  return t;
}

template <class G>
double t_rem_sf(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto E = RemSF(GA));
  return t;
}

template <class G>
double t_rem_cc(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto E = RemCC(GA));
  if (P.getOptionValue("-check")) {
    cc_check(correct, E);
  }
  return t;
}

template <class G>
double t_ndopt_cc(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto E = NdOptCC(GA));
  if (P.getOptionValue("-check")) {
    cc_check(correct, E);
  }
  return t;
}

template <class G>
double t_ndopt_sf(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto E = NdOptSF(GA));
  return t;
}

template <class G>
double t_gbbs_hybridcc(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = gbbs_hybridcc::CC(GA, beta, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class G>
double t_stergiou_shortcutting(G& GA, commandLine P, pbbs::sequence<uintE>& correct) {
  time(t, auto CC = stergiou_shortcut::CC_stergiou_shortcutting(GA));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}
