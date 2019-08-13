#pragma once

#include "ligra.h"
#include "utils/stats.h"

#include "Yoshida.h"

static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

template <class G>
double t_mis_yoshida(G& GA, commandLine P) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, YoshidaMIS(GA));
  return t;
}

