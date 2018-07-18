#pragma once

#include <limits.h>
#include "../lib/macros.h"

// size of edge-offsets.
// If the number of edges is more than sizeof(MAX_UINT),
// you should set the LONG flag on the command line.
#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX
#else
typedef int intT;
typedef unsigned int uintT;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX
#endif

// edge size macros.
// If the number of vertices is more than sizeof(MAX_UINT)
// you should set the EDGELONG flag on the command line.
#if defined(EDGELONG)
typedef long intE;
typedef unsigned long uintE;
#define INT_E_MAX LONG_MAX
#define UINT_E_MAX ULONG_MAX
#else
typedef int intE;
typedef unsigned int uintE;
#define INT_E_MAX INT_MAX
#define UINT_E_MAX UINT_MAX
#endif

// edgemap_sparse_blocked granularity macro
constexpr const size_t kEMBlockSize = 8000;

// ======= compression macros and constants =======
constexpr const size_t PARALLEL_DEGREE = 1000;
// Take care in pushing this threshold too high; vertices with degree <
// pack_threshold stack allocate these bytes.
constexpr const size_t PD_PACK_THRESHOLD = 10000;

// Each vertex larger than PD_PACK_THRESHOLD is allocated deg(v) / kTempSpaceConstant
constexpr const size_t kTemporarySpaceConstant = 10;
typedef unsigned char uchar;

#define LAST_BIT_SET(b) (b & (0x80))
#define EDGE_SIZE_PER_BYTE 7

#if !defined(PD) && !defined(AMORTIZEDPD)
#define compression byte
#else
#ifdef AMORTIZEDPD
#define compression bytepd_amortized
#else
#define compression bytepd
#endif
#endif
