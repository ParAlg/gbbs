#pragma once

// legacy memory allocation macro
#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))

// scan/filter macros; used by sequence implementations
#define _SCAN_LOG_BSIZE 10
#define _SCAN_BSIZE (1 << _SCAN_LOG_BSIZE)
#define _F_BSIZE (2*_SCAN_BSIZE)

