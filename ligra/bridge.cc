#include "bridge.h"

#include <stdio.h>

inline int t_to_stringlen(long a) { return 21; }
inline void type_to_string(char* s, long a) { sprintf(s, "%ld", a); }

inline int t_to_stringlen(unsigned long a) { return 21; }
inline void type_to_string(char* s, unsigned long a) { sprintf(s, "%lu", a); }

inline uint t_to_stringlen(uint a) { return 12; }
inline void type_to_string(char* s, uint a) { sprintf(s, "%u", a); }

inline int t_to_stringlen(int a) { return 12; }
inline void type_to_string(char* s, int a) { sprintf(s, "%d", a); }

inline int t_to_stringlen(double a) { return 18; }

inline int t_to_stringlen(char* a) { return strlen(a) + 1; }
inline void type_to_string(char* s, char* a) { sprintf(s, "%s", a); }

inline void type_to_string(char* s, double a) { sprintf(s, "%.11le", a); }
