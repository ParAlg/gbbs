#include "bridge.h"

#include <stdio.h>

namespace pbbslib {

int t_to_stringlen(long a) { return 21; }
void type_to_string(char* s, long a) { sprintf(s, "%ld", a); }

int t_to_stringlen(unsigned long a) { return 21; }
void type_to_string(char* s, unsigned long a) { sprintf(s, "%lu", a); }

uint t_to_stringlen(uint a) { return 12; }
void type_to_string(char* s, uint a) { sprintf(s, "%u", a); }

int t_to_stringlen(int a) { return 12; }
void type_to_string(char* s, int a) { sprintf(s, "%d", a); }

int t_to_stringlen(double a) { return 18; }

int t_to_stringlen(char* a) { return strlen(a) + 1; }
void type_to_string(char* s, char* a) { sprintf(s, "%s", a); }

void type_to_string(char* s, double a) { sprintf(s, "%.11le", a); }

}  // namespace pbbslib
