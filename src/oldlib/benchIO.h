// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
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

#include <math.h>
#include "lib/utilities.h"
#include "utils.h"

#include <algorithm>
#include <cstring>
#include <iostream>

namespace benchIO {

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long n;          // total number of characters
  char* Chars;     // array storing all strings
  long m;          // number of substrings
  char** Strings;  // pointers to strings (all should be null terminated)
  words() {}
  words(char* C, long nn, char** S, long mm)
      : n(nn), Chars(C), m(mm), Strings(S) {}
  void del() {
    pbbs::free_array(Chars);
    pbbs::free_array(Strings);
  }
};

inline bool isSpace(char c) {
  switch (c) {
    case '\r':
    case '\t':
    case '\n':
    case 0:
    case ' ':
      return true;
    default:
      return false;
  }
}

struct toLong {
  long operator()(bool v) { return (long)v; }
};

// parallel code for converting a string to words
inline words stringToWords(char* Str, long n) {
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    if (isSpace(Str[i])) Str[i] = 0;
  });

  // mark start of words
  bool* FL = pbbs::new_array_no_init<bool>(n);
  FL[0] = Str[0];
  parallel_for_bc(i, 1, n, (n > pbbs::kSequentialForThreshold),
                  { FL[i] = Str[i] && !Str[i - 1]; });

  // offset for each start of word
  ligra_utils::_seq<long> Off = ligra_utils::seq::packIndex<long>(FL, n);
  long m = Off.n;
  long* offsets = Off.A;

  // pointer to each start of word
  char** SA = newA(char*, m);
  parallel_for_bc(j, 0, m, (m > pbbs::kSequentialForThreshold),
                  { SA[j] = Str + offsets[j]; });

  pbbs::free_array(offsets);
  pbbs::free_array(FL);
  return words(Str, n, SA, m);
}

inline int writeStringToFile(char* S, long n, char* fileName) {
  std::ofstream file(fileName, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << "\n";
    return 1;
  }
  file.write(S, n);
  file.close();
  return 0;
}

inline int xToStringLen(long a) { return 21; }
inline void xToString(char* s, long a) { sprintf(s, "%ld", a); }

inline int xToStringLen(unsigned long a) { return 21; }
inline void xToString(char* s, unsigned long a) { sprintf(s, "%lu", a); }

inline uint xToStringLen(uint a) { return 12; }
inline void xToString(char* s, uint a) { sprintf(s, "%u", a); }

inline int xToStringLen(int a) { return 12; }
inline void xToString(char* s, int a) { sprintf(s, "%d", a); }

inline int xToStringLen(double a) { return 18; }
inline void xToString(char* s, double a) { sprintf(s, "%.11le", a); }

inline int xToStringLen(char* a) { return strlen(a) + 1; }
inline void xToString(char* s, char* a) { sprintf(s, "%s", a); }

template <class A, class B>
inline int xToStringLen(std::pair<A, B> a) {
  return xToStringLen(a.first) + xToStringLen(a.second) + 1;
}
template <class A, class B>
inline void xToString(char* s, std::pair<A, B> a) {
  int l = xToStringLen(a.first);
  xToString(s, a.first);
  s[l] = ' ';
  xToString(s + l + 1, a.second);
}

struct notZero {
  bool operator()(char A) { return A > 0; }
};

template <class T>
inline ligra_utils::_seq<char> arrayToString(T* A, long n) {
  long* L = pbbs::new_array_no_init<long>(n);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { L[i] = xToStringLen(A[i]) + 1; });
  long m = ligra_utils::seq::scan(L, L, n, ligra_utils::addF<long>(), (long)0);
  char* B = pbbs::new_array_no_init<char>(m);
  parallel_for_bc(j, 0, m, (m > pbbs::kSequentialForThreshold), { B[j] = 0; });
  parallel_for_bc(i, 0, n - 1, (n - 1 > pbbs::kSequentialForThreshold), {
    xToString(B + L[i], A[i]);
    B[L[i + 1] - 1] = '\n';
  });
  xToString(B + L[n - 1], A[n - 1]);
  B[m - 1] = '\n';
  pbbs::free_array(L);
  char* C = pbbs::new_array_no_init<char>(m + 1);
  long mm = ligra_utils::seq::filter(B, C, m, notZero());
  C[mm] = 0;
  pbbs::free_array(B);
  return ligra_utils::_seq<char>(C, mm);
}

template <class T>
inline void writeArrayToStream(std::ofstream& os, T* A, long n) {
  long BSIZE = 1000000;
  long offset = 0;
  while (offset < n) {
    // Generates a string for a sequence of size at most BSIZE
    // and then wrties it to the output stream
    ligra_utils::_seq<char> S =
        arrayToString(A + offset, std::min(BSIZE, n - offset));
    os.write(S.A, S.n);
    S.del();
    offset += BSIZE;
  }
}

template <class T>
inline void writeArrayToStream(std::ofstream& os, T* A, size_t n) {
  size_t BSIZE = 1000000;
  size_t offset = 0;
  while (offset < n) {
    // Generates a string for a sequence of size at most BSIZE
    // and then wrties it to the output stream
    std::cout << "Writing offset = " << offset << "\n";
    ligra_utils::_seq<char> S =
        arrayToString(A + offset, std::min(BSIZE, n - offset));
    os.write(S.A, S.n);
    S.del();
    offset += BSIZE;
  }
}

template <class T>
inline int writeArrayToFile(std::string header, T* A, long n, char* fileName) {
  std::ofstream file(fileName, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << "\n";
    return 1;
  }
  file << header << "\n";
  writeArrayToStream(file, A, n);
  file.close();
  return 0;
}

inline ligra_utils::_seq<char> readStringFromFile(char* fileName) {
  std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << "\n";
    abort();
  }
  long end = file.tellg();
  file.seekg(0, std::ios::beg);
  long n = end - file.tellg();
  char* bytes = pbbs::new_array_no_init<char>(n + 1);
  file.read(bytes, n);
  file.close();
  return ligra_utils::_seq<char>(bytes, n);
}
};  // namespace benchIO
