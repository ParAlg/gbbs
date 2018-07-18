#pragma once

#include <math.h>
#include "parse_command_line.h"
#include "utils.h"

#include <cstring>
#include <iostream>
#include <algorithm>

namespace benchIO {
  using namespace std;

  // A structure that keeps a sequence of strings all allocated from
  // the same block of memory
  struct words {
    long n; // total number of characters
    char* Chars;  // array storing all strings
    long m; // number of substrings
    char** Strings; // pointers to strings (all should be null terminated)
    words() {}
    words(char* C, long nn, char** S, long mm)
      : Chars(C), n(nn), Strings(S), m(mm) {}
    void del() {free(Chars); free(Strings);}
  };

  inline bool isSpace(char c) {
    switch (c)  {
    case '\r':
    case '\t':
    case '\n':
    case 0:
    case ' ' : return true;
    default : return false;
    }
  }

  struct toLong { long operator() (bool v) {return (long) v;} };

  // parallel code for converting a string to words
  words stringToWords(char *Str, long n) {
    parallel_for (long i=0; i < n; i++)
      if (isSpace(Str[i])) Str[i] = 0;

    // mark start of words
    bool *FL = newA(bool,n);
    FL[0] = Str[0];
    parallel_for (long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];

    // offset for each start of word
    _seq<long> Off = seq::packIndex<long>(FL, n);
    long m = Off.n;
    long *offsets = Off.A;

    // pointer to each start of word
    char **SA = newA(char*, m);
    parallel_for (long j=0; j < m; j++) SA[j] = Str+offsets[j];

    free(offsets); free(FL);
    return words(Str,n,SA,m);
  }

  int writeStringToFile(char* S, long n, char* fileName) {
    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file.write(S, n);
    file.close();
    return 0;
  }

  inline int xToStringLen(long a) { return 21;}
  inline void xToString(char* s, long a) { sprintf(s,"%ld",a);}

  inline int xToStringLen(unsigned long a) { return 21;}
  inline void xToString(char* s, unsigned long a) { sprintf(s,"%lu",a);}

  inline uint xToStringLen(uint a) { return 12;}
  inline void xToString(char* s, uint a) { sprintf(s,"%u",a);}

  inline int xToStringLen(int a) { return 12;}
  inline void xToString(char* s, int a) { sprintf(s,"%d",a);}

  inline int xToStringLen(double a) { return 18;}
  inline void xToString(char* s, double a) { sprintf(s,"%.11le", a);}

  inline int xToStringLen(char* a) { return strlen(a)+1;}
  inline void xToString(char* s, char* a) { sprintf(s,"%s",a);}

  template <class A, class B>
  inline int xToStringLen(pair<A,B> a) {
    return xToStringLen(a.first) + xToStringLen(a.second) + 1;
  }
  template <class A, class B>
  inline void xToString(char* s, pair<A,B> a) {
    int l = xToStringLen(a.first);
    xToString(s,a.first);
    s[l] = ' ';
    xToString(s+l+1,a.second);
  }

  struct notZero { bool operator() (char A) {return A > 0;}};

  template <class T>
  _seq<char> arrayToString(T* A, long n) {
    long* L = newA(long,n);
    {parallel_for(long i=0; i < n; i++) L[i] = xToStringLen(A[i])+1;}
    long m = seq::scan(L,L,n,addF<long>(),(long) 0);
    char* B = newA(char,m);
    parallel_for(long j=0; j < m; j++)
      B[j] = 0;
    parallel_for(long i=0; i < n-1; i++) {
      xToString(B + L[i],A[i]);
      B[L[i+1] - 1] = '\n';
    }
    xToString(B + L[n-1],A[n-1]);
    B[m-1] = '\n';
    free(L);
    char* C = newA(char,m+1);
    long mm = seq::filter(B,C,m,notZero());
    C[mm] = 0;
    free(B);
    return _seq<char>(C,mm);
  }

  template <class T>
  void writeArrayToStream(ofstream& os, T* A, long n) {
    long BSIZE = 1000000;
    long offset = 0;
    while (offset < n) {
      // Generates a string for a sequence of size at most BSIZE
      // and then wrties it to the output stream
      _seq<char> S = arrayToString(A+offset,min(BSIZE,n-offset));
      os.write(S.A, S.n);
      S.del();
      offset += BSIZE;
    }
  }

  template <class T>
  void writeArrayToStream(ofstream& os, T* A, size_t n) {
    size_t BSIZE = 1000000;
    size_t offset = 0;
    while (offset < n) {
      // Generates a string for a sequence of size at most BSIZE
      // and then wrties it to the output stream
      cout << "Writing offset = " << offset << endl;
      _seq<char> S = arrayToString(A + offset, min(BSIZE, n - offset));
      os.write(S.A, S.n);
      S.del();
      offset += BSIZE;
    }
  }

  template <class T>
    int writeArrayToFile(string header, T* A, long n, char* fileName) {
    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file << header << endl;
    writeArrayToStream(file, A, n);
    file.close();
    return 0;
  }

  _seq<char> readStringFromFile(char *fileName) {
    ifstream file (fileName, ios::in | ios::binary | ios::ate);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      abort();
    }
    long end = file.tellg();
    file.seekg (0, ios::beg);
    long n = end - file.tellg();
    char* bytes = newA(char,n+1);
    file.read (bytes,n);
    file.close();
    return _seq<char>(bytes,n);
  }
};
