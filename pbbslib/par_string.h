// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2016 Guy Blelloch and the PBBS team
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

#include <iostream>
#include <fstream>
#include <string>
#include "get_time.h"
#include "sequence_ops.h"

namespace pbbs {

  sequence<char> read_string_from_file(std::string fileName,
				       size_t start, size_t end) {
    std::ifstream file (fileName, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      exit(1);
    }
    size_t length = file.tellg();
    start = std::min(start,length);
    if (end == 0) end = length;
    else end = std::min(end,length);
    size_t n = end - start;
    file.seekg (start, std::ios::beg);
    char* bytes = new_array<char>(n+1);
    //char* bytes = new_array_no_init<char>(n+1);
    file.read (bytes,n);
    file.close();
    return sequence<char>(bytes,n);
  }

  // A structure that keeps a sequence of strings all allocated from
  // the same block of memory
  class words {
  public:
    words() {}

    template<typename Func>
    words(sequence<char>, Func);

    template<typename Func>
    words(std::string filename, Func is_separator, size_t begin=0, size_t end=0);

    ~words() { //if (n > 0) { free(Chars); free(Strings);}
    }

    char* operator[] (const int i) {return Strings[i];}
    size_t size() {return Strings.size();}
    size_t mem_size() {return Str.size() + size()*sizeof(char*);}

  private:
    sequence<char> Str;
    sequence<char*> Strings;
  };

  auto is_space = [] (char c) -> bool {
    switch (c)  {
    case '\r':
    case '\t':
    case '\n':
    case 0:
    case ' ' : return true;
    default : return false;
    }
  };

  auto is_not_alphanum = [] (char c) -> bool {
    return !((c >= 48 && c < 58) ||
	     (c >= 65 && c < 91) ||
	     (c >= 97 && c < 123));
  };

  auto is_not_alpha = [] (char c) -> bool {
    return !((c >= 65 && c < 91) ||
	     (c >= 97 && c < 123));
  };

  // parallel code for converting a string to words
  sequence<char*> string_to_words(sequence<char> Str) {
    size_t n = Str.size();

    // mark start of words
    sequence<bool> FL(n);
    FL[0] = Str[0];
    auto set_f = [&] (size_t i) {
      FL[i] = Str[i] && !Str[i-1];};
    parallel_for(0, n, set_f, 10000);

    auto f = [&] (size_t i) {return &Str[i];};
    auto offset = delayed_seq<char*>(n, f);
    return pack(offset, FL);
  }

  template <typename Func>
  words::words(sequence<char> Str, Func is_separator) : Str(Str) {
    auto set_f = [&] (size_t i) { if (is_separator(Str[i])) Str[i] = 0;};
    parallel_for (0, Str.size(), set_f, 10000);
    Strings = string_to_words(Str);
  }

  template <typename Func>
  words::words(std::string filename, Func is_separator, size_t start, size_t end) {
    sequence<char> S = read_string_from_file(filename, start, end);
    words(S, is_separator);
  }
};

