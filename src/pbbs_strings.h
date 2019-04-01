// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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

#include "pbbslib/strings/string_basics.h"

namespace pbbslib {

  // Reads a character sequence from a file
  //    if end is zero or larger than file, then returns full file
  //    if start past end of file then returns an empty string
  inline sequence<char> char_seq_from_file(std::string filename, size_t start=0, size_t end=0) {
    return pbbs::char_seq_from_file(filename, start, end);
  }

  // Writes a character sequence to a file, returns 0 if successful
  template <class CharSeq>
  int char_seq_to_file(CharSeq S, char const *fileName) {
    return pbbs::char_seq_to_file(S, fileName);
  }

  // Returns a sequence of character ranges, one per token
  // The tokens are the longest contiguous subsequences of non space characters.
  // The ranges are over the original sequence, so it should not be deleted
  template <class Seq, class UnaryPred>
  sequence<range<char*>> tokens(Seq const &S, UnaryPred const &is_space) {
    return pbbs::tokens(S, is_space);
  }

  // Zeros out all spaces, and returns a pointer to the start of each token
  // Can be used with c style char* functions on each token since they will be null
  // terminated.
  template <class Seq, class UnaryPred>
  sequence<char*> tokenize(Seq &S, UnaryPred const &is_space) {
    return pbbs::tokenize(S, is_space);
  }

  // Returns a sequence of character ranges, one per partition
  // The StartFlags sequence specifies the start of each partition
  // The two arguments must be of the same length
  // Location 0 is always a start
  template <class Seq, class BoolSeq>
  sequence<range<char*>> partition_at(Seq const &S, BoolSeq const &StartFlags) {
    return pbbs::partition_at(S, StartFlags);
  }

}
