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

#include <iostream>
#include <fstream>
#include "../sequence.h"

namespace pbbs {

  // Reads a character sequence from a file
  //    if end is zero or larger than file, then returns full file
  //    if start past end of file then returns an empty string
  sequence<char>
  char_seq_from_file(const std::string& filename, size_t start=0, size_t end=0);

  // Writes a character sequence to a file, returns 0 if successful
  template <class CharSeq>
  int char_seq_to_file(CharSeq const &S, char const *fileName);

  // Returns a sequence of character ranges, one per token
  // The tokens are the longest contiguous subsequences of non space characters.
  // The ranges are over the original sequence, so it should not be deleted
  template <class Seq, class UnaryPred>
  sequence<range<char*>> tokens(Seq const &S, UnaryPred const &is_space);

  // Zeros out all spaces, and returns a pointer to the start of each token
  // Can be used with c style char* functions on each token since they will be null
  // terminated.
  template <class Seq, class UnaryPred>
  sequence<char*> tokenize(Seq &S, UnaryPred const &is_space);

  // Returns a sequence of character ranges, one per partition
  // The StartFlags sequence specifies the start of each partition
  // The two arguments must be of the same length
  // Location 0 is always a start
  template <class Seq, class BoolSeq>
  sequence<range<char*>> partition_at(Seq const &S, BoolSeq const &StartFlags);

  // ********************************
  // Code Follows
  // ********************************

  template <class CharSeq>
  int char_seq_to_file(CharSeq const &S, char const *fileName) {
    size_t n = S.size();
    std::ofstream file (fileName, std::ios::out | std::ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file for writing: " << fileName << std::endl;
      return 1;
    }
    file.write(S.begin(), n);
    file.close();
    return 0;
  }

  // standard definition of a space character
  bool is_space(char c);

  template <class Seq, class UnaryPred>
  sequence<range<char*>> tokens(Seq const &S, UnaryPred const &is_space) {
    size_t n = S.size();

    sequence<bool> StartFlags(n+1);
    sequence<bool> EndFlags(n+1);
    parallel_for(1, n, [&] (long i) {
	bool is = is_space(S[i]);
	StartFlags[i] = !is && is_space(S[i-1]);
	EndFlags[i] = is && !is_space(S[i-1]);
      }, 10000);
    EndFlags[0] = StartFlags[n] = false;
    StartFlags[0] = !is_space(S[0]);
    EndFlags[n] = !is_space(S[n-1]);

    // offset for each start of word
    sequence<long> Starts = pbbs::pack_index<long>(StartFlags);
    sequence<long> Ends = pbbs::pack_index<long>(EndFlags);

    return sequence<range<char*>>(Starts.size(), [&] (size_t i) {
	return range<char*>(S.slice(Starts[i], Ends[i]));});
  }

  template <class Seq, class UnaryPred>
  sequence<char*> tokenize(Seq  &S, UnaryPred const &is_space) {
    size_t n = S.size();

    // clear spaces (side
    parallel_for (0, n, [&] (size_t i) {
	if (is_space(S[i])) S[i] = 0;}, 10000);

    auto StartFlags = delayed_seq<bool>(n, [&] (long i) {
	return (i==0) ? S[i] : S[i] && !S[i-1];});

    // offset for each start of word
    sequence<long> Starts = pbbs::pack_index<long>(StartFlags);

    return sequence<char*>(Starts.size(), [&] (size_t i) {
	return S.begin() + Starts[i];});
  }


  template <class Seq, class BoolSeq>
  sequence<range<char*>> partition_at(Seq const &S, BoolSeq const &StartFlags) {
    size_t n = S.size();
    if (StartFlags.size() != n)
      std::cout << "Unequal sizes in pbbs::partition_at" << std::endl;

    sequence<long> Starts = pbbs::pack_index<long>(StartFlags);
    return sequence<range<char*>>(Starts.size(), [&] (size_t i) {
	long end = (i==Starts.size()-1) ? n : Starts[i+1];
	return range<char*>(S.slice(Starts[i],end));});
  }

}  // namespace pbbs

