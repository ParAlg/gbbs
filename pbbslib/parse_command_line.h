// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#include <string>

struct commandLine {
  int argc;
  char** argv;
  std::string comLine;
  commandLine(int _c, char** _v, std::string _cl);

  commandLine(int _c, char** _v);

  void badArgument();

  // get an argument
  // i is indexed from the last argument = 0, second to last indexed 1, ..
  char* getArgument(int i);

  // looks for two filenames
  std::pair<char*,char*> IOFileNames();

  std::pair<size_t,char*> sizeAndFileName();

  bool getOption(const std::string& option);

  char* getOptionValue(const std::string& option);

  std::string
  getOptionValue(const std::string& option, const std::string& defaultValue);

  long getOptionLongValue(const std::string& option, long defaultValue);

  int getOptionIntValue(const std::string& option, int defaultValue);

  double getOptionDoubleValue(const std::string& option, double defaultValue);
};
