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

#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace gbbs {

struct commandLine {
  int argc;
  char** argv;
  std::string comLine;

  commandLine(int _c, char** _v, std::string _cl)
      : argc(_c), argv(_v), comLine(std::move(_cl)) {}

  commandLine(int _c, char** _v)
      : argc(_c), argv(_v), comLine("bad arguments") {}

  void badArgument() const {
    std::cout << "usage: " << argv[0] << " " << comLine << "\n";
    abort();
  }

  // get an argument
  // i is indexed from the last argument = 0, second to last indexed 1, ..
  char* getArgument(int i) const {
    if (argc < 2 + i) badArgument();
    return argv[argc - 1 - i];
  }

  // looks for two filenames
  std::pair<char*, char*> IOFileNames() const {
    if (argc < 3) badArgument();
    return std::pair<char*, char*>(argv[argc - 2], argv[argc - 1]);
  }

  std::pair<int, char*> sizeAndFileName() const {
    if (argc < 3) badArgument();
    return std::pair<int, char*>(atoi(argv[argc - 2]), (char*)argv[argc - 1]);
  }

  bool getOption(const std::string& option) const {
    for (int i = 1; i < argc; i++)
      if ((std::string)argv[i] == option) return true;
    return false;
  }

  char* getOptionValue(const std::string& option) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) return argv[i + 1];
    return NULL;
  }

  std::string getOptionValue(const std::string& option,
                             const std::string& defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) return (std::string)argv[i + 1];
    return defaultValue;
  }

  int getOptionIntValue(const std::string& option, int defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) {
        int r = atoi(argv[i + 1]);
        return r;
      }
    return defaultValue;
  }

  size_t getOptionLongValue(const std::string& option,
                            size_t defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) {
        long r = atol(argv[i + 1]);
        return r;
      }
    return defaultValue;
  }

  double getOptionDoubleValue(const std::string& option,
                              double defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) {
        double val;
        if (sscanf(argv[i + 1], "%lf", &val) == EOF) {
          badArgument();
        }
        return val;
      }
    return defaultValue;
  }
};

}  // namespace gbbs
