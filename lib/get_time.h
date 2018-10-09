// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and the PBBS team
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
#include <sys/time.h>
#include <iomanip>
#include <iostream>
#include <string>

struct timer {
  double total_time;
  double last_time;
  bool on;
  std::string name;
  struct timezone tzp;

  timer(std::string name = "PBBS time", bool _start = true)
      : total_time(0.0), on(false), name(name), tzp({0, 0}) {
    if (_start) start();
  }

  double get_time() {
    timeval now;
    gettimeofday(&now, &tzp);
    return ((double)now.tv_sec) + ((double)now.tv_usec) / 1000000.;
  }

  void start() {
    on = 1;
    last_time = get_time();
  }

  double stop() {
    on = 0;
    double d = (get_time() - last_time);
    total_time += d;
    return d;
  }

  void reset() {
    total_time = 0.0;
    on = 0;
  }

  double get_total() {
    if (on)
      return total_time + get_time() - last_time;
    else
      return total_time;
  }

  double get_next() {
    if (!on) return 0.0;
    double t = get_time();
    double td = t - last_time;
    total_time += td;
    last_time = t;
    return td;
  }

  void report(double time, std::string str) {
    std::cout << name << ": " << str << ": " << std::setprecision(3) << time
              << "\n";
  }

  void total() {
    report(get_total(), "total");
    total_time = 0.0;
  }

  void reportTotal(std::string str) {
    report(get_total(), str);
    total_time = 0.0;
  }

  void next(std::string str) {
    if (on) report(get_next(), str);
  }
};
