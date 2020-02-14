#pragma once

#include <sys/time.h>
#include <string>

struct timer {
  double total_time;
  double last_time;
  bool on;
  std::string name;
  struct timezone tzp;

  timer(std::string name = "# PBBS time", bool _start = true);

  double get_time();

  void start();

  double stop();

  void reset();

  double get_total();

  double get_next();

  void report(double time, const std::string& str);

  void total();

  void reportTotal(const std::string& str);

  void next(const std::string& str);
};

static timer _tm;
#define startTime() _tm.start();
#define nextTime(_string) _tm.next(_string);
