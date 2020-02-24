#include "get_time.h"

#include <iomanip>
#include <iostream>

timer::timer(std::string _name, bool _start)
    : total_time(0.0), on(false), name(std::move(_name)), tzp({0, 0}) {
  if (_start) start();
}

double timer::get_time() {
  timeval now;
  gettimeofday(&now, &tzp);
  return ((double)now.tv_sec) + ((double)now.tv_usec) / 1000000.;
}

void timer::start() {
  on = 1;
  last_time = get_time();
}

double timer::stop() {
  on = 0;
  double d = (get_time() - last_time);
  total_time += d;
  return d;
}

void timer::reset() {
  total_time = 0.0;
  on = 0;
}

double timer::get_total() {
  if (on)
    return total_time + get_time() - last_time;
  else
    return total_time;
}

double timer::get_next() {
  if (!on) return 0.0;
  double t = get_time();
  double td = t - last_time;
  total_time += td;
  last_time = t;
  return td;
}

void timer::report(double time, const std::string& str) {
  std::ios::fmtflags cout_settings = std::cout.flags();
  std::cout.precision(4);
  std::cout << std::fixed;
  std::cout << name << ": ";
  if (str.length() > 0) std::cout << str << ": ";
  std::cout << time << std::endl;
  std::cout.flags(cout_settings);
}

void timer::total() {
  report(get_total(), "total");
  total_time = 0.0;
}

void timer::reportTotal(const std::string& str) { report(get_total(), str); }

void timer::next(const std::string& str) {
  if (on) report(get_next(), str);
}
