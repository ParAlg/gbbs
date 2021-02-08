#include <algorithm>
#include <iostream>
#include <fstream>

#include <parlay/io.h>
#include <parlay/primitives.h>
#include <parlay/random.h>

#ifdef USE_PAM
#include <pam/get_time.h>
#include <pam/parse_command_line.h>
#else
#include <cpam/get_time.h>
#include <cpam/parse_command_line.h>
#endif

#include "aspen.h"


int main(int argc, char** argv) {
  std::cout << "Hello!" << std::endl;

  return 0;
}
