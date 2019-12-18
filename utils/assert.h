// Copyright 2019 Thomas Tseng
#pragma once

#include <exception>
#include <iostream>

#define ASSERT(condition, message) \
  do { \
    if (!(condition)) { \
      std::cerr << __FILE__ << ":" << __LINE__ << ": Failed assertion `"  \
          #condition "`: " << message << std::endl; \
      std::terminate(); \
    } \
  } while (false)
