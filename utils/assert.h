#pragma once

#include <exception>
#include <iostream>

#include <boost/preprocessor/facilities/overload.hpp>

// Asserts on a condition, printing an error and terminating if the condition is
// false.
//
// This is overloaded and can take either one or two arguments.
//
// Arguments
// ---------
//   condition: bool
//     The condition on which to assert.
//   (optional) message: string
//     An explanatory message to print out when the condition fails.
#define ASSERT(...) BOOST_PP_OVERLOAD(__ASSERT, __VA_ARGS__)(__VA_ARGS__)

// Internal. Implementation of ASSERT with one argument.
#define __ASSERT1(condition) \
  do { \
    if (!(condition)) { \
      std::cerr << __FILE__ << ":" << __LINE__ << ": Failed assertion `"  \
          #condition "`" << std::endl; \
      std::terminate(); \
    } \
  } while (false)

// Internal. Implementation of ASSERT with two arguments.
#define __ASSERT2(condition, message) \
  do { \
    if (!(condition)) { \
      std::cerr << __FILE__ << ":" << __LINE__ << ": Failed assertion `"  \
          #condition "`: " << message << std::endl; \
      std::terminate(); \
    } \
  } while (false)
