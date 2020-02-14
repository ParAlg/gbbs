#include "parse_command_line.h"

#include <cstring>
#include <fstream>
#include <iostream>

commandLine::commandLine(int _c, char** _v, std::string _cl)
  : argc(_c), argv(_v), comLine(std::move(_cl)) {
    if (getOption("-h") || getOption("-help"))
      badArgument();
  }

commandLine::commandLine(int _c, char** _v)
  : argc(_c), argv(_v), comLine("bad arguments") { }

void commandLine::badArgument() {
  std::cout << "usage: " << argv[0] << " " << comLine << std::endl;
  exit(0);
}

char* commandLine::getArgument(int i) {
  if (argc < 2+i) badArgument();
  return argv[argc-1-i];
}

std::pair<char*,char*> commandLine::IOFileNames() {
  if (argc < 3) badArgument();
  return std::pair<char*,char*>(argv[argc-2],argv[argc-1]);
}

std::pair<size_t,char*> commandLine::sizeAndFileName() {
  if (argc < 3) badArgument();
  return std::pair<size_t,char*>(std::atoi(argv[argc-2]),(char*) argv[argc-1]);
}

bool commandLine::getOption(const std::string& option) {
  for (int i = 1; i < argc; i++)
    if ((std::string) argv[i] == option) return true;
  return false;
}

char* commandLine::getOptionValue(const std::string& option) {
  for (int i = 1; i < argc-1; i++)
    if ((std::string) argv[i] == option) return argv[i+1];
  return NULL;
}

std::string commandLine::getOptionValue(
    const std::string& option, const std::string& defaultValue) {
  for (int i = 1; i < argc-1; i++)
    if ((std::string) argv[i] == option) return (std::string) argv[i+1];
  return defaultValue;
}

long
commandLine::getOptionLongValue(const std::string& option, long defaultValue) {
  for (int i = 1; i < argc-1; i++)
    if ((std::string) argv[i] == option) {
      long r = atol(argv[i+1]);
      if (r < 0) badArgument();
      return r;
    }
  return defaultValue;
}

int
commandLine::getOptionIntValue(const std::string& option, int defaultValue) {
  for (int i = 1; i < argc-1; i++)
    if ((std::string) argv[i] == option) {
      int r = atoi(argv[i+1]);
      if (r < 0) badArgument();
      return r;
    }
  return defaultValue;
}

double commandLine::getOptionDoubleValue(
    const std::string& option, double defaultValue) {
  for (int i = 1; i < argc-1; i++)
    if ((std::string) argv[i] == option) {
      double val;
      if (sscanf(argv[i+1], "%lf",  &val) == EOF) {
        badArgument();
      }
      return val;
    }
  return defaultValue;
}
