#include "parse_command_line.h"

commandLine::commandLine(int _c, char** _v, std::string _cl)
    : argc(_c), argv(_v), comLine(std::move(_cl)) {}

commandLine::commandLine(int _c, char** _v)
    : argc(_c), argv(_v), comLine("bad arguments") {}

void commandLine::badArgument() {
  std::cout << "usage: " << argv[0] << " " << comLine << "\n";
  abort();
}

// get an argument
// i is indexed from the last argument = 0, second to last indexed 1, ..
char* commandLine::getArgument(int i) {
  if (argc < 2 + i) badArgument();
  return argv[argc - 1 - i];
}

// looks for two filenames
std::pair<char*, char*> commandLine::IOFileNames() {
  if (argc < 3) badArgument();
  return std::pair<char*, char*>(argv[argc - 2], argv[argc - 1]);
}

std::pair<int, char*> commandLine::sizeAndFileName() {
  if (argc < 3) badArgument();
  return std::pair<int, char*>(atoi(argv[argc - 2]), (char*)argv[argc - 1]);
}

bool commandLine::getOption(const std::string& option) {
  for (int i = 1; i < argc; i++)
    if ((std::string)argv[i] == option) return true;
  return false;
}

char* commandLine::getOptionValue(const std::string& option) {
  for (int i = 1; i < argc - 1; i++)
    if ((std::string)argv[i] == option) return argv[i + 1];
  return NULL;
}

std::string commandLine::getOptionValue(
    const std::string& option, const std::string& defaultValue) {
  for (int i = 1; i < argc - 1; i++)
    if ((std::string)argv[i] == option) return (std::string)argv[i + 1];
  return defaultValue;
}

int commandLine::getOptionIntValue(
    const std::string& option, int defaultValue) {
  for (int i = 1; i < argc - 1; i++)
    if ((std::string)argv[i] == option) {
      int r = atoi(argv[i + 1]);
      return r;
    }
  return defaultValue;
}

size_t commandLine::getOptionLongValue(
    const std::string& option, size_t defaultValue) {
  for (int i = 1; i < argc - 1; i++)
    if ((std::string)argv[i] == option) {
      long r = atol(argv[i + 1]);
      return r;
    }
  return defaultValue;
}

double commandLine::getOptionDoubleValue(
    const std::string& option, double defaultValue) {
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
