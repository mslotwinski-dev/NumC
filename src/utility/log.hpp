#pragma once
#include "../inc.hpp"

enum Colors { NONE, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, ORANGE, LIGHT };

class Log {
 public:
  Log() {}
  ~Log() {}

  static void Info(std::string message);
  static void Warn(std::string message);
  static void Error(std::string message);

  static std::string Color(Colors);
};