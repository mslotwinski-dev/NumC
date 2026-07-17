#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "log.hpp"

namespace numc {
namespace utility {

/// @brief Reads numerical data from a text file into a 2D vector.
/// @param filename The path to the text file.
/// @return A vector of rows, where each row is a vector of values.
template <typename T = double>
std::vector<std::vector<T>> read_data_file(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    Log::Error("Failed to open file: " + filename);
    throw std::runtime_error("Could not open file " + filename);
  }

  std::vector<std::vector<T>> data;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue; // Skip empty lines and comments
    
    std::stringstream ss(line);
    std::vector<T> row;
    T val;
    while (ss >> val) {
      row.push_back(val);
    }
    if (!row.empty()) {
      data.push_back(row);
    }
  }

  return data;
}

}  // namespace utility
}  // namespace numc
