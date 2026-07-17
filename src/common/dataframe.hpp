#pragma once

#include "vector.hpp"
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <functional>

namespace numc {

/// @brief DataFrame class for columnar data handling (like Pandas/ROOT RDataFrame).
template <typename T = double>
class DataFrame {
 private:
  std::map<std::string, vector<T>> data_;
  size_t row_count_;

 public:
  DataFrame() : row_count_(0) {}

  /// @brief Loads data from a CSV file.
  /// First row is assumed to be column names.
  static DataFrame<T> from_csv(const std::string& filename, char delimiter = ',') {
    DataFrame<T> df;
    std::ifstream file(filename);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    if (!std::getline(file, line)) return df;

    std::vector<std::string> columns;
    std::stringstream ss(line);
    std::string col_name;
    while (std::getline(ss, col_name, delimiter)) {
      // Trim whitespace from col_name if needed
      columns.push_back(col_name);
      df.data_[col_name] = vector<T>();
    }

    while (std::getline(file, line)) {
      std::stringstream row_ss(line);
      std::string val_str;
      size_t col_idx = 0;
      while (std::getline(row_ss, val_str, delimiter) && col_idx < columns.size()) {
        try {
          T val = static_cast<T>(std::stod(val_str));
          df.data_[columns[col_idx]].push_back(val);
        } catch (...) {
          // If parsing fails, push 0.0 or nan
          df.data_[columns[col_idx]].push_back(T(0.0));
        }
        col_idx++;
      }
      df.row_count_++;
    }
    return df;
  }

  /// @brief Add a new column.
  void add_column(const std::string& name, const vector<T>& column_data) {
    if (row_count_ > 0 && column_data.size() != row_count_) {
      throw std::invalid_argument("Column size mismatch");
    }
    data_[name] = column_data;
    if (row_count_ == 0) {
      row_count_ = column_data.size();
    }
  }

  /// @brief Get a column by name.
  vector<T>& operator[](const std::string& name) {
    return data_.at(name);
  }

  const vector<T>& operator[](const std::string& name) const {
    return data_.at(name);
  }

  size_t get_row_count() const { return row_count_; }

  /// @brief Returns a new DataFrame containing only rows that satisfy the predicate.
  /// The predicate takes a map of column_name -> value for a single row.
  DataFrame<T> filter(const std::function<bool(const std::map<std::string, T>&)>& predicate) const {
    DataFrame<T> result;
    for (const auto& pair : data_) {
      result.data_[pair.first] = vector<T>();
    }

    for (size_t i = 0; i < row_count_; ++i) {
      std::map<std::string, T> row_data;
      for (const auto& pair : data_) {
        row_data[pair.first] = pair.second[i];
      }
      if (predicate(row_data)) {
        for (const auto& pair : data_) {
          result.data_[pair.first].push_back(pair.second[i]);
        }
        result.row_count_++;
      }
    }
    return result;
  }
  /// @brief Returns a new DataFrame containing only rows where mask[i] is true.
  DataFrame<T> filter(const std::vector<bool>& mask) const {
    if (mask.size() != row_count_) {
      throw std::invalid_argument("Mask size does not match DataFrame row count.");
    }
    DataFrame<T> result;
    for (const auto& pair : data_) {
      result.data_[pair.first] = vector<T>();
    }

    for (size_t i = 0; i < row_count_; ++i) {
      if (mask[i]) {
        for (const auto& pair : data_) {
          result.data_[pair.first].push_back(pair.second[i]);
        }
        result.row_count_++;
      }
    }
    return result;
  }
};

}  // namespace numc
