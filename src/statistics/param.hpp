#pragma once

#include <string>
#include <iostream>

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief A parameter class for statistical distributions or curve fitting.
/// It wraps a value and provides additional context such as name and uncertainty.
template <typename T = double>
class param {
 private:
  std::string name_;
  T value_;
  T error_;

 public:
  param(T val = T(0.0)) : name_(""), value_(val), error_(T(0.0)) {}
  param(const std::string& name, T val = T(0.0), T error = T(0.0))
      : name_(name), value_(val), error_(error) {}

  /// @brief Implicit conversion to the underlying numeric type.
  operator T() const { return value_; }

  /// @brief Assign a new value.
  param& operator=(T val) {
    value_ = val;
    return *this;
  }

  T value() const { return value_; }
  T error() const { return error_; }
  std::string name() const { return name_; }

  void set_value(T val) { value_ = val; }
  void set_error(T err) { error_ = err; }
  void set_name(const std::string& name) { name_ = name; }

  friend std::ostream& operator<<(std::ostream& os, const param& p) {
    if (!p.name_.empty()) {
      os << p.name_ << " = ";
    }
    os << p.value_;
    if (p.error_ > T(0.0)) {
      os << " +/- " << p.error_;
    }
    return os;
  }
};

/// @}

}  // namespace statistics
}  // namespace numc
