#pragma once

#include <string>
#include <iostream>
#include <memory>
#include "../common/function.hpp"

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief A parameter class for statistical distributions or curve fitting.
/// It wraps a value inside a shared_ptr, allowing its value to be dynamically
/// evaluated within any mathematical expression (func).
template <typename T = double>
class param : public func<T> {
 private:
  std::string name_;
  std::shared_ptr<T> value_ptr_;
  T error_;

 public:
  param(T val = T(0.0)) 
      : func<T>(), name_(""), value_ptr_(std::make_shared<T>(val)), error_(T(0.0)) {
    // Override the base func to return the dynamic shared value
    auto ptr = value_ptr_;
    this->f = [ptr](const vector<T>&) { return *ptr; };
  }

  param(const std::string& name, T val = T(0.0), T error = T(0.0))
      : func<T>(), name_(name), value_ptr_(std::make_shared<T>(val)), error_(error) {
    auto ptr = value_ptr_;
    this->f = [ptr](const vector<T>&) { return *ptr; };
  }

  /// @brief Implicit conversion to the underlying numeric type.
  operator T() const { return *value_ptr_; }

  /// @brief Assign a new value.
  param& operator=(T val) {
    *value_ptr_ = val;
    return *this;
  }

  T value() const { return *value_ptr_; }
  T error() const { return error_; }
  std::string name() const { return name_; }

  void set_value(T val) { *value_ptr_ = val; }
  void set_error(T err) { error_ = err; }
  void set_name(const std::string& name) { name_ = name; }

  friend std::ostream& operator<<(std::ostream& os, const param& p) {
    if (!p.name_.empty()) {
      os << p.name_ << " = ";
    }
    os << *(p.value_ptr_);
    if (p.error_ > T(0.0)) {
      os << " +/- " << p.error_;
    }
    return os;
  }
};

/// @}

}  // namespace statistics
}  // namespace numc
