#pragma once

#include "../common/function.hpp"
#include "../common/vector.hpp"
#include "param.hpp"
#include <cmath>
#include <memory>

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief Base class for statistical distributions.
/// Inherits from func to allow evaluation of PDF out-of-the-box.
template <typename T = double>
class dist : public func<T> {
 public:
  /// @brief Construct a distribution with a specific PDF.
  dist(std::function<T(T)> pdf) : func<T>(pdf) {}

  virtual ~dist() = default;

  /// @brief Cumulative Distribution Function.
  /// @param x Point at which to evaluate the CDF.
  virtual T cdf(T x) const = 0;

  /// @brief Mean of the distribution.
  virtual T mean() const = 0;

  /// @brief Variance of the distribution.
  virtual T variance() const = 0;

  /// @brief Standard deviation of the distribution.
  virtual T std_dev() const {
    return std::sqrt(variance());
  }
};

/// @brief Normal (Gaussian) distribution.
template <typename T = double>
class normal_dist : public dist<T> {
 private:
  param<T> mu_;
  param<T> sigma_;

 public:
  normal_dist(T mean = T(0.0), T std = T(1.0))
      : dist<T>([mean, std](T x) {
          T diff = x - mean;
          return std::exp(-diff * diff / (T(2.0) * std * std)) / (std * std::sqrt(T(2.0) * T(3.14159265358979323846)));
        }),
        mu_("mu", mean),
        sigma_("sigma", std) {}

  T cdf(T x) const override {
    return T(0.5) * (T(1.0) + std::erf((x - mu_.value()) / (sigma_.value() * std::sqrt(T(2.0)))));
  }

  T mean() const override { return mu_.value(); }
  T variance() const override { return sigma_.value() * sigma_.value(); }
};

/// @brief Uniform distribution.
template <typename T = double>
class uniform_dist : public dist<T> {
 private:
  param<T> a_;
  param<T> b_;

 public:
  uniform_dist(T a = T(0.0), T b = T(1.0))
      : dist<T>([a, b](T x) {
          if (x >= a && x <= b) return T(1.0) / (b - a);
          return T(0.0);
        }),
        a_("a", a),
        b_("b", b) {}

  T cdf(T x) const override {
    if (x < a_.value()) return T(0.0);
    if (x > b_.value()) return T(1.0);
    return (x - a_.value()) / (b_.value() - a_.value());
  }

  T mean() const override { return (a_.value() + b_.value()) / T(2.0); }
  T variance() const override { 
    T diff = b_.value() - a_.value();
    return (diff * diff) / T(12.0); 
  }
};

/// @brief Exponential distribution.
template <typename T = double>
class exponential_dist : public dist<T> {
 private:
  param<T> lambda_;

 public:
  exponential_dist(T lambda = T(1.0))
      : dist<T>([lambda](T x) {
          if (x >= T(0.0)) return lambda * std::exp(-lambda * x);
          return T(0.0);
        }),
        lambda_("lambda", lambda) {}

  T cdf(T x) const override {
    if (x < T(0.0)) return T(0.0);
    return T(1.0) - std::exp(-lambda_.value() * x);
  }

  T mean() const override { return T(1.0) / lambda_.value(); }
  T variance() const override { return T(1.0) / (lambda_.value() * lambda_.value()); }
};

/// @}

}  // namespace statistics
}  // namespace numc
