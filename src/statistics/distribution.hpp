#pragma once

#include "../common/function.hpp"
#include "../common/vector.hpp"
#include "param.hpp"
#include "sampling.hpp"
#include "../analysis/multivariable/optimization.hpp"
#include <cmath>
#include <memory>
#include <stdexcept>

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief Base class for 1D statistical distributions.
/// Inherits from func to allow evaluation of PDF out-of-the-box.
template <typename T = double>
class dist : public func<T> {
 protected:
  T min_val_ = T(0.0);
  T max_val_ = T(1.0);
  T normalization_factor_ = T(1.0);
  T pdf_max_ = T(0.0); // For rejection sampling

 public:
  dist() : func<T>() {}
  
  /// @brief Construct from an existing func<T>
  dist(const func<T>& base_func) : func<T>(base_func) {}

  virtual ~dist() = default;

  /// @brief Set the domain range for the distribution.
  dist& range(T min, T max) {
    min_val_ = min;
    max_val_ = max;
    return *this;
  }

  /// @brief Override evaluation to include the normalization factor and bounds check.
  T operator()(T x) const {
    if (x < min_val_ || x > max_val_) return T(0.0);
    return normalization_factor_ * func<T>::operator()(x);
  }

  /// @brief Find maximum of PDF automatically (heuristic grid search).
  void compute_pdf_max() {
    T max_val = T(0.0);
    T step = (max_val_ - min_val_) / T(1000.0);
    for (T x = min_val_; x <= max_val_; x += step) {
      T val = this->operator()(x);
      if (val > max_val) max_val = val;
    }
    pdf_max_ = max_val;
  }

  /// @brief Normalize the distribution so that its integral over the range is 1.
  dist& normalize() {
    T integral_val = this->integral(min_val_, max_val_);
    if (integral_val > T(0.0)) {
      normalization_factor_ = T(1.0) / integral_val;
    }
    compute_pdf_max();
    return *this;
  }

  /// @brief Cumulative Distribution Function from min_val_ to x.
  virtual T cdf(T x) const {
    if (x <= min_val_) return T(0.0);
    if (x >= max_val_) return T(1.0);
    return this->integral(min_val_, x);
  }

  /// @brief Rejection sampling from this distribution.
  virtual vector<T> sample(size_t n_samples, unsigned seed = 42) {
    if (pdf_max_ <= T(0.0)) {
      compute_pdf_max();
    }
    std::vector<T> raw = rejection_sample_1d<T>(*this, n_samples, min_val_, max_val_, pdf_max_, seed);
    return vector<T>(raw); // Convert to numc::vector
  }

  virtual T mean() const {
    func<T> x_f = func<T>::x();
    func<T> expected = x_f * (*this);
    return expected.integral(min_val_, max_val_);
  }

  virtual T variance() const {
    T m = mean();
    func<T> x_f = func<T>::x();
    func<T> var_f = ((x_f - m) ^ 2.0) * (*this);
    return var_f.integral(min_val_, max_val_);
  }

  virtual T std_dev() const {
    return std::sqrt(variance());
  }
};

/// @brief Base class for N-Dimensional statistical distributions.
template <typename T = double>
class ndist : public func<T> {
 protected:
  vector<T> min_bounds_;
  vector<T> max_bounds_;
  T normalization_factor_ = T(1.0);
  T pdf_max_ = T(0.0);

 public:
  ndist() : func<T>() {}
  ndist(const func<T>& base_func) : func<T>(base_func) {}

  /// @brief Set ND bounds
  ndist& range(const vector<T>& min_b, const vector<T>& max_b) {
    min_bounds_ = min_b;
    max_bounds_ = max_b;
    return *this;
  }

  T operator()(const vector<T>& v) const {
    for (size_t i = 0; i < v.size(); ++i) {
      if (v[i] < min_bounds_[i] || v[i] > max_bounds_[i]) return T(0.0);
    }
    return normalization_factor_ * func<T>::operator()(v);
  }

  T operator()(T x, T y) const { return this->operator()({x, y}); }
  T operator()(T x, T y, T z) const { return this->operator()({x, y, z}); }

  void compute_pdf_max() {
    // For 2D, we can do a simple grid search
    if (min_bounds_.size() == 2) {
      T max_val = T(0.0);
      T step_x = (max_bounds_[0] - min_bounds_[0]) / T(100.0);
      T step_y = (max_bounds_[1] - min_bounds_[1]) / T(100.0);
      for (T x = min_bounds_[0]; x <= max_bounds_[0]; x += step_x) {
        for (T y = min_bounds_[1]; y <= max_bounds_[1]; y += step_y) {
          T val = this->operator()({x, y});
          if (val > max_val) max_val = val;
        }
      }
      pdf_max_ = max_val;
    } else {
      pdf_max_ = T(1.0); // Fallback
    }
  }

  ndist& normalize() {
    if (min_bounds_.size() == 2) {
      T integral_val = this->integral_2d(min_bounds_[0], max_bounds_[0], min_bounds_[1], max_bounds_[1]);
      if (integral_val > T(0.0)) {
        normalization_factor_ = T(1.0) / integral_val;
      }
    } else {
      auto [int_val, err] = this->integral_mc(min_bounds_, max_bounds_, 100000);
      if (int_val > T(0.0)) {
        normalization_factor_ = T(1.0) / int_val;
      }
    }
    compute_pdf_max();
    return *this;
  }

  /// @brief Sample for 2D. Returns a pair of vectors (x_samples, y_samples).
  std::pair<vector<T>, vector<T>> sample_2d(size_t n_samples, unsigned seed = 42) {
    if (min_bounds_.size() != 2) throw std::invalid_argument("Must be 2D to sample_2d");
    if (pdf_max_ <= T(0.0)) compute_pdf_max();
    return rejection_sample_2d<T>(*this, n_samples, min_bounds_[0], max_bounds_[0], min_bounds_[1], max_bounds_[1], pdf_max_, seed);
  }
};

/// @}

}  // namespace statistics
}  // namespace numc
