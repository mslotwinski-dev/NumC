#pragma once

#include "../inc.hpp"
#include "../common/vector.hpp"
#include <vector>
#include <random>
#include <functional>
#include <stdexcept>

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief Generate samples from a 1D probability density function (PDF) using Rejection Sampling.
/// @param pdf The probability density function.
/// @param n_samples The number of samples to generate.
/// @param x_min The minimum x bound of the sampling domain.
/// @param x_max The maximum x bound of the sampling domain.
/// @param pdf_max The maximum value of the PDF in the domain (used as a bounding box).
/// @param seed Random seed for reproducibility.
template <typename T = double, typename F>
std::vector<T> rejection_sample_1d(const F& pdf, size_t n_samples, T x_min, T x_max, T pdf_max, unsigned seed = 42) {
  if (x_min >= x_max || pdf_max <= T(0.0)) {
    throw std::invalid_argument("Invalid sampling boundaries or pdf_max.");
  }
  
  vector<T> samples;
  samples.reserve(n_samples);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<T> dist_x(x_min, x_max);
  std::uniform_real_distribution<T> dist_y(T(0.0), pdf_max);

  while (samples.size() < n_samples) {
    T x_rand = dist_x(rng);
    T y_rand = dist_y(rng);
    if (y_rand <= pdf(x_rand)) {
      samples.push_back(x_rand);
    }
  }

  return samples;
}

/// @brief Generate samples from a 2D probability density function (PDF) using Rejection Sampling.
/// @param pdf The 2D probability density function f(x, y).
/// @param n_samples The number of samples to generate.
/// @param x_min The minimum x bound.
/// @param x_max The maximum x bound.
/// @param y_min The minimum y bound.
/// @param y_max The maximum y bound.
/// @param pdf_max The maximum value of the PDF in the 2D domain.
/// @param seed Random seed.
template <typename T = double, typename F>
std::pair<vector<T>, vector<T>> rejection_sample_2d(const F& pdf, size_t n_samples, 
                                                 T x_min, T x_max, T y_min, T y_max, T pdf_max, unsigned seed = 42) {
  if (x_min >= x_max || y_min >= y_max || pdf_max <= T(0.0)) {
    throw std::invalid_argument("Invalid 2D sampling boundaries or pdf_max.");
  }
  
  vector<T> samples_x;
  vector<T> samples_y;
  samples_x.reserve(n_samples);
  samples_y.reserve(n_samples);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<T> dist_x(x_min, x_max);
  std::uniform_real_distribution<T> dist_y(y_min, y_max);
  std::uniform_real_distribution<T> dist_z(T(0.0), pdf_max);

  while (samples_x.size() < n_samples) {
    T x_rand = dist_x(rng);
    T y_rand = dist_y(rng);
    T z_rand = dist_z(rng);
    if (z_rand <= pdf(x_rand, y_rand)) {
      samples_x.push_back(x_rand);
      samples_y.push_back(y_rand);
    }
  }

  return {samples_x, samples_y};
}

/// @}

}  // namespace statistics
}  // namespace numc
