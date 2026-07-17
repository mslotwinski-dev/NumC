#pragma once

#include "../inc.hpp"
#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include <random>

namespace numc {

/// @defgroup random Random Number Generation
/// @ingroup statistics
/// @brief High-level API for generating numc::vectors filled with random numbers.
/// @{

namespace random {

  /// @brief Generates a vector filled with uniformly distributed random numbers.
  /// @tparam T Floating point type.
  /// @param min Minimum value.
  /// @param max Maximum value.
  /// @param size Size of the generated vector.
  /// @param seed Random seed.
  /// @return A numc::vector filled with uniformly distributed values.
  template <typename T = double>
  vector<T> uniform(T min, T max, size_t size, unsigned seed = std::random_device{}()) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<T> dist(min, max);
    vector<T> res(size);
    for (size_t i = 0; i < size; ++i) {
      res[i] = dist(rng);
    }
    return res;
  }

  /// @brief Generates a tensor filled with uniformly distributed random numbers.
  /// @tparam T Floating point type.
  /// @param rows Number of rows.
  /// @param cols Number of columns.
  /// @param min Minimum value.
  /// @param max Maximum value.
  /// @param seed Random seed.
  /// @return A numc::tensor filled with uniformly distributed values.
  template <typename T = double>
  tensor<T> uniform_tensor(size_t rows, size_t cols, T min, T max, unsigned seed = std::random_device{}()) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<T> dist(min, max);
    tensor<T> res({rows, cols});
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        res(i, j) = dist(rng);
      }
    }
    return res;
  }

  /// @brief Generates a vector filled with normally (Gaussian) distributed random numbers.
  /// @tparam T Floating point type.
  /// @param mean Mean of the distribution.
  /// @param stddev Standard deviation of the distribution.
  /// @param size Size of the generated vector.
  /// @param seed Random seed.
  /// @return A numc::vector filled with normally distributed values.
  template <typename T = double>
  vector<T> normal(T mean, T stddev, size_t size, unsigned seed = std::random_device{}()) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<T> dist(mean, stddev);
    vector<T> res(size);
    for (size_t i = 0; i < size; ++i) {
      res[i] = dist(rng);
    }
    return res;
  }

} // namespace random

/// @}

} // namespace numc
