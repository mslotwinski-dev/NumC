#pragma once

#include "../common/complex.hpp"
#include "../common/vector.hpp"
#include "../common/function.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Computes the Discrete Fourier Transform (DFT) of a real-valued signal.
template <typename T = double>
std::vector<complex<T>> dft(const vector<T>& x) {
  size_t N = x.size();
  std::vector<complex<T>> X(N);
  T two_pi = T(2.0) * T(numc::PI);
  for (size_t k = 0; k < N; ++k) {
    T re = T(0.0), im = T(0.0);
    for (size_t n = 0; n < N; ++n) {
      T angle = two_pi * T(k) * T(n) / T(N);
      re += x[n] * std::cos(angle);
      im -= x[n] * std::sin(angle);
    }
    X[k] = complex<T>(re, im);
  }
  return X;
}

/// @brief Computes the Inverse Discrete Fourier Transform (IDFT).
template <typename T = double>
vector<T> idft(const std::vector<complex<T>>& X) {
  size_t N = X.size();
  vector<T> x(N);
  T two_pi = T(2.0) * T(numc::PI);
  for (size_t n = 0; n < N; ++n) {
    T re = T(0.0);
    for (size_t k = 0; k < N; ++k) {
      T angle = two_pi * T(k) * T(n) / T(N);
      re += X[k].real() * std::cos(angle) - X[k].imag() * std::sin(angle);
    }
    x[n] = re / T(N);
  }
  return x;
}

/// @brief In-place Cooley-Tukey radix-2 FFT.
template <typename T = double>
void fft_inplace(std::vector<complex<T>>& X) {
  size_t N = X.size();
  if (N <= 1) return;

  // Bit-reversal permutation
  size_t j = 0;
  for (size_t i = 1; i < N; ++i) {
    size_t bit = N >> 1;
    while (j & bit) { j ^= bit; bit >>= 1; }
    j ^= bit;
    if (i < j) std::swap(X[i], X[j]);
  }

  // Butterfly stages
  for (size_t len = 2; len <= N; len <<= 1) {
    T angle = -T(2.0) * T(numc::PI) / T(len);
    complex<T> wlen(std::cos(angle), std::sin(angle));
    for (size_t i = 0; i < N; i += len) {
      complex<T> w(T(1.0), T(0.0));
      for (size_t jj = 0; jj < len / 2; ++jj) {
        complex<T> u = X[i + jj];
        complex<T> v = w * X[i + jj + len / 2];
        X[i + jj] = u + v;
        X[i + jj + len / 2] = u - v;
        w = w * wlen;
      }
    }
  }
}

/// @brief Computes the FFT of a real-valued signal. Zero-pads to the next power of 2.
template <typename T = double>
std::vector<complex<T>> fft(const vector<T>& x) {
  size_t N = x.size();
  // Find next power of 2
  size_t M = 1;
  while (M < N) M <<= 1;
  std::vector<complex<T>> X(M);
  for (size_t i = 0; i < N; ++i) X[i] = complex<T>(x[i], T(0.0));
  for (size_t i = N; i < M; ++i) X[i] = complex<T>(T(0.0), T(0.0));
  fft_inplace(X);
  return X;
}

/// @brief In-place Inverse FFT.
template <typename T = double>
void ifft_inplace(std::vector<complex<T>>& X) {
  size_t N = X.size();
  // Conjugate
  for (size_t i = 0; i < N; ++i) {
    X[i] = complex<T>(X[i].real(), -X[i].imag());
  }
  fft_inplace(X);
  // Conjugate and scale
  for (size_t i = 0; i < N; ++i) {
    X[i] = complex<T>(X[i].real() / T(N), -X[i].imag() / T(N));
  }
}

/// @brief Computes the IFFT and returns the real part as a vector.
template <typename T = double>
vector<T> ifft(std::vector<complex<T>> X) {
  ifft_inplace(X);
  vector<T> result(X.size());
  for (size_t i = 0; i < X.size(); ++i) result[i] = X[i].real();
  return result;
}

/// @brief Computes the power spectral density |X(k)|^2 / N.
template <typename T = double>
vector<T> power_spectrum(const vector<T>& x) {
  auto X = fft(x);
  size_t N = X.size();
  vector<T> P(N);
  for (size_t k = 0; k < N; ++k) {
    P[k] = (X[k].real() * X[k].real() + X[k].imag() * X[k].imag()) / T(N);
  }
  return P;
}

/// @brief Computes the frequency bins for an FFT of size N with sample rate fs.
template <typename T = double>
vector<T> fft_freqs(size_t N, T fs = T(1.0)) {
  vector<T> freqs(N);
  for (size_t i = 0; i <= N / 2; ++i) {
    freqs[i] = T(i) * fs / T(N);
  }
  for (size_t i = N / 2 + 1; i < N; ++i) {
    freqs[i] = (T(i) - T(N)) * fs / T(N);
  }
  return freqs;
}

/// @brief Convolves two real-valued signals using FFT.
template <typename T = double>
vector<T> convolve(const vector<T>& a, const vector<T>& b) {
  size_t outLen = a.size() + b.size() - 1;
  size_t M = 1;
  while (M < outLen) M <<= 1;

  std::vector<complex<T>> A(M), B(M);
  for (size_t i = 0; i < a.size(); ++i) A[i] = complex<T>(a[i], T(0.0));
  for (size_t i = 0; i < b.size(); ++i) B[i] = complex<T>(b[i], T(0.0));

  fft_inplace(A);
  fft_inplace(B);

  std::vector<complex<T>> C(M);
  for (size_t i = 0; i < M; ++i) C[i] = A[i] * B[i];

  ifft_inplace(C);

  vector<T> result(outLen);
  for (size_t i = 0; i < outLen; ++i) result[i] = C[i].real();
  return result;
}

/// @brief Computes the cross-correlation of two signals.
template <typename T = double>
vector<T> correlate(const vector<T>& a, const vector<T>& b) {
  size_t n = a.size(), m = b.size();
  size_t outLen = n + m - 1;
  // Reverse b and convolve
  vector<T> b_rev(m);
  for (size_t i = 0; i < m; ++i) b_rev[i] = b[m - 1 - i];
  return convolve(a, b_rev);
}

/// @}

}  // namespace analysis
}  // namespace numc
