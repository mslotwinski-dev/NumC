#pragma once

#include <cstdint>
#include <tuple>
#include <stdexcept>

namespace numc {
namespace number_theory {

/// @addtogroup number_theory
/// @{

/// @brief Computes (base^exp) % mod efficiently using binary exponentiation.
inline uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t mod) {
  if (mod == 0) throw std::invalid_argument("Modulo cannot be zero.");
  uint64_t res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp % 2 == 1) {
      // Use 128-bit integer to prevent overflow during multiplication
      res = static_cast<uint64_t>((static_cast<unsigned __int128>(res) * base) % mod);
    }
    base = static_cast<uint64_t>((static_cast<unsigned __int128>(base) * base) % mod);
    exp /= 2;
  }
  return res;
}

/// @brief Computes the greatest common divisor (GCD) of a and b.
inline uint64_t gcd(uint64_t a, uint64_t b) {
  while (b != 0) {
    a %= b;
    std::swap(a, b);
  }
  return a;
}

/// @brief Computes the least common multiple (LCM) of a and b.
inline uint64_t lcm(uint64_t a, uint64_t b) {
  if (a == 0 || b == 0) return 0;
  return (a / gcd(a, b)) * b;
}

/// @brief Extended Euclidean Algorithm.
/// Returns a tuple {gcd(a, b), x, y} such that a*x + b*y = gcd(a, b).
inline std::tuple<int64_t, int64_t, int64_t> extended_gcd(int64_t a, int64_t b) {
  if (a == 0) {
    return {b, 0, 1};
  }
  auto [g, x1, y1] = extended_gcd(b % a, a);
  int64_t x = y1 - (b / a) * x1;
  int64_t y = x1;
  return {g, x, y};
}

/// @brief Computes the modular multiplicative inverse of `a` modulo `m`.
/// Returns `x` such that (a * x) % m == 1. Throws if inverse doesn't exist (i.e. gcd(a, m) != 1).
inline uint64_t mod_inverse(uint64_t a, uint64_t m) {
  auto [g, x, y] = extended_gcd(static_cast<int64_t>(a), static_cast<int64_t>(m));
  if (g != 1) {
    throw std::invalid_argument("Modular inverse does not exist (a and m are not coprime).");
  }
  int64_t m_signed = static_cast<int64_t>(m);
  return static_cast<uint64_t>((x % m_signed + m_signed) % m_signed);
}

/// @brief Chinese Remainder Theorem.
/// Solves the system of congruences: x ≡ a_i (mod m_i).
/// The moduli m_i must be pairwise coprime.
/// Returns the smallest non-negative solution x modulo M (where M = product of m_i).
inline uint64_t chinese_remainder_theorem(const std::vector<uint64_t>& a, const std::vector<uint64_t>& m) {
  if (a.size() != m.size() || a.empty()) {
    throw std::invalid_argument("Arrays a and m must be of the same non-zero size.");
  }

  uint64_t M = 1;
  for (uint64_t mi : m) M *= mi;

  uint64_t result = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    uint64_t Mi = M / m[i];
    uint64_t yi = mod_inverse(Mi, m[i]);
    // Use 128-bit integer for intermediate multiplication
    unsigned __int128 term = static_cast<unsigned __int128>(a[i]) * Mi * yi;
    result = static_cast<uint64_t>((result + term) % M);
  }

  return result;
}

/// @}

}  // namespace number_theory
}  // namespace numc
