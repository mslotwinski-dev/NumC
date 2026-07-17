#pragma once

#include "modular.hpp"
#include <vector>
#include <random>

namespace numc {
namespace number_theory {

/// @addtogroup number_theory
/// @{

/// @brief Generates all prime numbers up to a given limit using the Sieve of Eratosthenes.
inline std::vector<uint64_t> sieve_of_eratosthenes(uint64_t limit) {
  if (limit < 2) return {};
  std::vector<bool> is_prime(limit + 1, true);
  is_prime[0] = is_prime[1] = false;

  for (uint64_t p = 2; p * p <= limit; ++p) {
    if (is_prime[p]) {
      for (uint64_t i = p * p; i <= limit; i += p)
        is_prime[i] = false;
    }
  }

  std::vector<uint64_t> primes;
  for (uint64_t p = 2; p <= limit; ++p) {
    if (is_prime[p]) primes.push_back(p);
  }
  return primes;
}

/// @brief Miller-Rabin primality test. Deterministic for n < 2^64.
/// @param n Number to test.
/// @return True if n is prime, false otherwise.
inline bool is_prime_miller_rabin(uint64_t n) {
  if (n < 2) return false;
  if (n == 2 || n == 3) return true;
  if (n % 2 == 0) return false;

  uint64_t d = n - 1;
  int s = 0;
  while (d % 2 == 0) {
    d /= 2;
    s++;
  }

  // Bases that provide deterministic results for n < 2^64
  static const uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

  for (uint64_t a : bases) {
    if (n <= a) break;
    
    uint64_t x = mod_pow(a, d, n);
    if (x == 1 || x == n - 1) continue;

    bool composite = true;
    for (int r = 1; r < s; ++r) {
      x = static_cast<uint64_t>((static_cast<unsigned __int128>(x) * x) % n);
      if (x == n - 1) {
        composite = false;
        break;
      }
    }

    if (composite) return false;
  }
  return true;
}

/// @brief Helper function for Pollard's rho algorithm.
inline uint64_t pollard_rho(uint64_t n) {
  if (n % 2 == 0) return 2;
  if (is_prime_miller_rabin(n)) return n;

  std::mt19937_64 rng(1337);
  std::uniform_int_distribution<uint64_t> dist(1, n - 1);
  
  uint64_t x = dist(rng);
  uint64_t c = dist(rng);
  uint64_t y = x;
  uint64_t g = 1;

  auto f = [&](uint64_t val) {
    return static_cast<uint64_t>((static_cast<unsigned __int128>(val) * val + c) % n);
  };

  while (g == 1) {
    x = f(x);
    y = f(f(y));
    uint64_t diff = x > y ? x - y : y - x;
    g = gcd(diff, n);
  }

  if (g == n) return pollard_rho(n); // Failed, retry (recursion is fine with random c)
  return g;
}

/// @brief Computes the prime factorization of a number using Pollard's rho algorithm.
/// Returns a sorted list of prime factors.
inline std::vector<uint64_t> prime_factors(uint64_t n) {
  if (n <= 1) return {};
  if (is_prime_miller_rabin(n)) return {n};

  uint64_t divisor = pollard_rho(n);
  std::vector<uint64_t> factors1 = prime_factors(divisor);
  std::vector<uint64_t> factors2 = prime_factors(n / divisor);
  
  factors1.insert(factors1.end(), factors2.begin(), factors2.end());
  std::sort(factors1.begin(), factors1.end());
  return factors1;
}

/// @brief Euler's totient function φ(n).
/// Computes the number of positive integers up to n that are relatively prime to n.
inline uint64_t euler_totient(uint64_t n) {
  if (n == 0) return 0;
  uint64_t result = n;
  std::vector<uint64_t> factors = prime_factors(n);
  
  // Unique factors
  if (!factors.empty()) {
    uint64_t last = 0;
    for (uint64_t p : factors) {
      if (p != last) {
        result -= result / p;
        last = p;
      }
    }
  }
  return result;
}

/// @}

}  // namespace number_theory
}  // namespace numc
