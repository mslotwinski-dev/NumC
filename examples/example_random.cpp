#include "../src/numc.h"
#include <iostream>

using namespace numc;
using namespace numc::statistics;

int main() {
  std::cout << "=== Random, Dist, and Param Example ===" << std::endl;

  // 1. Random numbers
  std::cout << "\n--- Random module ---" << std::endl;
  std::cout << "Uniform [0, 1] (5 values): " << random::uniform(0.0, 1.0, 5) << std::endl;
  std::cout << "Uniform [10, 20] (3 values): " << random::uniform(10.0, 20.0, 3) << std::endl;
  std::cout << "Normal (mu=0, std=1) (4 values): " << random::normal(0.0, 1.0, 4) << std::endl;
  std::cout << "Normal (mu=5, std=2) (4 values): " << random::normal(5.0, 2.0, 4) << std::endl;

  // 2. param and func
  std::cout << "\n--- Dynamic param in func ---" << std::endl;
  param<double> m("m", 2.0);
  param<double> b("b", 1.0);
  
  // Clean declarative mathematical syntax: f(x) = m*x + b
  auto x = funcv::x;
  
  func<double> f = m * x + b;
  std::cout << "f(x) = m*x + b" << std::endl;
  std::cout << "f(5) = " << f(5.0) << " (with m=" << m << ", b=" << b << ")" << std::endl;
  
  m = 3.0; // Instantly changes the behavior of 'f'
  b = 0.0;
  std::cout << "Changed m to 3.0 and b to 0.0" << std::endl;
  std::cout << "f(5) = " << f(5.0) << " (Dynamic update)" << std::endl;

  // 3. dist
  std::cout << "\n--- Dist and Rejection Sampling ---" << std::endl;
  // Let's create a custom distribution: f(x) = x^2 on [0, 2]
  dist<double> my_dist = x^2;
  
  my_dist.range(0.0, 2.0).normalize();
  
  std::cout << "PDF value at x=1 (normalized): " << my_dist(1.0) << std::endl;
  std::cout << "CDF value at x=1 (should be 1/8): " << my_dist.cdf(1.0) << std::endl;
  
  std::cout << "Sampling 5 numbers..." << std::endl;
  auto samples = my_dist.sample(5);
  for (double val : samples) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  return 0;
}
