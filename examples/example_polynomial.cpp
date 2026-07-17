#include "../src/numc.h"
#include <iostream>

using namespace numc;

int main() {
  std::cout << "=== Polynomial Example ===" << std::endl;

  // p(x) = 1.0 + 2.0*x + 3.0*x^2
  polynomial<double> p1({1.0, 2.0, 3.0});
  std::cout << "p1(x) = " << p1 << std::endl;

  // p2(x) = -1.0 + 5.0*x
  polynomial<double> p2{-1.0, 5.0};
  std::cout << "p2(x) = " << p2 << std::endl;

  polynomial<double> sum = p1 + p2;
  std::cout << "p1(x) + p2(x) = " << sum << std::endl;

  polynomial<double> prod = p1 * p2;
  std::cout << "p1(x) * p2(x) = " << prod << std::endl;

  polynomial<double> deriv = p1.derivative();
  std::cout << "p1'(x) = " << deriv << std::endl;

  polynomial<double> integ = p1.antiderivative(0.0);
  std::cout << "Integral of p1(x) = " << integ << std::endl;

  double eval = p1(2.0);
  std::cout << "p1(2.0) = " << eval << std::endl;

  // Convert to func for analysis features
  func<double> f = p1.to_func();
  std::cout << "p1 evaluated as func at 2.0 = " << f(2.0) << std::endl;

  return 0;
}
