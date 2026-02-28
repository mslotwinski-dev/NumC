#include <iostream>
#include <numc>

using namespace numc;

int main() {
  vector v1 = {1.0, 2.0, 3.0};
  vector v2 = {4.0, 5.0, 6.0};

  vector v3 = v1 + v2;   // Vector addition
  vector v4 = v1 * 2.0;  // Scalar multiplication
  double x5 = v1 * v2;   // Dot product
  vector v6 = v1 ^ v2;   // Cross product

  std::cout << "Vector created! v1 = " << v1 << std::endl;  // Vector created! v1 = (1, 2, 3)
  std::cout << "Vector created! v2 = " << v2 << std::endl;  // Vector created! v2 = (4, 5, 6)
  std::cout << "Vector created! v3 = " << v3 << std::endl;  // Vector created! v3 = (5, 7, 9)
  std::cout << "Vector created! v4 = " << v4 << std::endl;  // Vector created! v4 = (2, 4, 6)
  std::cout << "Scalar created! x5 = " << x5 << std::endl;  // Vector created! x5 = 32
  std::cout << "Vector created! v6 = " << v6 << std::endl;  // Vector created! v6 = (-3, 6, -3)

  return 0;
}