#include <iostream>
#include <numc>

using namespace numc;

int main() {
  vector v1 = {1.0, 2.0, 3.0};
  vector v2 = {4.0, 5.0, 6.0};
  vector v3(3);

  // vector v3 = v1 + v2;   // Vector addition
  // vector v4 = v1 * 2.0;  // Scalar multiplication

  std::cout << "Vector created! v1 = " << v1 << std::endl;  // Vector created! v1 = (1, 2, 3)
  std::cout << "Vector created! v2 = " << v2 << std::endl;  // Vector created! v2 = (4, 5, 6)
  std::cout << "Vector created! v3 = " << v3 << std::endl;  // Vector created! v3 = (0, 0, 0)
  // std::cout << "Vector created! v4 = " << v4 << std::endl;  // Vector created! v4 = (2, 4, 6)

  return 0;
}