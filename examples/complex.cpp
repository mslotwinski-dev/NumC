#include <iostream>
#include <numc>  // NumC library

using namespace numc;

int main() {
  auto i = complex<>::i();  // == complex<double>::i().

  complex z1 = 3 + 4 * i;
  complex z2 = 1 + 2 * i;

  complex z3 = z1 + z2;         // Addition
  complex z4 = ~z3;             // Conjugate
  complex z5 = z1 ^ z2;         // Exponentiation
  complex z6 = (E ^ (PI * i));  // Euler's formula

  std::cout << "Complex number created! z1 = " << z1 << std::endl;  // Complex number created! z1 = 3 + 4i
  std::cout << "Complex number created! z2 = " << z2 << std::endl;  // Complex number created! z2 = 1 + 2i
  std::cout << "Complex number created! z3 = " << z3 << std::endl;  // Complex number created! z3 = 4 + 6i
  std::cout << "Complex number created! z4 = " << z4 << std::endl;  // Complex number created! z4 = 4 - 6i
  std::cout << "Complex number created! z5 = " << z5 << std::endl;  // Complex number created! z5 = -0.419813 - 0.660452i
  std::cout << "Complex number created! z6 = " << z6 << std::endl;  // Complex number created! z6 = -1

  return 0;
}
