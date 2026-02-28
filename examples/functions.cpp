#include <iostream>
#include <numc>

using namespace numc;

int main() {
  auto x = func<>::x();  // Defualt funkcja f(x) = x na typie double

  func f = x ^ x;   // Definiowanie funkcji f(x) = x^x
  func g = sin(x);  // Definiowanie funkcji g(x) = sin(x)

  // func g = [](double x) { return std::sin(x); };
  // Postaram się, żeby definiowanie ręczne było rzadko potrzebne.

  double y = f(2);

  std::cout << "f(2) = " << y << std::endl;                                                   // f(2) = 4
  std::cout << "g(pi) = " << g(PI) << std::endl;                                              // g(pi) = 1.22461e-16 (mało)
  std::cout << "(f + g)(1) = " << (f + g)(1) << std::endl;                                    // f(1) = 1.84147
  std::cout << "f(1) = " << f(1) << ", f(2) = " << f(2) << ", f(3) = " << f(3) << std::endl;  // f(1) = 1, f(2) = 4, f(3) = 27

  return 0;
}
