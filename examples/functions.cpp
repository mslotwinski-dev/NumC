#include <iostream>
#include <numc>

using namespace numc;

int main() {
  auto x = func<>::x();  // Defualt funkcja f(x) = x na typie double

  func f = x ^ 2;               // Definiowanie funkcji f(x) = x^x
  func g = 2 * sin(x / 2) + x;  // Definiowanie funkcji g(x) = sin(x)

  plt::show({f, g}, -10.0, 10.0, "Funkcje", "t [s]", "A [cm]");

  // func g = [](double x) { return std::sin(x); };
  // Postaram się, żeby definiowanie ręczne było rzadko potrzebne.

  func h = [](double x) {
    if (x <= 0) return 2137.0;
    if (x > 0) return 420.0 * x;
    return 0.0;
  };

  plt::show(h, -10, 10, "Funkcja", "t [s]", "A [cm]");

  double y = f(2);

  std::cout << "f(2) = " << y << std::endl;                                                   // f(2) = 4
  std::cout << "g(pi) = " << g(PI) << std::endl;                                              // g(pi) = 1.22461e-16 (mało)
  std::cout << "(f + g)(1) = " << (f + g)(1) << std::endl;                                    // f(1) = 1.84147
  std::cout << "f(1) = " << f(1) << ", f(2) = " << f(2) << ", f(3) = " << f(3) << std::endl;  // f(1) = 1, f(2) = 4, f(3) = 27

  return 0;
}
