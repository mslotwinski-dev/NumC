#include "../src/numc.h"
#include <iostream>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 9 Replication: Szukanie Pierwiastków ==========\n";

  auto x = func<>::x();
  
  // Funkcja 1: f(x) = exp(-x) * sin(5*x)
  func f1 = exp(-x) * sin(5.0 * x);
  
  // Funkcja 2: f(x) = 0.1 * x^2 - 0.5
  func f2 = 0.1 * pow(x, 2.0) - 0.5;

  // Równanie transcedentalne: chcemy znaleźć przecięcie f1(x) == f2(x)
  func diff = f1 - f2;

  // Znajdźmy miejsce zerowe funkcją z API
  // Newton-Raphson wymaga pochodnej, Bisection nie. Użyjmy Secant lub Bisection.
  double root_val = optimalization::bisection<double>(diff, 0.0, 3.0);
  
  std::cout << "Znalezione przecięcie w punkcie x = " << root_val << "\n";
  std::cout << "Wartość f1(x) = " << f1(root_val) << "\n";
  std::cout << "Wartość f2(x) = " << f2(root_val) << "\n";

  plot dashboard("ROOT Lab 9: Szukanie przecięcia krzywych");
  dashboard.set_xlabel("x").set_ylabel("f(x)");
  
  dashboard.add(f1, 0.0, 4.0, "f1(x) = e^{-x} sin(5x)", 300);
  dashboard.add(f2, 0.0, 4.0, "f2(x) = 0.1x^2 - 0.5", 300);
  
  vector<double> rx = {root_val};
  vector<double> ry = {f1(root_val)};
  dashboard.add_scatter(rx, ry, "Punkt przecięcia (Root)");

  dashboard.show();

  return 0;
}
