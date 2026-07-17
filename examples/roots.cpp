#include <iostream>
#include <iomanip>

#include <numc>

using namespace numc;
using namespace numc::optimalization;

int main() {
  auto x = func<>::x();

  // 1. Funkcja testowa do poszukiwania miejsc zerowych: f(x) = x^3 - 2x - 5
  // Spodziewany pierwiastek: x ≈ 2.09455

  func f = (x ^ 3) - 2.0 * x - 5.0;

  // 2. Funkcja do metody punktu stałego: g(x) = cos(x)
  // Szuka punktu, w którym x = cos(x)
  // Spodziewany pierwiastek: x ≈ 0.739085
  
  func g = cos(x);

  std::cout << std::fixed << std::setprecision(6);
  
  // --- Metody przedziałowe (wymagają podania [a, b], gdzie funkcja zmienia znak) ---
  
  // Metoda równego podziału (Bisection)
  double root_bisection = bisection(f, 2.0, 3.0);
  std::cout << "[Bisection]           x = " << root_bisection 
            << " | f(x) = " << f(root_bisection) << "\n";

  // Metoda fałszywego założenia (Regula Falsi)
  double root_regula = regula_falsi(f, 2.0, 3.0);
  std::cout << "[Regula Falsi]        x = " << root_regula 
            << " | f(x) = " << f(root_regula) << "\n";


  // --- Metody otwarte (wymagają punktów startowych, zazwyczaj szybsza zbieżność) ---

  // Metoda punktu stałego (Fixed-point iteration) - rozwiązuje x = g(x)
  double root_fixed = fixed_point(g, 1.0);
  std::cout << "[Fixed-Point (cos)]   x = " << root_fixed 
            << " | cos(x) - x = " << (g(root_fixed) - root_fixed) << "\n";

  // Metoda Newtona-Raphsona (wymaga 1 punktu startowego, wewn. numeryczna pochodna)
  double root_newton = newton(f, 2.0);
  std::cout << "[Newton-Raphson]      x = " << root_newton 
            << " | f(x) = " << f(root_newton) << "\n";

  // Metoda Halleya (zbieżność sześcienna, wykorzystuje 1 i 2 pochodną numerycznie)
  double root_halley = halley(f, 2.0);
  std::cout << "[Halley's Method]     x = " << root_halley 
            << " | f(x) = " << f(root_halley) << "\n";

  // Metoda siecznych (wymaga 2 punktów startowych)
  double root_secant = secant(f, 2.0, 3.0);
  std::cout << "[Secant Method]       x = " << root_secant 
            << " | f(x) = " << f(root_secant) << "\n";

  // Odwrotna interpolacja kwadratowa (Inverse Quadratic, wymaga 3 punktów startowych)
  double root_iqi = inverse_quadratic(f, 1.0, 2.0, 3.0);
  std::cout << "[Inverse Quadratic]   x = " << root_iqi 
            << " | f(x) = " << f(root_iqi) << "\n";

  return 0;
}