#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::analysis;

int main() {
  std::cout << "=== Projekt Analysis: Elektromagnetyzm ===\n" << std::endl;

  auto x = numc::x;
  auto y = numc::y;
  auto z = numc::z;

  // 1. Definicja potencjału elektrycznego V(x, y, z)
  // Załóżmy dipol w centrum lub inny ładunek, dla uproszczenia weźmiemy V(x,y,z) = 1 / sqrt(x^2 + y^2 + z^2)
  // A tutaj zdefiniujemy prostszy potencjał analityczny: V = x^2 + 2*y^2 + 3*z^2
  func<double> V = (x ^ 2.0) + 2.0 * (y ^ 2.0) + 3.0 * (z ^ 2.0);
  
  std::cout << "1. Pole Potencjalu Elektrycznego V(x, y, z):" << std::endl;
  std::cout << "   V(1, 2, 3) = " << V(1.0, 2.0, 3.0) << " V" << std::endl;

  // 2. Natezenie pola elektrycznego E = -nabla V
  auto E_field = nabla(V); // To zwraca pole wektorowe
  
  vector<double> point = {1.0, 2.0, 3.0};
  vector<double> E_val = E_field(point);
  
  // Ponieważ E = -gradient(V), musimy odwrócić znak wektora
  E_val = E_val * (-1.0);

  std::cout << "\n2. Natezenie pola elektrycznego E = -nabla V w punkcie (1, 2, 3):" << std::endl;
  std::cout << "   E = [" << E_val[0] << ", " << E_val[1] << ", " << E_val[2] << "] V/m" << std::endl;

  // 3. Dywergencja i Prawo Gaussa: nabla * E = rho / epsilon_0
  // Ponieważ nabla operuje na funkcjach, weźmy z powrotem negatywny gradient jako funkcję
  // Aby obliczyć dywergencję na polu, możemy przepisać pole do wektora funkcji
  std::vector<func<double>> E_func_vec = {
    -2.0 * x,       // -dV/dx
    -4.0 * y,       // -dV/dy
    -6.0 * z        // -dV/dz
  };
  
  func<double> div_E = nabla * E_func_vec;
  std::cout << "\n3. Dywergencja Pola E (Prawo Gaussa):" << std::endl;
  std::cout << "   nabla * E w (1, 2, 3) = " << div_E(1.0, 2.0, 3.0) << " V/m^2" << std::endl;

  // 4. Rotacja Pola E: nabla ^ E (zgodnie z Faradayem dla pola elektrostatycznego powinno być 0)
  auto rot_E = nabla ^ E_func_vec;
  vector<double> rot_E_val = rot_E(point);
  std::cout << "\n4. Rotacja Pola E (nabla ^ E):" << std::endl;
  std::cout << "   nabla ^ E w (1, 2, 3) = [" << rot_E_val[0] << ", " << rot_E_val[1] << ", " << rot_E_val[2] << "]" << std::endl;

  // 5. Calkowanie numeryczne - Potencjał na ścieżce (1D)
  // Załóżmy, że poruszamy się wzdłuż osi X od 0 do 5
  func<double> Vx = (x ^ 2.0); // Potencjał tylko wzdluz x
  double praca = adaptive_simpson<double>(Vx, 0.0, 5.0);
  std::cout << "\n5. Calkowanie: Praca pola na odcinku [0, 5] osi X:" << std::endl;
  std::cout << "   W = " << praca << " J" << std::endl;

  return 0;
}
