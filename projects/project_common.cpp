#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;

int main() {
  std::cout << "=== Projekt Common: Termodynamika i Wielomiany ===\n" << std::endl;

  // 1. DataFrame - Wczytanie i symulacja danych pomiarowych (Temperatura w K vs Objętość w L)
  // V = V0 * (1 + alpha * T + beta * T^2)
  DataFrame df;
  df.add_column("Temperatura", {273.15, 298.15, 323.15, 348.15, 373.15});
  df.add_column("Objetosc", {22.4, 24.5, 26.8, 29.3, 31.9}); // Symulowane dane
  
  std::cout << "1. Pomiary termodynamiczne (Gaz doskonały i nieliniowości):" << std::endl;
  // df.print() is not available, iterate and print manually if needed
  
  // 2. Polynomial - Dopasowanie wielomianu (zwykła inicjalizacja przykładowych współczynników dla uproszczenia)
  // Zakładamy, że wyliczyliśmy V(T) = c0 + c1*T + c2*T^2
  polynomial<double> V_T({-0.1, 0.082, 0.0001}); 
  std::cout << "\n2. Równanie stanu jako wielomian V(T) = c0 + c1*T + c2*T^2:" << std::endl;
  std::cout << "   V(298.15 K) = " << V_T(298.15) << " L" << std::endl;

  // 3. Pochodna numeryczna do wyznaczenia nieliniowego współczynnika rozszerzalności cieplnej: alpha(T) = (1/V) * dV/dT
  std::cout << "\n3. Różniczkowanie numeryczne do wyznaczenia dV/dT:" << std::endl;
  double T_eval = 298.15;
  func<double> V_func = V_T.to_func();
  
  double dV_dT = numc::analysis::derivative_central_O4(V_func, T_eval);
  double alpha = (1.0 / V_T(T_eval)) * dV_dT;
  
  std::cout << "   T = " << T_eval << " K" << std::endl;
  std::cout << "   dV/dT = " << dV_dT << " L/K" << std::endl;
  std::cout << "   Wspolczynnik rozszerzalnosci cieplnej alpha(T) = " << alpha << " 1/K" << std::endl;

  // 4. Sparse Matrix - System równań stanu dla wielu gazów
  std::cout << "\n4. Macierz Rzadka - Układ równań przepływu ciepła w sieci:" << std::endl;
  sparse_matrix<double> heat_transfer(4, 4, {
    {0, 0, 2.0}, {0, 1, -1.0},
    {1, 0, -1.0}, {1, 1, 2.0}, {1, 2, -1.0},
    {2, 1, -1.0}, {2, 2, 2.0}, {2, 3, -1.0},
    {3, 2, -1.0}, {3, 3, 2.0}
  });
  
  std::cout << heat_transfer << std::endl;

  return 0;
}
