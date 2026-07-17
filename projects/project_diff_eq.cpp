#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::differential_equations;

int main() {
  std::cout << "=== Projekt Differential Equations: Atraktor Lorenza ===\n" << std::endl;

  // Parametry układu Lorenza
  double sigma = 10.0;
  double rho = 28.0;
  double beta = 8.0 / 3.0;

  // Układ 3 nieliniowych równań różniczkowych zwyczajnych (ODE)
  // dy/dt = F(t, y)
  // y = [X, Y, Z]
  auto lorenz = [sigma, rho, beta](double t, const vector<double>& y) -> vector<double> {
    double X = y[0];
    double Y = y[1];
    double Z = y[2];

    double dX_dt = sigma * (Y - X);
    double dY_dt = X * (rho - Z) - Y;
    double dZ_dt = X * Y - beta * Z;

    return {dX_dt, dY_dt, dZ_dt};
  };

  // Warunki początkowe
  double t0 = 0.0;
  double t_end = 50.0;
  double dt = 0.01;
  vector<double> y0 = {1.0, 1.0, 1.0};

  std::cout << "1. Symulacja chaotycznego Atraktora Lorenza" << std::endl;
  std::cout << "   Rozwiazywanie za pomoca metody Runge-Kutty 4. rzędu (RK4)..." << std::endl;

  auto [t_vals, y_vals] = run_kut4<double>(lorenz, t0, y0, t_end, dt);

  std::cout << "   Symulacja zakonczona. Wygenerowano " << t_vals.size() << " krokow czasowych." << std::endl;
  std::cout << "\n2. Ostatnie 5 stanow wektora stanu [X, Y, Z]:" << std::endl;
  
  size_t n = t_vals.size();
  for (size_t i = n > 5 ? n - 5 : 0; i < n; ++i) {
    std::cout << "   t = " << std::fixed << std::setprecision(2) << t_vals[i] 
              << " | X = " << std::setw(8) << y_vals[i][0] 
              << ", Y = " << std::setw(8) << y_vals[i][1] 
              << ", Z = " << std::setw(8) << y_vals[i][2] << std::endl;
  }

  std::cout << "\n3. Zapis do wykresu (wykorzystanie modulu wizualizacji)..." << std::endl;
  // Przygotowanie wektorów do wykresu szybkiego (Z vs X - motyl Lorenza)
  std::vector<double> X_plt, Z_plt;
  for (const auto& y : y_vals) {
    X_plt.push_back(y[0]);
    Z_plt.push_back(y[2]);
  }

  // Odkomentuj ponizsze, aby wygenerowac szybki scatter plot w terminalu
  // plt::show(vector<double>(X_plt), vector<double>(Z_plt), "Atraktor Lorenza (XZ)", "X", "Z");
  
  std::cout << "   [OK] Wyniki sa gotowe do wyswietlenia." << std::endl;

  return 0;
}
