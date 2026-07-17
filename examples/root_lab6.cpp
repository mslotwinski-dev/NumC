#include "../src/numc.h"
#include <iostream>
#include <cmath>
#include <random>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 6 Replication: Całkowanie Monte Carlo ==========\n";

  // Cel: Obliczyć objętość bryły w 3D określonej przez przenikające się walce.
  // x^2 + y^2 <= 1  oraz  x^2 + z^2 <= 1
  // Analityczna objętość (Stejnmetza bryła) to V = 16/3 ~= 5.33333...

  // Objętość pudła otaczającego: [-1, 1]x[-1, 1]x[-1, 1]
  
  analysis::mfunc<double> inside_volume = std::function<double(const numc::vector<double>&)>([](const numc::vector<double>& v) {
    if ((v[0] * v[0] + v[1] * v[1] <= 1.0) && (v[0] * v[0] + v[2] * v[2] <= 1.0)) return 1.0;
    return 0.0;
  });

  std::vector<int> N_values = {1000, 10000, 100000, 1000000};
  int num_vals = N_values.size();
  vector<double> errors(num_vals);
  vector<double> N_plot(num_vals);

  int idx = 0;
  vector<double> lower_bound = {-1.0, -1.0, -1.0};
  vector<double> upper_bound = {1.0, 1.0, 1.0};

  for (int N : N_values) {
    auto [v_mc, error] = inside_volume.integral_mc(lower_bound, upper_bound, N);
    
    std::cout << "Dla N = " << N << " :\n";
    std::cout << "  Zmierzone V: " << v_mc << " +/- " << error << "\n";
    std::cout << "  Błąd wzgledny prawdziwy: " << std::abs(v_mc - (16.0/3.0)) / (16.0/3.0) * 100.0 << "%\n";
    
    N_plot[idx] = std::log10(static_cast<double>(N));
    errors[idx] = std::log10(error);
    idx++;
  }

  // WIZUALIZACJA: błąd spadający wg 1/sqrt(N)
  plot dashboard("ROOT Lab 6: Zależność błędu MC od N");
  dashboard.set_xlabel("log10(N)").set_ylabel("log10(Error)");
  
  auto x = func<>::x();
  // Theoretical line: log10(C / sqrt(N)) = log10(C) - 0.5 * log10(N)
  func theoretical_error = 0.576 - 0.5 * x;

  dashboard.add_scatter(N_plot, errors, "Błędy pomiarowe MC");
  dashboard.add(theoretical_error, 2.0, 7.0, "Teoretyczne tempo zbieżności -0.5", 100);

  dashboard.show();

  return 0;
}
