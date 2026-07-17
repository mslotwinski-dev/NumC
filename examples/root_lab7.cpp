#include "../src/numc.h"
#include <iostream>
#include <cmath>
#include <random>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 7 Replication: Centralne Twierdzenie Graniczne ==========\n";

  const int NUM_SAMPLES = 100000;
  
  // Testujemy sumy dla K=1, K=2, K=10 zmiennych losowych z użyciem arytmetyki wektorowej
  vector<double> sum_K1 = random::uniform(-0.5, 0.5, NUM_SAMPLES);
  vector<double> sum_K2 = sum_K1 + random::uniform(-0.5, 0.5, NUM_SAMPLES);
  
  vector<double> sum_K10 = sum_K2;
  for (int i = 0; i < 8; ++i) {
    sum_K10 = sum_K10 + random::uniform(-0.5, 0.5, NUM_SAMPLES);
  }

  statistics::Histogram1D<double> hist1(50, -1.0, 1.0);
  statistics::Histogram1D<double> hist2(50, -2.0, 2.0);
  statistics::Histogram1D<double> hist10(50, -5.0, 5.0);

  for (int i = 0; i < NUM_SAMPLES; ++i) {
    hist1.fill(sum_K1[i]);
    hist2.fill(sum_K2[i]);
    hist10.fill(sum_K10[i]);
  }
  hist1.normalize();
  hist2.normalize();
  hist10.normalize();

  // WIZUALIZACJA
  plot dashboard("ROOT Lab 7: Centralne Twierdzenie Graniczne");
  dashboard.set_xlabel("Suma K zmiennych jednostajnych").set_ylabel("Znormalizowana gęstość");

  dashboard.add_scatter(hist1.get_bin_centers(), hist1.get_bin_contents(), "K = 1");
  dashboard.add_scatter(hist2.get_bin_centers(), hist2.get_bin_contents(), "K = 2");
  dashboard.add_scatter(hist10.get_bin_centers(), hist10.get_bin_contents(), "K = 10");

  // Teoretyczny rozkład Gaussa dla K=10.
  // Uniform(-0.5, 0.5) ma wariancję 1/12, średnią 0.
  // Suma 10 ma średnią 0, wariancję 10/12 = 5/6
  double sigma = std::sqrt(5.0 / 6.0);
  auto x = func<>::x();
  func gaussian = (1.0 / (sigma * std::sqrt(2.0 * std::acos(-1.0)))) * exp(-0.5 * pow(x / sigma, 2.0));
  
  dashboard.add(gaussian, -5.0, 5.0, "Rozkład Normalny Teoretyczny (K=10)");

  dashboard.show();

  return 0;
}
