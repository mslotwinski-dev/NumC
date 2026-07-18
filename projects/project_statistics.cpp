#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::statistics;

int main() {
  std::cout << "=== Projekt Statistics: Monte Carlo Rozpadu Promieniotworczego ===\n" << std::endl;

  auto x = funcv::x;

  // 1. Modelowanie czasu rozpadu jądra jako rozkładu wykładniczego
  // Prawdopodobieństwo rozpadu w czasie t: f(t) = lambda * exp(-lambda * t)
  double lambda = 0.5; // Stała rozpadu
  
  std::cout << "1. Modelowanie rozpadu promieniotworczego dla lambda = " << lambda << std::endl;
  
  // Definiujemy dystrybuantę (PDF) i korzystamy z pięknego inicjatora `dist`
  dist<double> pdf_rozpad = lambda * exp(-lambda * x);
  
  // Próbkowanie z odrzuceniem dla t w zakresie [0, 10] s
  double max_pdf = lambda; // Maximum of lambda * e^(-lambda*t) is at t=0
  int liczba_czastek = 50000;
  
  std::cout << "   Symulacja dla " << liczba_czastek << " cząstek w przedziale czasowym [0, 10] s..." << std::endl;
  std::vector<double> czasy_rozpadu = rejection_sample_1d(pdf_rozpad, liczba_czastek, 0.0, 10.0, max_pdf);

  // 2. Budowanie i analiza Histogramu
  std::cout << "\n2. Budowa histogramu czasów rozpadu (10 koszykow)" << std::endl;
  Histogram1D<double> hist(10, 0.0, 10.0, czasy_rozpadu);
  
  hist.print();

  // 3. Regresja liniowa dla zlogarytmowanych danych N(t)
  // Teoretycznie: N(t) = N0 * exp(-lambda * t) -> ln(N(t)) = ln(N0) - lambda * t
  std::cout << "\n3. Aproksymacja stalej rozpadu (Regresja Liniowa na podstawie histogramu):" << std::endl;
  
  vector<double> t_punkty;
  vector<double> ln_N_punkty;
  
  auto bin_centers = hist.get_bin_centers();
  auto bin_contents = hist.get_bin_contents();

  for (size_t i = 0; i < hist.get_num_bins(); ++i) {
    if (bin_contents[i] > 0) {
      t_punkty.push_back(bin_centers[i]);
      ln_N_punkty.push_back(std::log(bin_contents[i]));
    }
  }

  auto [m, b, r2] = linear_regression<double>(t_punkty, ln_N_punkty);
  
  std::cout << "   Dopasowana prosta: ln(N) = " << m << " * t + " << b << std::endl;
  std::cout << "   Estymowana stala rozpadu lambda (powinno byc ~" << lambda << "): " << -m << std::endl;

  return 0;
}
