#include "../src/numc.h"
#include <iostream>
#include <cmath>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 1 Replication (Beautiful API) ==========\n";
  const double PI = std::acos(-1.0);

  // a) Normalize PDF: f(x) = a * x * sin(x)
  auto x = func<>::x();
  func raw_pdf = x * sin(x);
  
  double integral_val = analysis::adaptive_simpson<double>(raw_pdf, 0.0, PI);
  double a = 1.0 / integral_val;
  std::cout << "Normalizacja stała 'a': " << a << "\n";
  
  func f = a * x * sin(x);

  // b) Dystrybuanta analityczna: F(x) = \int_0^x f(t) dt
  auto cdf = analysis::integral<double>(f, 0.0);
  std::cout << "Dystrybuanta w x=PI/2: " << cdf(PI/2.0) << "\n";

  // c) Wartość oczekiwana analitycznie
  func expected_func = x * f;
  double expected_val_analytical = analysis::adaptive_simpson<double>(expected_func, 0.0, PI);
  std::cout << "Wartość oczekiwana analitycznie: " << expected_val_analytical << "\n";

  // Numerycznie (za pomocą losowania Rejection Sampling)
  double pdf_max = a * PI;
  vector<double> samples = statistics::rejection_sample_1d<double>(f, 100000, 0.0, PI, pdf_max);
  
  statistics::Histogram1D<double> hist(100, 0.0, PI);
  for (double val : samples) hist.fill(val);
  hist.normalize();

  double expected_val_numerical = hist.get_mean();
  std::cout << "Wartość oczekiwana numerycznie: " << expected_val_numerical << "\n";

  // d) Numeryczne wlasciwosci
  double p_1_2 = analysis::adaptive_simpson<double>(f, 1.0, 2.0);
  std::cout << "P(1 <= x <= 2): " << p_1_2 << "\n";
  std::cout << "Wariancja z histogramu: " << hist.get_variance() << "\n";

  std::cout << "Mediana: " << statistics::percentile(samples, 50.0) << "\n";
  std::cout << "Kwartyl dolny: " << statistics::percentile(samples, 25.0) << "\n";
  std::cout << "Decyl 0.7: " << statistics::percentile(samples, 70.0) << "\n";

  // PIĘKNE RYSOWANIE:
  plot dashboard("ROOT Lab 1: PDF i CDF");
  dashboard.set_xlabel("x").set_ylabel("Wartość");
  dashboard.add(f, 0.0, PI, "f(x) - Gęstość prawdopodobieństwa", 200);
  dashboard.add(cdf, 0.0, PI, "F(x) - Dystrybuanta", 200);

  // Dodajemy wylosowany histogram jako scatter
  dashboard.add_scatter(hist.get_bin_centers(), hist.get_bin_contents(), "Histogram z Rejection Sampling");
  
  dashboard.show(); // Zapisze wykres jako plik HTML

  return 0;
}
