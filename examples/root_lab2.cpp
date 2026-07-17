#include "../src/numc.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 2 Replication (Beautiful API) ==========\n";
  
  statistics::param<double> m("m", 1.0);
  statistics::param<double> n("n", 2.0);

  // f(x, y) = c * (e^(mx) + sin(ny))
  analysis::mfunc<double> raw_pdf = std::function<double(const numc::vector<double>&)>([m, n](const numc::vector<double>& v) {
    return std::exp(m.value() * v[0]) + std::sin(n.value() * v[1]);
  });

  // a) Normalize
  double integral_val = raw_pdf.integral_2d(0.0, 1.0, 0.0, 1.0);
  double c = 1.0 / integral_val;
  std::cout << "Normalizacja stała 'c': " << c << "\n";

  analysis::mfunc<double> f = std::function<double(const numc::vector<double>&)>([c, m, n](const numc::vector<double>& v) {
    return c * (std::exp(m.value() * v[0]) + std::sin(n.value() * v[1]));
  });

  double pdf_max = c * (std::exp(1.0) + 1.0);

  std::cout << "Pobieranie 100,000 próbek...\n";
  auto samples = statistics::rejection_sample_2d<double>(f, 100000, 0.0, 1.0, 0.0, 1.0, pdf_max);

  // New beautiful histogram constructor!
  statistics::Histogram2D<double> hist2d(50, 0.0, 1.0, 50, 0.0, 1.0, samples);
  hist2d.normalize();

  // b) Dystrybuanta z histogramu (numerycznie)
  auto cdf2d = hist2d.get_cdf();
  std::cout << "Dystrybuanta CDF 2D w punkcie (0.5, 0.5): " << cdf2d.get_bin_content(24, 24) << "\n";

  // c) Gęstość brzegowa g(x)
  auto hist_gx = hist2d.projection_x();
  
  auto x = func<>::x();
  func ga = (c / n.value()) * (n.value() * exp(m.value() * x) - std::cos(n.value()) + 1.0);

  std::cout << "\nPorównanie gęstości brzegowej g(x) przy x=0.5:\n";
  std::cout << " - Analitycznie g_a(0.5): " << ga(0.5) << "\n";
  std::cout << " - Numerycznie z Histogramu: " << hist_gx.get_bin_content(24) << "\n"; 

  // d) Gęstość brzegowa h(y)
  auto hist_hy = hist2d.projection_y();
  std::cout << "Numerycznie gęstość brzegowa h(y) przy y=0.5: " << hist_hy.get_bin_content(24) << "\n";

  // RYSOWANIE: Porównanie modelu z danymi numerycznymi
  plot dashboard("ROOT Lab 2: Gęstość brzegowa g(x)");
  dashboard.set_xlabel("x").set_ylabel("g(x)");
  
  // Natively plotting histogram
  dashboard.add(hist_gx, "Histogram g(x)");
  dashboard.add(ga, 0.0, 1.0, "Analityczna gęstość g_a(x)", 100);
  
  dashboard.show();

  return 0;
}
