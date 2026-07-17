#include "../src/numc.h"
#include <iostream>
#include <cmath>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 3 Replication (Beautiful API) ==========\n";
  const double PI_2 = std::acos(0.0); // pi/2

  // f(x, y) = c * sin(x+y) / sqrt(x^2 + y)
  analysis::mfunc<double> raw_pdf = std::function<double(const numc::vector<double>&)>([](const numc::vector<double>& v) {
    if (v[0] == 0.0 && v[1] == 0.0) return 0.0;
    return std::sin(v[0] + v[1]) / std::sqrt(v[0] * v[0] + v[1]);
  });

  // Normalize
  double integral_val = raw_pdf.integral_2d(0.0, PI_2, 0.0, PI_2);
  double c = 1.0 / integral_val;
  std::cout << "Normalizacja stała 'c': " << c << "\n";

  analysis::mfunc<double> f = std::function<double(const numc::vector<double>&)>([c](const numc::vector<double>& v) {
    if (v[0] == 0.0 && v[1] == 0.0) return 0.0;
    return c * std::sin(v[0] + v[1]) / std::sqrt(v[0] * v[0] + v[1]);
  });

  double pdf_max = c * 2.0;

  std::cout << "Losowanie 200,000 próbek...\n";
  auto samples = statistics::rejection_sample_2d<double>(f, 200000, 0.0, PI_2, 0.0, PI_2, pdf_max);

  // Auto-fill Histogram 2D
  statistics::Histogram2D<double> hist2d(100, 0.0, PI_2, 100, 0.0, PI_2, samples);
  hist2d.normalize();

  // Marginal projections
  auto g_x = hist2d.projection_x();
  auto h_y = hist2d.projection_y();

  std::cout << "\nParametry wyliczone na podstawie Histogramu 2D:\n";
  std::cout << "E(X): " << hist2d.mean_x() << "\n";
  std::cout << "E(Y): " << hist2d.mean_y() << "\n";
  std::cout << "sigma(X): " << std::sqrt(hist2d.variance_x()) << "\n";
  std::cout << "sigma(Y): " << std::sqrt(hist2d.variance_y()) << "\n";
  std::cout << "cov(X,Y): " << hist2d.covariance() << "\n";
  std::cout << "rho(X,Y): " << hist2d.correlation() << "\n";

  // RYSOWANIE: Gęstości brzegowe z 2D
  plot dashboard("ROOT Lab 3: Gęstości brzegowe X oraz Y");
  dashboard.set_xlabel("Wartość").set_ylabel("Gęstość");
  
  dashboard.add(g_x, "Histogram g(x)");
  dashboard.add(h_y, "Histogram h(y)");
  
  dashboard.show();

  return 0;
}
