#include "../src/numc.h"
#include <iostream>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 5 Replication: Curve Fitting ==========\n";

  double true_a = 2.0, true_b = -1.5, true_c = 4.0;
  
  auto x = func<>::x();
  func true_model = true_a * pow(x, 2.0) + true_b * x + true_c;

  // Generowanie danych z modelu analitycznego
  vector<double> x_data = linspace(-10.0, 10.0, 50);
  vector<double> y_data = true_model(x_data) + random::normal(0.0, 5.0, 50);

  // Definicja parametrów (z wartościami początkowymi)
  statistics::param<double> a("a", 1.0);
  statistics::param<double> b("b", 1.0);
  statistics::param<double> c("c", 1.0);

  // Model dopasowania jako funkcja lambda (lub func)
  auto model = [&](double x_val) {
    return a.value() * x_val * x_val + b.value() * x_val + c.value();
  };

  // Dopasowanie modelu do danych
  statistics::fit<double>(model, x_data, y_data, {&a, &b, &c});
  
  std::cout << "Prawdziwe parametry: a=" << true_a << " b=" << true_b << " c=" << true_c << "\n";
  std::cout << "Dopasowane parametry: a=" << a.value() << " b=" << b.value() << " c=" << c.value() << "\n";

  func fitted_model = a.value() * pow(x, 2.0) + b.value() * x + c.value();

  // WIZUALIZACJA:
  plot dashboard("ROOT Lab 5: Dopasowanie krzywych (Least Squares)");
  dashboard.set_xlabel("x").set_ylabel("y");
  
  dashboard.add_scatter(x_data, y_data, "Zaszumione dane pomiarowe");
  dashboard.add(true_model, -10.0, 10.0, "Prawdziwy Model", 200);
  dashboard.add(fitted_model, -10.0, 10.0, "Dopasowany Model (LM)", 200);
  
  dashboard.show();

  return 0;
}
