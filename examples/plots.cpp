#include <numc>

using namespace numc;

int main() {
  // Używamy globalnej zmiennej x zdefiniowanej w numc::x
  auto x = numc::x;

  func<double> f1 = sin(x) * (x / 2.0);
  func<double> f2 = exp(-0.1 * x) * cos(x) * 5.0;

  vector<double> data_x = {-8.0, -5.5, -2.1, 0.0, 3.14, 6.28, 9.0};
  vector<double> data_y = {3.5, -2.1, -1.0, 5.2, -4.8, 5.1, -4.0};

  plot dashboard("Analiza drgań");
  
  dashboard.set_xlabel("Czas [ms]").set_ylabel("Amplituda [V]");
  dashboard.add(f1, -10.0, 10.0, "Model teoretyczny");
  dashboard.add(f2, -10.0, 10.0, "Model z tłumieniem");
  dashboard.add_scatter(data_x, data_y, "Pomiary z czujnika");
  dashboard.show();

  // FAST PLOT: Szybkie rysowanie wykresu bez tworzenia obiektu plot.
  plt::show<double>({f1, f2}, -10.0, 10.0, "Szybki wykres", "Czas [ms]", "Amplituda [V]");

  return 0;
}