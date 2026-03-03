#include <numc>

using namespace numc;

int main() {
  auto x = func<>::x();

  func f1 = sin(x) * (x / 2.0);
  func f2 = exp(-0.1 * x) * cos(x) * 5.0;

  vector<double> data_x = {-8.0, -5.5, -2.1, 0.0, 3.14, 6.28, 9.0};
  vector<double> data_y = {3.5, -2.1, -1.0, 5.2, -4.8, 5.1, -4.0};

  plot dashboard("Analiza drgań");
  dashboard.set_xlabel("Czas [ms]").set_ylabel("Amplituda [V]");
  dashboard.add(f1, -10.0, 10.0, "Model teoretyczny");
  dashboard.add(f2, -10.0, 10.0, "Model z tłumieniem");
  dashboard.add_scatter(data_x, data_y, "Pomiary z czujnika");
  dashboard.show("raport_z_eksperymentu.html");

  return 0;
}