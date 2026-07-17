#include <iostream>
#include <iomanip>
#include <numc>

using namespace numc;

int main() {
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "╔══════════════════════════════════════════════════╗\n";
  std::cout << "║           NumC — Showcase pięknego API           ║\n";
  std::cout << "╚══════════════════════════════════════════════════╝\n\n";

  // ── 1. Wyrażenia funkcyjne ──────────────────────────────────────────────
  auto x = numc::x();

  func f1 = sin(x) * (x / 2.0);
  func f2 = exp(-0.1 * x) * cos(x) * 5.0;
  func f3 = derivative(f2);

  std::cout << "f1(PI/2)  = " << f1(PI / 2) << "\n";
  std::cout << "f2(0)     = " << f2(0.0) << "\n";
  std::cout << "f2'(0)    = " << f3(0.0) << "  (numeryczna pochodna)\n";
  std::cout << "∫f1 dx    = " << f1.integral(-PI, PI) << "  na [-π, π]\n\n";

  // ── 2. Wektory ──────────────────────────────────────────────────────────
  vector v1 = {1.0, 2.0, 3.0};
  vector v2 = {4.0, 5.0, 6.0};

  vector v3 = v1 + v2;
  vector v4 = v1 * 2.0;
  double dot = v1 * v2;
  vector cross = v1 ^ v2;

  std::cout << "v1 + v2   = " << v3 << "\n";
  std::cout << "v1 · v2   = " << dot << "\n";
  std::cout << "v1 × v2   = " << cross << "\n\n";

  // ── 3. Pomiary + statystyka ─────────────────────────────────────────────
  vector data_x = {-8.0, -5.5, -2.1, 0.0, 3.14, 6.28, 9.0};
  vector data_y = {3.5, -2.1, -1.0, 5.2, -4.8, 5.1, -4.0};

  std::cout << "Średnia pomiarów y = " << mean(data_y) << "\n";
  std::cout << "Odch. std. y       = " << std_dev(data_y) << "\n";

  auto [a, b, r2] = linear_regression(data_x, data_y);
  std::cout << "Regresja liniowa: y = " << a << " + " << b << "·x  (R² = " << r2 << ")\n\n";

  // ── 4. Miejsca zerowe ───────────────────────────────────────────────────
  func poly = (x ^ 3) - 2.0 * x - 5.0;
  double root = newton(poly, 2.0);
  std::cout << "Newton: pierwiastek x³−2x−5 ≈ " << root << "  (f(x) = " << poly(root) << ")\n\n";

  // ── 5. Wizualizacja ─────────────────────────────────────────────────────
  std::cout << "Otwieram wykres interaktywny w przeglądarce…\n";

  plot dashboard("Analiza drgań");

  dashboard.set_xlabel("Czas [ms]")
           .set_ylabel("Amplituda [V]")
           .add(f1, -10.0, 10.0, "Model teoretyczny")
           .add(f2, -10.0, 10.0, "Model z tłumieniem")
           .add(data_x, data_y, "Pomiary z czujnika");
  dashboard.show();

  // FAST PLOT: jedna linia, bez tworzenia obiektu plot.
  plt::show({f1, f2}, -10.0, 10.0, "Szybki wykres", "Czas [ms]", "Amplituda [V]",
            data_x, data_y, "Pomiary z czujnika");

  std::cout << "\nGotowe. Miłej pracy z NumC!\n";
  return 0;
}
