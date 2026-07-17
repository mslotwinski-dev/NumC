#include "../src/numc.h"
#include <iostream>
#include <random>

using namespace numc;

// Prosty (zły) LCG dla celów dydaktycznych: X_{n+1} = (A * X_n + C) mod M
class BadLCG {
  unsigned long long state;
public:
  BadLCG(unsigned long long seed) : state(seed) {}
  double operator()() {
    state = (65 * state + 17) % 2048; // Słabe parametry, widać wyraźne wzorce
    return static_cast<double>(state) / 2048.0;
  }
};

int main() {
  std::cout << "========== ROOT Lab 8 Replication: Generatory Liczb Losowych ==========\n";

  BadLCG bad_rng(42);
  int NUM_POINTS = 5000;
  vector<double> bad_x(NUM_POINTS), bad_y(NUM_POINTS);
  vector<double> good_x = random::uniform(0.0, 1.0, NUM_POINTS);
  vector<double> good_y = random::uniform(0.0, 1.0, NUM_POINTS);

  for (int i = 0; i < NUM_POINTS; ++i) {
    bad_x[i] = bad_rng();
    bad_y[i] = bad_rng();
  }

  // WIZUALIZACJA PORÓWNAWCZA
  plot dashboard("ROOT Lab 8: LCG vs Mersenne Twister");
  dashboard.set_xlabel("X_i").set_ylabel("X_{i+1}");
  
  dashboard.add_scatter(bad_x, bad_y, "Bad LCG - widoczne wzorce i linie");
  dashboard.add_scatter(good_x, good_y, "Mersenne Twister - jednorodny szum");

  dashboard.show();

  return 0;
}
