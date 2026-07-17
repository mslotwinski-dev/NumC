#include "../src/numc.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace numc;
using namespace numc::optimalization;

int main() {
  std::cout << "=== Projekt Optimalization: Optymalizacja Trajektorii Rakiety ===\n" << std::endl;

  // Cel: Zminimalizować całkowite zużycie paliwa rakiety potrzebne do wejścia na orbitę,
  // zachowując jednocześnie ostateczną prędkość i pułap (wymuszone przez karę w funkcji kosztu).

  // Zmienne decyzyjne (wektor x):
  // x[0] - kąt pochylenia startowego (pitch angle) w radianach (np. 1.0 rad)
  // x[1] - początkowy ciąg silnika w kN (np. 500.0)
  
  // Modelujemy złożoną, nieliniową funkcję kosztu
  auto funkcja_kosztu = [](const vector<double>& v) -> double {
    double kat = v[0];
    double ciag = v[1];
    
    // Złapmy proste granice z karą (kara za ujemny kąt lub ciąg poniżej 100)
    double kara = 0.0;
    if (kat < 0.0) kara += 1e6 * (0.0 - kat);
    if (kat > 1.57) kara += 1e6 * (kat - 1.57); // 90 stopni
    if (ciag < 100.0) kara += 1e6 * (100.0 - ciag);
    if (ciag > 1000.0) kara += 1e6 * ciag;

    // Zużycie paliwa (model)
    double paliwo = ciag * 0.05 + 100.0 / std::cos(kat + 0.1);

    // Błąd orbity (symulowana różnica od wymaganej predkosci orbity zalezna od kata i ciagu)
    double predkosc_docelowa = 7.8; // km/s
    double predkosc_osiagnieta = ciag * 0.015 * std::sin(kat); 
    
    double blad_orbity = std::pow(predkosc_docelowa - predkosc_osiagnieta, 2);

    // Całkowity koszt: paliwo + duża waga dla błędu orbity + kary granic
    return paliwo + 500.0 * blad_orbity + kara;
  };

  vector<double> start_guess = {1.0, 500.0}; // Początkowe przybliżenie
  double h = 10.0; // Rozmiar sympleksu

  std::cout << "1. Uruchamianie algorytmu Downhill Simplex (Nelder-Mead)..." << std::endl;
  std::cout << "   Poczatkowy punkt: kat = " << start_guess[0] << " rad, ciag = " << start_guess[1] << " kN" << std::endl;
  std::cout << "   Poczatkowy koszt: " << funkcja_kosztu(start_guess) << std::endl;

  auto [optimum, iterations] = nelder_mead<double>(funkcja_kosztu, start_guess, h);

  std::cout << "\n2. Znaleziono optymalne parametry po " << iterations << " iteracjach:" << std::endl;
  std::cout << "   Optymalny kat pochylenia: " << optimum[0] << " rad (" << optimum[0] * 180.0 / 3.14159 << " deg)" << std::endl;
  std::cout << "   Optymalny ciag silnika: " << optimum[1] << " kN" << std::endl;
  std::cout << "   Minimalny koszt (paliwo + kary): " << funkcja_kosztu(optimum) << std::endl;

  return 0;
}
