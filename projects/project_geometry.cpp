#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::geometry;

int main() {
  std::cout << "=== Projekt Geometry: Omijanie Przeszkód i Kolizje (Robotyka) ===\n" << std::endl;

  // Reprezentacja przeszkód jako zbior punktów w przestrzeni roboczej 2D
  std::vector<Point<double>> chmura_punktow = {
    {1.0, 1.0}, {2.0, 0.5}, {3.0, 1.2}, 
    {2.5, 2.5}, {1.5, 2.0}, {1.8, 1.5},
    {2.2, 1.0}, {2.8, 0.8}
  };

  std::cout << "1. Algorytm Convex Hull - wyznaczanie zewnetrznej granicy strefy zakazanej" << std::endl;
  auto otoczka = convex_hull(chmura_punktow);
  
  std::cout << "   Otoczka wypukla uksztaltowana z " << otoczka.size() << " wierzcholkow:" << std::endl;
  for (const auto& pt : otoczka) {
    std::cout << "   (" << pt[0] << ", " << pt[1] << ")" << std::endl;
  }

  // Symulacja poruszania sie ramienia robota lub drona
  Point<double> pozycja_robota = {2.0, 1.5}; // Punkt wewnatrz
  Point<double> pozycja_celu = {4.0, 3.0};   // Punkt na zewnatrz
  
  std::cout << "\n2. Algorytm Point In Polygon - sprawdzanie kolizji ze strefa" << std::endl;
  bool czy_kolizja = point_in_polygon(pozycja_robota, otoczka);
  std::cout << "   Czy robot w (2.0, 1.5) jest w strefie zakazanej? " << (czy_kolizja ? "TAK (KOLIZJA)" : "NIE") << std::endl;
  
  bool cel_kolizja = point_in_polygon(pozycja_celu, otoczka);
  std::cout << "   Czy cel w (4.0, 3.0) jest w strefie zakazanej? " << (cel_kolizja ? "TAK (KOLIZJA)" : "NIE") << std::endl;

  // Odleglosc promienia lasera do przeszkody
  Point<double> laser_start = {0.0, 0.0};
  Point<double> laser_end = {3.0, 3.0};
  Point<double> wierzcholek_przeszkody = {1.5, 2.0}; // Jeden z punktow otoczki
  
  std::cout << "\n3. Detekcja bliskosci - Ray Casting do przeszkody" << std::endl;
  double odleglosc = distance_point_to_segment(wierzcholek_przeszkody, laser_start, laser_end);
  std::cout << "   Odleglosc przeszkody w (1.5, 2.0) od promienia lasera: " << odleglosc << " m" << std::endl;

  return 0;
}
