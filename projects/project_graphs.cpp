#include "../src/numc.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace numc;
using namespace numc::graphs;

int main() {
  std::cout << "=== Projekt Graphs: Analiza Sieci Energetycznej ===\n" << std::endl;

  // Reprezentacja węzłów (miasta, elektrownie, stacje przesyłowe)
  std::vector<std::string> wezly = {
    "Elektrownia Glowna",   // 0
    "Stacja Transformatorowa A", // 1
    "Stacja Transformatorowa B", // 2
    "Miasto X",             // 3
    "Miasto Y",             // 4
    "Miasto Z"              // 5
  };

  // Graf skierowany, wagi krawędzi oznaczają stratę przesyłu (np. opór linii, spadek napięcia)
  Graph<double> power_grid(wezly.size());
  
  // Elektrownia wysyła prąd do stacji
  power_grid.add_edge(0, 1, 0.5, true); 
  power_grid.add_edge(0, 2, 0.8, true);
  
  // Stacje rozsyłają do miast i między sobą
  power_grid.add_edge(1, 2, 0.2, true);
  power_grid.add_edge(1, 3, 1.5, true);
  power_grid.add_edge(1, 4, 2.0, true);
  
  power_grid.add_edge(2, 4, 0.5, true);
  power_grid.add_edge(2, 5, 1.2, true);
  
  power_grid.add_edge(4, 3, 0.3, true); // Zapasowa linia z Y do X

  std::cout << "1. Model Sieci Przesylowej zbudowany. Liczba wezlow: " << power_grid.size() << std::endl;
  
  std::cout << "\n2. Algorytm Dijkstry - Znajdowanie najmniej stratnej sciezki przesylowej" << std::endl;
  size_t start_node = 0; // Elektrownia
  auto [distances, prev] = power_grid.dijkstra(start_node);

  for (size_t i = 1; i < wezly.size(); ++i) {
    std::cout << "   Minimalna strata energii do " << wezly[i] << ": " << distances[i] << " (jednostek oporu)" << std::endl;
    
    // Rekonstrukcja trasy
    std::vector<size_t> path;
    int curr = i;
    while (curr != -1) {
      path.push_back(curr);
      curr = prev[curr];
    }
    
    std::cout << "      Trasa przesyłu: ";
    for (int p = path.size() - 1; p >= 0; --p) {
      std::cout << wezly[path[p]] << (p > 0 ? " -> " : "");
    }
    std::cout << "\n" << std::endl;
  }

  return 0;
}
