#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::linear_algebra;

int main() {
  std::cout << "=== Projekt Linear Algebra: Kwantowy Oscylator Harmoniczny ===\n" << std::endl;

  // W modelowaniu oscylatora harmonicznego w 1D metodą Różnic Skończonych
  // Równanie Schrödingera: - (hbar^2 / 2m) d^2psi/dx^2 + V(x)psi = E*psi
  // Uprościmy: hbar = 1, m = 1 -> -0.5 * d^2psi/dx^2 + 0.5 * x^2 * psi = E*psi
  
  int N = 5; // Dla małego przykładu wymiar macierzy N x N
  double dx = 1.0; // Uproszczony krok siatki przestrzennej

  tensor<double> H({(size_t)N, (size_t)N});

  // Tworzenie macierzy Hamiltoniana H (T + V)
  // Elementy kinetyczne (T) z centralnej różnicy skończonej i potencjalne (V) na diagonali
  for (int i = 0; i < N; ++i) {
    double x = (i - N/2.0) * dx;
    double V_x = 0.5 * x * x;
    
    // Diagonala (energia potencjalna + 2 z pochodnej z drugiego rzędu)
    H(i, i) = V_x + (1.0 / (dx * dx));
    
    // Pozadiagonalne (kinetyczne z pochodnej drugiego rzędu)
    if (i > 0) H(i, i - 1) = -0.5 / (dx * dx);
    if (i < N - 1) H(i, i + 1) = -0.5 / (dx * dx);
  }

  std::cout << "1. Skonstruowano macierz Hamiltoniana H:\n" << H << std::endl;

  std::cout << "2. Obliczanie poziomów energetycznych (Eigenvalues) za pomocą algorytmu Jacobiego..." << std::endl;
  
  auto [eigenvalues, eigenvectors] = jacobi(H);
  
  // Sortowanie wartości własnych od najniższej energii (stan podstawowy)
  sort_jacobi(eigenvalues, eigenvectors);

  std::cout << "\n3. Dozwolone stany energetyczne (Eigenvalues) E_n:" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << "   Poziom E_" << i << " = " << std::fixed << std::setprecision(4) << eigenvalues[i] 
              << "   (Teoretycznie: ~" << (i + 0.5) << " w jednostkach hbar*omega)" << std::endl;
  }

  return 0;
}
