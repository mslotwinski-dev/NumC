#include "../src/numc.h"
#include <iostream>
#include <random>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 11 Replication: Rozkład PCA ==========\n";

  int N = 500;
  vector<double> rx = random::normal(0.0, 5.0, N);
  vector<double> ry = random::normal(0.0, 1.0, N);
  
  // Transformacja obrotu o kąt theta = 45 stopni
  double theta = std::acos(-1.0) / 4.0;
  double cos_t = std::cos(theta);
  double sin_t = std::sin(theta);

  vector<double> X = rx * cos_t - ry * sin_t;
  vector<double> Y = rx * sin_t + ry * cos_t;

  std::vector<vector<double>> dataset = {X, Y};

  // 1. Wyznaczenie macierzy kowariancji
  tensor<double> cov = statistics::cov_matrix(dataset);
  std::cout << "Macierz Kowariancji:\n" << cov << "\n";

  // 2. Wartości i wektory własne
  // używając iteracji Jacobiego (zakładamy symetryczną macierz)
  auto eigens = linear_algebra::jacobi(cov);
  vector<double> ev = eigens.first;  // wartości własne
  tensor<double> V = eigens.second;  // wektory własne (kolumny)
  linear_algebra::sort_jacobi(ev, V); // Opcjonalne posortowanie


  std::cout << "Wartości własne: " << ev << "\n";
  std::cout << "Wektory własne (kolumny):\n" << V << "\n";

  // 3. Transformacja danych do bazy wektorów własnych (PCA)
  // X_pca = V^T * X_orig
  tensor<double> V_T = linear_algebra::transpose(V);
  
  vector<double> X_pca(N), Y_pca(N);
  for (int i = 0; i < N; ++i) {
    vector<double> point = {X[i], Y[i]};
    vector<double> projected = linear_algebra::matvec(V_T, point);
    X_pca[i] = projected[0];
    Y_pca[i] = projected[1];
  }

  // WIZUALIZACJA
  plot dashboard("ROOT Lab 11: Analiza Głównych Składowych (PCA)");
  dashboard.set_xlabel("x").set_ylabel("y");

  dashboard.add_scatter(X, Y, "Oryginalne zaszumione dane (skorelowane)");
  dashboard.add_scatter(X_pca, Y_pca, "Zrzutowane dane PCA (nieskorelowane)");

  // Linie głównych osi na z oryginalnymi danymi
  // PC1
  vector<double> pc1_x = {-ev[0]*V(0,0), ev[0]*V(0,0)};
  vector<double> pc1_y = {-ev[0]*V(1,0), ev[0]*V(1,0)};
  dashboard.add_scatter(pc1_x, pc1_y, "Główna oś 1 (PC1)");

  dashboard.show();

  return 0;
}
