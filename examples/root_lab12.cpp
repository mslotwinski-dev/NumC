#include "../src/numc.h"
#include <iostream>
#include <random>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 12 Replication: Machine Learning (K-NN) ==========\n";

  int N_per_class = 50;
  
  vector<double> class1_x = random::normal(2.0, 1.0, N_per_class);
  vector<double> class1_y = random::normal(2.0, 1.0, N_per_class);
  vector<double> class2_x = random::normal(-2.0, 1.5, N_per_class);
  vector<double> class2_y = random::normal(-2.0, 1.5, N_per_class);

  std::vector<vector<double>> X_train;
  std::vector<double> y_train;

  for (int i = 0; i < N_per_class; ++i) {
    X_train.push_back({class1_x[i], class1_y[i]});
    y_train.push_back(1.0);
    
    X_train.push_back({class2_x[i], class2_y[i]});
    y_train.push_back(2.0);
  }

  // Trenowanie klasyfikatora K-NN
  machine_learning::KNNClassifier<double> knn(5); // K = 5
  knn.train(X_train, y_train);

  // Generowanie siatki "decision boundary" na plot
  int grid_size = 50;
  double min_val = -6.0, max_val = 6.0;
  double step = (max_val - min_val) / (grid_size - 1);

  vector<double> bg_c1_x, bg_c1_y;
  vector<double> bg_c2_x, bg_c2_y;

  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      double px = min_val + i * step;
      double py = min_val + j * step;
      vector<double> pt = {px, py};
      double pred = knn.predict(pt);
      
      if (pred == 1.0) {
        bg_c1_x.push_back(px);
        bg_c1_y.push_back(py);
      } else {
        bg_c2_x.push_back(px);
        bg_c2_y.push_back(py);
      }
    }
  }

  // WIZUALIZACJA
  plot dashboard("ROOT Lab 12: Decision Boundary (K-NN)");
  dashboard.set_xlabel("x").set_ylabel("y");

  // Rysujemy tło (bardzo drobne punkty, dla klastra 1 i 2)
  dashboard.add_scatter(bg_c1_x, bg_c1_y, "Decision Region 1");
  dashboard.add_scatter(bg_c2_x, bg_c2_y, "Decision Region 2");

  // Rysujemy rzeczywiste punkty na wierzchu
  dashboard.add_scatter(class1_x, class1_y, "Klasa 1 (Train)");
  dashboard.add_scatter(class2_x, class2_y, "Klasa 2 (Train)");

  dashboard.show();

  return 0;
}
