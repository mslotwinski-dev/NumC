#include "../src/numc.h"
#include <iostream>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 10 Replication: Równania Różniczkowe ==========\n";

  // System Oscylatora Tłumionego: 
  // d^2x/dt^2 + 2*gamma*dx/dt + omega^2*x = 0
  // Stan y = [x, v]
  // dy/dt = [v, -2*gamma*v - omega^2*x]

  double gamma = 0.15;
  double omega = 1.0;

  auto oscillator_deriv = [gamma, omega](double t, const vector<double>& y) {
    vector<double> dy(2);
    dy[0] = y[1];
    dy[1] = -2.0 * gamma * y[1] - omega * omega * y[0];
    return dy;
  };

  // Warunki początkowe
  double t_start = 0.0;
  double t_stop = 20.0;
  double h = 0.05;
  vector<double> y0 = {1.0, 0.0}; // Wychylenie początkowe x=1.0, prędkość v=0.0

  // Rozwiązywanie metodą Runge-Kutta 4 rzędu (API numc::differential_equations::run_kut4)
  auto result = differential_equations::run_kut4<double>(oscillator_deriv, t_start, y0, t_stop, h);
  
  vector<double> T_out = result.first;
  std::vector<vector<double>> Y = result.second;

  vector<double> x_pos(T_out.size());
  vector<double> v_vel(T_out.size());
  for (size_t i = 0; i < T_out.size(); ++i) {
    x_pos[i] = Y[i][0];
    v_vel[i] = Y[i][1];
  }

  // WIZUALIZACJA
  plot dashboard("ROOT Lab 10: Tłumiony Oscylator Harmoniczny (RK4)");
  dashboard.set_xlabel("Czas [t]").set_ylabel("Amplituda");
  dashboard.add_scatter(T_out, x_pos, "Położenie x(t)");
  dashboard.add_scatter(T_out, v_vel, "Prędkość v(t)");

  // Drugi wykres (Phase Space)
  plot phase_space("ROOT Lab 10: Portret fazowy (Phase Space)");
  phase_space.set_xlabel("x(t)").set_ylabel("v(t)");
  phase_space.add_scatter(x_pos, v_vel, "Trajektoria fazowa");

  dashboard.show();
  phase_space.show();

  return 0;
}
