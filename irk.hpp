#ifndef IRK_HPP
#define IRK_HPP

#include <vector>
#include <functional>
#include <string>

namespace irk {

//defines the type alias for the ODE function: f(t, y, dydt)
using ODE_Func = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;

//loads the Butcher tableaus for the given family and stage count
void get_tableau(const std::string& family, int s,
                 std::vector<std::vector<double>>& A,
                 std::vector<double>& b,
                 std::vector<double>& c);

//defines the IRK step with Gauss-Seidel relaxation
void step_collocation(ODE_Func f, double t, const std::vector<double>& y, double h,
                      const std::vector<std::vector<double>>& A,
                      const std::vector<double>& b,
                      const std::vector<double>& c,
                      int s, int n, std::vector<double>& y_next,
                      int sweeps = 12, double tol = 1e-10);

//main solver for any collocation IRK method using Gauss-Seidel relaxation
void solve_collocation(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                       double h, int n, std::vector<std::vector<double>>& Yout,
                       std::vector<double>& tgrid, const std::string& family, int s,
                       int sweeps = 12, double tol = 1e-10);

}

#endif
