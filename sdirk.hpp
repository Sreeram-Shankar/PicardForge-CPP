#ifndef SDIRK_HPP
#define SDIRK_HPP

#include <vector>
#include <functional>
#include <cmath>

namespace sdirk {

//defines the type alias for the ODE function: f(t, y, dydt)
using ODE_Func = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;

//defines the SDIRK step with Gauss-Seidel relaxation
void step_sdirk(ODE_Func f, double t, const std::vector<double>& y, double h,
                const std::vector<std::vector<double>>& A,
                const std::vector<double>& b,
                const std::vector<double>& c,
                int s, int n, std::vector<double>& y_next,
                int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK2
void solve_sdirk2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                  double h, int n, std::vector<std::vector<double>>& Yout,
                  std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK3
void solve_sdirk3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                  double h, int n, std::vector<std::vector<double>>& Yout,
                  std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK4
void solve_sdirk4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                  double h, int n, std::vector<std::vector<double>>& Yout,
                  std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

}

#endif