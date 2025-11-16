#ifndef AM_HPP
#define AM_HPP

#include <vector>
#include <functional>

namespace am {

//defines the type alias for the ODE function: f(t, y, dydt)
using ODE_Func = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM2
void solve_am2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM3
void solve_am3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM4
void solve_am4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM5
void solve_am5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

}

#endif