#ifndef BDF_HPP
#define BDF_HPP

#include <vector>
#include <functional>

namespace bdf {

//defines the type alias for the ODE function: f(t, y, dydt)
using ODE_Func = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BE (BDF1)
void solve_be(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
              double h, int n, std::vector<std::vector<double>>& Yout,
              std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF2
void solve_bdf2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF3
void solve_bdf3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF4
void solve_bdf4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF5
void solve_bdf5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF6
void solve_bdf6(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps = 12, double tol = 1e-10);

}

#endif

