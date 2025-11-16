#ifndef RK_HPP
#define RK_HPP

#include <vector>
#include <functional>

namespace rk {

//defines the type alias for the ODE function: f(t, y, dydt)    
using ODE_Func = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;

//defines the step for RK1
void step_rk1(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next);

//main solver for RK1
void solve_rk1(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid);

//defines the step for RK2
void step_rk2(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next);

//main solver for RK2
void solve_rk2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid);

//defines the step for RK3
void step_rk3(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next);

//main solver for RK3
void solve_rk3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid);

//defines the step for RK4
void step_rk4(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next);

//main solver for RK4
void solve_rk4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid);

//defines the step for RK5
void step_rk5(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next);

//main solver for RK5
void solve_rk5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid);

}

#endif