#ifndef AB_HPP
#define AB_HPP

#include <vector>
#include <functional>

namespace ab {

//defines the type alias for the ODE function: f(t, y, dydt)
using ODE_Func = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;

//defines the step for AB2
void step_ab2(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, double h, int n,
              std::vector<double>& y_next);

//main solver for AB2
void solve_ab2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid);

//defines the step for AB3
void step_ab3(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, const std::vector<double>& f_prev2,
              double h, int n, std::vector<double>& y_next);

//main solver for AB3
void solve_ab3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid);

//defines the step for AB4
void step_ab4(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, const std::vector<double>& f_prev2,
              const std::vector<double>& f_prev3, double h, int n,
              std::vector<double>& y_next);

//main solver for AB4
void solve_ab4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid);

//defines the step for AB5
void step_ab5(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, const std::vector<double>& f_prev2,
              const std::vector<double>& f_prev3, const std::vector<double>& f_prev4,
              double h, int n, std::vector<double>& y_next);

//main solver for AB5
void solve_ab5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid);

}

#endif