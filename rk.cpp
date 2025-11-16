#include "rk.hpp"
#include <cmath>
#include <algorithm>

namespace rk {

//defines the step for RK1
void step_rk1(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next) {
    std::vector<double> k1(n);
    f(t, y, k1);
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y[i] + h * k1[i];
    }
}

//main solver for RK1
void solve_rk1(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_rk1(f, tgrid[k], Yout[k], h, n, Yout[k+1]);
    }
}

//defines the step for RK2
void step_rk2(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next) {
    std::vector<double> k1(n), k2(n);
    f(t, y, k1);
    std::vector<double> y_temp(n);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + h * k1[i];
    }
    f(t + h, y_temp, k2);
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y[i] + (h / 2.0) * (k1[i] + k2[i]);
    }
}

//main solver for RK2
void solve_rk2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_rk2(f, tgrid[k], Yout[k], h, n, Yout[k+1]);
    }
}

//defines the step for RK3
void step_rk3(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next) {
    std::vector<double> k1(n), k2(n), k3(n);
    f(t, y, k1);
    std::vector<double> y_temp(n);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/2.0) * k1[i];
    }
    f(t + h/2.0, y_temp, k2);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + h * (-k1[i] + 2.0*k2[i]);
    }
    f(t + h, y_temp, k3);
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y[i] + (h / 6.0) * (k1[i] + 4.0*k2[i] + k3[i]);
    }
}

//main solver for RK3
void solve_rk3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_rk3(f, tgrid[k], Yout[k], h, n, Yout[k+1]);
    }
}

//defines the step for RK4
void step_rk4(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next) {
    std::vector<double> k1(n), k2(n), k3(n), k4(n);
    f(t, y, k1);
    std::vector<double> y_temp(n);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/2.0) * k1[i];
    }
    f(t + h/2.0, y_temp, k2);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/2.0) * k2[i];
    }
    f(t + h/2.0, y_temp, k3);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + h * k3[i];
    }
    f(t + h, y_temp, k4);
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}

//main solver for RK4
void solve_rk4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_rk4(f, tgrid[k], Yout[k], h, n, Yout[k+1]);
    }
}

//defines the step for RK5
void step_rk5(ODE_Func f, double t, const std::vector<double>& y, double h, int n,
              std::vector<double>& y_next) {
    std::vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n);
    f(t, y, k1);
    std::vector<double> y_temp(n);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/4.0) * k1[i];
    }
    f(t + h/4.0, y_temp, k2);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/8.0) * (k1[i] + k2[i]);
    }
    f(t + h/4.0, y_temp, k3);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/2.0) * k3[i];
    }
    f(t + h/2.0, y_temp, k4);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/16.0) * (3.0*k1[i] + 9.0*k4[i]);
    }
    f(t + 3.0*h/4.0, y_temp, k5);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + (h/7.0) * (2.0*k1[i] + 3.0*k2[i] + 4.0*k4[i] - 12.0*k3[i]);
    }
    f(t + h, y_temp, k6);
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y[i] + (h / 90.0) * (7.0*k1[i] + 32.0*k3[i] + 12.0*k4[i] + 32.0*k5[i] + 7.0*k6[i]);
    }
}

//main solver for RK5
void solve_rk5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_rk5(f, tgrid[k], Yout[k], h, n, Yout[k+1]);
    }
}
}