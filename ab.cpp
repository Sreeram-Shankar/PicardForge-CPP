#include "ab.hpp"
#include <cmath>
#include <algorithm>

namespace ab {

//defines the step for AB2
void step_ab2(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, double h, int n,
              std::vector<double>& y_next) {
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y_n[i] + h * ((3.0/2.0)*f_n[i] - (1.0/2.0)*f_prev[i]);
    }
}

//main solver for AB2
void solve_ab2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    std::vector<double> k1(n), k2(n);
    
    f(tgrid[0], Yout[0], F_vals[0]);
    
    //uses RK2 for first step
    f(tgrid[0], Yout[0], k1);
    std::vector<double> y_temp(n);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = Yout[0][i] + h * k1[i];
    }
    f(tgrid[0] + h, y_temp, k2);
    for (int i = 0; i < n; ++i) {
        Yout[1][i] = Yout[0][i] + 0.5 * h * (k1[i] + k2[i]);
    }
    f(tgrid[1], Yout[1], F_vals[1]);
    
    for (int k = 1; k < nsteps; ++k) {
        step_ab2(Yout[k], F_vals[k], F_vals[k-1], h, n, Yout[k+1]);
        f(tgrid[k+1], Yout[k+1], F_vals[k+1]);
    }
}

//defines the step for AB3
void step_ab3(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, const std::vector<double>& f_prev2,
              double h, int n, std::vector<double>& y_next) {
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y_n[i] + h * ((23.0/12.0)*f_n[i] - (16.0/12.0)*f_prev[i] + (5.0/12.0)*f_prev2[i]);
    }
}

//main solver for AB3
void solve_ab3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    std::vector<double> k1(n), k2(n), k3(n);
    
    f(tgrid[0], Yout[0], F_vals[0]);
    
    //uses RK3 for first two steps
    f(tgrid[0], Yout[0], k1);
    std::vector<double> y_temp(n);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = Yout[0][i] + (h/2.0) * k1[i];
    }
    f(tgrid[0] + h/2.0, y_temp, k2);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = Yout[0][i] + h * (-k1[i] + 2.0*k2[i]);
    }
    f(tgrid[0] + h, y_temp, k3);
    for (int i = 0; i < n; ++i) {
        Yout[1][i] = Yout[0][i] + (h/6.0) * (k1[i] + 4.0*k2[i] + k3[i]);
    }
    f(tgrid[1], Yout[1], F_vals[1]);
    
    f(tgrid[1], Yout[1], k1);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = Yout[1][i] + (h/2.0) * k1[i];
    }
    f(tgrid[1] + h/2.0, y_temp, k2);
    for (int i = 0; i < n; ++i) {
        y_temp[i] = Yout[1][i] + h * (-k1[i] + 2.0*k2[i]);
    }
    f(tgrid[1] + h, y_temp, k3);
    for (int i = 0; i < n; ++i) {
        Yout[2][i] = Yout[1][i] + (h/6.0) * (k1[i] + 4.0*k2[i] + k3[i]);
    }
    f(tgrid[2], Yout[2], F_vals[2]);
    
    for (int k = 2; k < nsteps; ++k) {
        step_ab3(Yout[k], F_vals[k], F_vals[k-1], F_vals[k-2], h, n, Yout[k+1]);
        f(tgrid[k+1], Yout[k+1], F_vals[k+1]);
    }
}

//defines the step for AB4
void step_ab4(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, const std::vector<double>& f_prev2,
              const std::vector<double>& f_prev3, double h, int n,
              std::vector<double>& y_next) {
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y_n[i] + h * ((55.0/24.0)*f_n[i] - (59.0/24.0)*f_prev[i] +
                                   (37.0/24.0)*f_prev2[i] - (9.0/24.0)*f_prev3[i]);
    }
}

//main solver for AB4
void solve_ab4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    std::vector<double> k1(n), k2(n), k3(n), k4(n);
    
    f(tgrid[0], Yout[0], F_vals[0]);
    
    //uses RK4 for first three steps
    for (int i = 0; i < 3; ++i) {
        f(tgrid[i], Yout[i], k1);
        std::vector<double> y_temp(n);
        for (int j = 0; j < n; ++j) {
            y_temp[j] = Yout[i][j] + (h/2.0) * k1[j];
        }
        f(tgrid[i] + h/2.0, y_temp, k2);
        for (int j = 0; j < n; ++j) {
            y_temp[j] = Yout[i][j] + (h/2.0) * k2[j];
        }
        f(tgrid[i] + h/2.0, y_temp, k3);
        for (int j = 0; j < n; ++j) {
            y_temp[j] = Yout[i][j] + h * k3[j];
        }
        f(tgrid[i] + h, y_temp, k4);
        for (int j = 0; j < n; ++j) {
            Yout[i+1][j] = Yout[i][j] + (h/6.0) * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }
        f(tgrid[i+1], Yout[i+1], F_vals[i+1]);
    }
    
    for (int k = 3; k < nsteps; ++k) {
        step_ab4(Yout[k], F_vals[k], F_vals[k-1], F_vals[k-2], F_vals[k-3], h, n, Yout[k+1]);
        f(tgrid[k+1], Yout[k+1], F_vals[k+1]);
    }
}

//defines the step for AB5
void step_ab5(const std::vector<double>& y_n, const std::vector<double>& f_n,
              const std::vector<double>& f_prev, const std::vector<double>& f_prev2,
              const std::vector<double>& f_prev3, const std::vector<double>& f_prev4,
              double h, int n, std::vector<double>& y_next) {
    y_next.resize(n);
    for (int i = 0; i < n; ++i) {
        y_next[i] = y_n[i] + h * ((1901.0/720.0)*f_n[i] - (2774.0/720.0)*f_prev[i] +
                                   (2616.0/720.0)*f_prev2[i] - (1274.0/720.0)*f_prev3[i] +
                                   (251.0/720.0)*f_prev4[i]);
    }
}

//main solver for AB5
void solve_ab5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    std::vector<double> k1(n), k2(n), k3(n), k4(n);
    
    f(tgrid[0], Yout[0], F_vals[0]);
    
    //uses RK4 for first four steps
    for (int i = 0; i < 4; ++i) {
        f(tgrid[i], Yout[i], k1);
        std::vector<double> y_temp(n);
        for (int j = 0; j < n; ++j) {
            y_temp[j] = Yout[i][j] + (h/2.0) * k1[j];
        }
        f(tgrid[i] + h/2.0, y_temp, k2);
        for (int j = 0; j < n; ++j) {
            y_temp[j] = Yout[i][j] + (h/2.0) * k2[j];
        }
        f(tgrid[i] + h/2.0, y_temp, k3);
        for (int j = 0; j < n; ++j) {
            y_temp[j] = Yout[i][j] + h * k3[j];
        }
        f(tgrid[i] + h, y_temp, k4);
        for (int j = 0; j < n; ++j) {
            Yout[i+1][j] = Yout[i][j] + (h/6.0) * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }
        f(tgrid[i+1], Yout[i+1], F_vals[i+1]);
    }
    
    for (int k = 4; k < nsteps; ++k) {
        step_ab5(Yout[k], F_vals[k], F_vals[k-1], F_vals[k-2], F_vals[k-3], F_vals[k-4], h, n, Yout[k+1]);
        f(tgrid[k+1], Yout[k+1], F_vals[k+1]);
    }
}
}