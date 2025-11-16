#include "sdirk.hpp"
#include <algorithm>
#include <cmath>

namespace sdirk {

//defines the SDIRK step with Gauss-Seidel relaxation
void step_sdirk(ODE_Func f, double t, const std::vector<double>& y, double h,
                const std::vector<std::vector<double>>& A,
                const std::vector<double>& b,
                const std::vector<double>& c,
                int s, int n, std::vector<double>& y_next,
                int sweeps, double tol) {
    
    //initializes all stages with the initial guess y
    std::vector<std::vector<double>> Y_stages(s, y);
    std::vector<std::vector<double>> Y_old(s, std::vector<double>(n));
    std::vector<double> rhs(n);
    std::vector<double> fval(n);
    
    //implements Gauss-Seidel relaxation
    for (int k = 0; k < sweeps; ++k) {
        Y_old = Y_stages;
        for (int i = 0; i < s; ++i) {
            std::fill(rhs.begin(), rhs.end(), 0.0);
            for (int j = 0; j < s; ++j) {
                f(t + c[j]*h, Y_stages[j], fval);
                for (int idx = 0; idx < n; ++idx) {
                    rhs[idx] += A[i][j] * fval[idx];
                }
            }
            for (int idx = 0; idx < n; ++idx) {
                Y_stages[i][idx] = y[idx] + h * rhs[idx];
            }
        }
        
        //computes L2 norm of all stage differences
        double diff_norm = 0.0;
        for (int i = 0; i < s; ++i) {
            for (int idx = 0; idx < n; ++idx) {
                double diff = Y_stages[i][idx] - Y_old[i][idx];
                diff_norm += diff * diff;
            }
        }
        diff_norm = std::sqrt(diff_norm);
        if (diff_norm < tol) break;
    }
    
    //computes the final state update
    y_next = y;
    for (int i = 0; i < s; ++i) {
        f(t + c[i]*h, Y_stages[i], fval);
        for (int idx = 0; idx < n; ++idx) {
            y_next[idx] += h * b[i] * fval[idx];
        }
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK2
void solve_sdirk2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                  double h, int n, std::vector<std::vector<double>>& Yout,
                  std::vector<double>& tgrid, int sweeps, double tol) {
    
    double gamma = 1.0 - 1.0/std::sqrt(2.0);
    std::vector<std::vector<double>> A = {{gamma, 0.0}, {1.0 - gamma, gamma}};
    std::vector<double> b = {1.0 - gamma, gamma};
    std::vector<double> c = {gamma, 1.0};
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    std::vector<double> y = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_sdirk(f, tgrid[k], y, h, A, b, c, 2, n, Yout[k+1], sweeps, tol);
        y = Yout[k+1];
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK3
void solve_sdirk3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                  double h, int n, std::vector<std::vector<double>>& Yout,
                  std::vector<double>& tgrid, int sweeps, double tol) {
    
    double gamma = 0.435866521508459;
    std::vector<std::vector<double>> A = {
        {gamma, 0.0, 0.0},
        {0.2820667395, gamma, 0.0},
        {1.208496649, -0.644363171, gamma}
    };
    std::vector<double> b = {1.208496649, -0.644363171, gamma};
    std::vector<double> c = {gamma, 0.7179332605, 1.0};
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    std::vector<double> y = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_sdirk(f, tgrid[k], y, h, A, b, c, 3, n, Yout[k+1], sweeps, tol);
        y = Yout[k+1];
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK4
void solve_sdirk4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                  double h, int n, std::vector<std::vector<double>>& Yout,
                  std::vector<double>& tgrid, int sweeps, double tol) {
    
    double gamma = 0.572816062482135;
    std::vector<std::vector<double>> A = {
        {gamma, 0.0, 0.0, 0.0},
        {-0.6557110092, gamma, 0.0, 0.0},
        {0.757184241, 0.237758128, gamma, 0.0},
        {0.155416858, 0.701913790, 0.142669351, gamma}
    };
    std::vector<double> b = {0.155416858, 0.701913790, 0.142669351, gamma};
    std::vector<double> c = {gamma, 0.344, 0.995, 1.0};
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    std::vector<double> y = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_sdirk(f, tgrid[k], y, h, A, b, c, 4, n, Yout[k+1], sweeps, tol);
        y = Yout[k+1];
    }
}
}