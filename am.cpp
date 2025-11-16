#include "am.hpp"
#include <cmath>
#include <algorithm>

namespace am {

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM2
void solve_am2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps, double tol) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    std::vector<double> y(n), y_old(n), f_next(n);
    
    f(tgrid[0], Yout[0], F_vals[0]);
    
    //uses BE bootstrap
    f(tgrid[1], Yout[0], f_next);
    for (int i = 0; i < n; ++i) {
        Yout[1][i] = Yout[0][i] + h * f_next[i];
    }
    f(tgrid[1], Yout[1], F_vals[1]);
    
    for (int k = 1; k < nsteps; ++k) {
        y = Yout[k];
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_next);
            for (int j = 0; j < n; ++j) {
                y[j] = Yout[k][j] + h * (0.5*f_next[j] + 0.5*F_vals[k][j]);
            }
            double diff_norm = 0.0;
            for (int j = 0; j < n; ++j) {
                double diff = y[j] - y_old[j];
                diff_norm += diff * diff;
            }
            diff_norm = std::sqrt(diff_norm);
            if (diff_norm < tol) break;
        }
        
        Yout[k+1] = y;
        F_vals[k+1] = f_next;
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM3
void solve_am3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps, double tol) {
    
    //uses AM2 to bootstrap 1 step
    std::vector<std::vector<double>> Y_boot;
    std::vector<double> tgrid_boot;
    solve_am2(f, t0, t0+2.0*h, y0, h, n, Y_boot, tgrid_boot, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    for (int i = 0; i < 3 && i < static_cast<int>(Y_boot.size()); ++i) {
        Yout[i] = Y_boot[i];
    }
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    
    //compute initial F
    for (int i = 0; i < 3; ++i) {
        f(tgrid[i], Yout[i], F_vals[i]);
    }
    
    //defines the main AM3 solver
    std::vector<double> y(n), y_old(n), f_next(n);
    for (int k = 2; k < nsteps; ++k) {
        y = Yout[k];
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_next);
            for (int j = 0; j < n; ++j) {
                y[j] = Yout[k][j] + h * ((5.0/12.0)*f_next[j] + (2.0/3.0)*F_vals[k][j] - (1.0/12.0)*F_vals[k-1][j]);
            }
            double diff_norm = 0.0;
            for (int j = 0; j < n; ++j) {
                double diff = y[j] - y_old[j];
                diff_norm += diff * diff;
            }
            diff_norm = std::sqrt(diff_norm);
            if (diff_norm < tol) break;
        }
        
        Yout[k+1] = y;
        F_vals[k+1] = f_next;
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM4
void solve_am4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps, double tol) {
    
    //uses AM3 to bootstrap 2 steps
    std::vector<std::vector<double>> Y_boot;
    std::vector<double> tgrid_boot;
    solve_am3(f, t0, t0+3.0*h, y0, h, n, Y_boot, tgrid_boot, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    for (int i = 0; i < 4 && i < static_cast<int>(Y_boot.size()); ++i) {
        Yout[i] = Y_boot[i];
    }
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    
    //compute initial F
    for (int i = 0; i < 4; ++i) {
        f(tgrid[i], Yout[i], F_vals[i]);
    }
    
    //defines the main AM4 solver
    std::vector<double> y(n), y_old(n), f_next(n);
    for (int k = 3; k < nsteps; ++k) {
        y = Yout[k];
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_next);
            for (int j = 0; j < n; ++j) {
                y[j] = Yout[k][j] + h * ((3.0/8.0)*f_next[j] + (19.0/24.0)*F_vals[k][j] -
                                         (5.0/24.0)*F_vals[k-1][j] + (1.0/24.0)*F_vals[k-2][j]);
            }
            double diff_norm = 0.0;
            for (int j = 0; j < n; ++j) {
                double diff = y[j] - y_old[j];
                diff_norm += diff * diff;
            }
            diff_norm = std::sqrt(diff_norm);
            if (diff_norm < tol) break;
        }
        
        Yout[k+1] = y;
        F_vals[k+1] = f_next;
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation AM5
void solve_am5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
               double h, int n, std::vector<std::vector<double>>& Yout,
               std::vector<double>& tgrid, int sweeps, double tol) {
    
    //uses AM4 to bootstrap 3 steps of history
    std::vector<std::vector<double>> Y_boot;
    std::vector<double> tgrid_boot;
    solve_am4(f, t0, t0+4.0*h, y0, h, n, Y_boot, tgrid_boot, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    for (int i = 0; i < 5 && i < static_cast<int>(Y_boot.size()); ++i) {
        Yout[i] = Y_boot[i];
    }
    
    std::vector<std::vector<double>> F_vals(nsteps + 1, std::vector<double>(n));
    
    //computes initial F
    for (int i = 0; i < 5; ++i) {
        f(tgrid[i], Yout[i], F_vals[i]);
    }
    
    //defines the main AM5 solver
    std::vector<double> y(n), y_old(n), f_next(n);
    for (int k = 4; k < nsteps; ++k) {
        y = Yout[k];
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_next);
            for (int j = 0; j < n; ++j) {
                y[j] = Yout[k][j] + h * ((251.0/720.0)*f_next[j] + (646.0/720.0)*F_vals[k][j] -
                                         (264.0/720.0)*F_vals[k-1][j] + (106.0/720.0)*F_vals[k-2][j] -
                                         (19.0/720.0)*F_vals[k-3][j]);
            }
            double diff_norm = 0.0;
            for (int j = 0; j < n; ++j) {
                double diff = y[j] - y_old[j];
                diff_norm += diff * diff;
            }
            diff_norm = std::sqrt(diff_norm);
            if (diff_norm < tol) break;
        }
        
        Yout[k+1] = y;
        F_vals[k+1] = f_next;
    }
}
}
