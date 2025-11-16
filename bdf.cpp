#include "bdf.hpp"
#include <cmath>
#include <algorithm>

namespace bdf {

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BE (BDF1)
void solve_be(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
              double h, int n, std::vector<std::vector<double>>& Yout,
              std::vector<double>& tgrid, int sweeps, double tol) {
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    
    std::vector<double> y(n), y_old(n), f_val(n);
    
    //defines the main BE solver
    for (int k = 0; k < nsteps; ++k) {
        y = Yout[k];
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_val);
            for (int j = 0; j < n; ++j) {
                y[j] = Yout[k][j] + h * f_val[j];
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
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF2
void solve_bdf2(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps, double tol) {
    
    //uses backward euler to bootstrap
    std::vector<std::vector<double>> Y_BE;
    std::vector<double> tgrid_BE;
    solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    if (Y_BE.size() > 1) {
        Yout[1] = Y_BE[1];
    }
    
    std::vector<double> y(n), y_old(n), f_val(n), rhs(n);
    
    for (int k = 1; k < nsteps; ++k) {
        y = Yout[k];
        
        //defines the rhs
        for (int j = 0; j < n; ++j) {
            rhs[j] = (-4.0*Yout[k][j] + Yout[k-1][j]) / (2.0*h);
        }
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_val);
            for (int j = 0; j < n; ++j) {
                y[j] = (rhs[j] + f_val[j]) / (3.0/(2.0*h));
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
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF3
void solve_bdf3(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps, double tol) {
    
    //bootstrap using BE then BDF2
    std::vector<std::vector<double>> Y_BE, Y_BDF2;
    std::vector<double> tgrid_BE, tgrid_BDF2;
    solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol);
    solve_bdf2(f, t0, t0+2.0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    if (Y_BE.size() > 1) {
        Yout[1] = Y_BE[1];
    }
    if (Y_BDF2.size() > 2) {
        Yout[2] = Y_BDF2[2];
    }
    
    std::vector<double> y(n), y_old(n), f_val(n), rhs(n);
    
    for (int k = 2; k < nsteps; ++k) {
        y = Yout[k];
        
        //defines the rhs
        for (int j = 0; j < n; ++j) {
            rhs[j] = (-11.0*Yout[k][j] + 18.0*Yout[k-1][j] - 9.0*Yout[k-2][j] + 2.0*Yout[k-3][j]) / (6.0*h);
        }
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_val);
            for (int j = 0; j < n; ++j) {
                y[j] = (rhs[j] + f_val[j]) / (11.0/(6.0*h));
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
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF4
void solve_bdf4(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps, double tol) {
    
    //bootstrap using BDF1–3
    std::vector<std::vector<double>> Y_BE, Y_BDF2, Y_BDF3;
    std::vector<double> tgrid_BE, tgrid_BDF2, tgrid_BDF3;
    solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol);
    solve_bdf2(f, t0, t0+2.0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol);
    solve_bdf3(f, t0, t0+3.0*h, y0, h, n, Y_BDF3, tgrid_BDF3, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    if (Y_BE.size() > 1) {
        Yout[1] = Y_BE[1];
    }
    if (Y_BDF2.size() > 2) {
        Yout[2] = Y_BDF2[2];
    }
    if (Y_BDF3.size() > 3) {
        Yout[3] = Y_BDF3[3];
    }
    
    std::vector<double> y(n), y_old(n), f_val(n), rhs(n);
    
    for (int k = 3; k < nsteps; ++k) {
        y = Yout[k];
        
        //defines the rhs
        for (int j = 0; j < n; ++j) {
            rhs[j] = (-25.0*Yout[k][j] + 48.0*Yout[k-1][j] - 36.0*Yout[k-2][j] +
                     16.0*Yout[k-3][j] - 3.0*Yout[k-4][j]) / (12.0*h);
        }
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_val);
            for (int j = 0; j < n; ++j) {
                y[j] = (rhs[j] + f_val[j]) / (25.0/(12.0*h));
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
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF5
void solve_bdf5(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps, double tol) {
    
    //bootstrap up to BDF4
    std::vector<std::vector<double>> Y_BE, Y_BDF2, Y_BDF3, Y_BDF4;
    std::vector<double> tgrid_BE, tgrid_BDF2, tgrid_BDF3, tgrid_BDF4;
    solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol);
    solve_bdf2(f, t0, t0+2.0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol);
    solve_bdf3(f, t0, t0+3.0*h, y0, h, n, Y_BDF3, tgrid_BDF3, sweeps, tol);
    solve_bdf4(f, t0, t0+4.0*h, y0, h, n, Y_BDF4, tgrid_BDF4, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    if (Y_BE.size() > 1) {
        Yout[1] = Y_BE[1];
    }
    if (Y_BDF2.size() > 2) {
        Yout[2] = Y_BDF2[2];
    }
    if (Y_BDF3.size() > 3) {
        Yout[3] = Y_BDF3[3];
    }
    if (Y_BDF4.size() > 4) {
        Yout[4] = Y_BDF4[4];
    }
    
    std::vector<double> y(n), y_old(n), f_val(n), rhs(n);
    
    for (int k = 4; k < nsteps; ++k) {
        y = Yout[k];
        
        //defines the rhs
        for (int j = 0; j < n; ++j) {
            rhs[j] = (-137.0*Yout[k][j] + 300.0*Yout[k-1][j] - 300.0*Yout[k-2][j] +
                     200.0*Yout[k-3][j] - 75.0*Yout[k-4][j] + 12.0*Yout[k-5][j]) / (60.0*h);
        }
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_val);
            for (int j = 0; j < n; ++j) {
                y[j] = (rhs[j] + f_val[j]) / (137.0/(60.0*h));
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
    }
}

//solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF6
void solve_bdf6(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                double h, int n, std::vector<std::vector<double>>& Yout,
                std::vector<double>& tgrid, int sweeps, double tol) {
    
    //bootstrap up to BDF5
    std::vector<std::vector<double>> Y_BE, Y_BDF2, Y_BDF3, Y_BDF4, Y_BDF5;
    std::vector<double> tgrid_BE, tgrid_BDF2, tgrid_BDF3, tgrid_BDF4, tgrid_BDF5;
    solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol);
    solve_bdf2(f, t0, t0+2.0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol);
    solve_bdf3(f, t0, t0+3.0*h, y0, h, n, Y_BDF3, tgrid_BDF3, sweeps, tol);
    solve_bdf4(f, t0, t0+4.0*h, y0, h, n, Y_BDF4, tgrid_BDF4, sweeps, tol);
    solve_bdf5(f, t0, t0+5.0*h, y0, h, n, Y_BDF5, tgrid_BDF5, sweeps, tol);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    if (Y_BE.size() > 1) {
        Yout[1] = Y_BE[1];
    }
    if (Y_BDF2.size() > 2) {
        Yout[2] = Y_BDF2[2];
    }
    if (Y_BDF3.size() > 3) {
        Yout[3] = Y_BDF3[3];
    }
    if (Y_BDF4.size() > 4) {
        Yout[4] = Y_BDF4[4];
    }
    if (Y_BDF5.size() > 5) {
        Yout[5] = Y_BDF5[5];
    }
    
    std::vector<double> y(n), y_old(n), f_val(n), rhs(n);
    
    for (int k = 5; k < nsteps; ++k) {
        y = Yout[k];
        
        //defines the rhs
        for (int j = 0; j < n; ++j) {
            rhs[j] = (-147.0*Yout[k][j] + 360.0*Yout[k-1][j] - 450.0*Yout[k-2][j] +
                     400.0*Yout[k-3][j] - 225.0*Yout[k-4][j] + 72.0*Yout[k-5][j] -
                     10.0*Yout[k-6][j]) / (60.0*h);
        }
        
        //implements Gauss-Seidel relaxation
        for (int i = 0; i < sweeps; ++i) {
            y_old = y;
            f(tgrid[k+1], y, f_val);
            for (int j = 0; j < n; ++j) {
                y[j] = (rhs[j] + f_val[j]) / (147.0/(60.0*h));
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
    }
}
}