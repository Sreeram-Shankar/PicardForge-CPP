#include "irk.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace irk {

//function that converts a string to lowercase
std::string to_lower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

//loads the Butcher tableaus for the given family and stage count
void get_tableau(const std::string& family, int s,
                 std::vector<std::vector<double>>& A,
                 std::vector<double>& b,
                 std::vector<double>& c) {
    
    std::string family_lower = to_lower(family);
    
    if (family_lower.substr(0, 5) == "gauss") {
        switch (s) {
            case 1: {
                A = {{0.5}};
                b = {1.0};
                c = {0.5};
                break;
            }
            case 2: {
                A = {{0.25, 0.5386751345948129},
                     {-0.03867513459481288, 0.25}};
                b = {0.5, 0.5};
                c = {0.7886751345948129, 0.21132486540518712};
                break;
            }
            case 3: {
                A = {{0.1388888888888889, 0.48042111196938335, 0.26798833376246945},
                     {-0.022485417203086815, 0.2222222222222222, 0.3002631949808646},
                     {0.009789444015308326, -0.03597666752493890, 0.1388888888888889}};
                b = {0.2777777777777778, 0.4444444444444444, 0.2777777777777778};
                c = {0.8872983346207417, 0.5, 0.1127016653792583};
                break;
            }
            case 4: {
                A = {{0.08696371128436346, 0.35267675751627186, 0.31344511474186835, 0.1774825722545226},
                     {-0.014190694931141143, 0.16303628871563654, 0.3539530060337440, 0.16719192197418877},
                     {0.006735500594538155, -0.027880428602470895, 0.16303628871563654, 0.18811811749986807},
                     {-0.003555149685795683, 0.012627462689404725, -0.026604180084998793, 0.08696371128436346}};
                b = {0.17392742256872693, 0.32607257743127307, 0.32607257743127307, 0.17392742256872693};
                c = {0.9305681557970263, 0.6699905217924281, 0.33000947820757187, 0.06943184420297371};
                break;
            }
            case 5: {
                A = {{0.05923172126404727, 0.25888469960875927, 0.2731900436258015, 0.24490812891049542, 0.11687532956022855},
                     {-0.00968756314195074, 0.11965716762484162, 0.30903655906408665, 0.22899605457899988, 0.12123243692686415},
                     {0.004687154523869941, -0.020690316430958285, 0.14222222222222222, 0.2600046516806415, 0.11377628800422460},
                     {-0.002768994398769603, 0.010318280670683357, -0.024592114619642200, 0.11965716762484162, 0.12815100567004528},
                     {0.001588112967865998, -0.005593793660812185, 0.011254400818642956, -0.019570364359076037, 0.05923172126404727}};
                b = {0.11846344252809454, 0.23931433524968323, 0.28444444444444444, 0.23931433524968323, 0.11846344252809454};
                c = {0.9530899229693320, 0.7692346550528416, 0.5, 0.23076534494715845, 0.04691007703066800};
                break;
            }
            default:
                throw std::runtime_error("Gauss-Legendre only implemented for s = 1...5");
        }
    }
    else if (family_lower.substr(0, 5) == "radau") {
        switch (s) {
            case 2: {
                A = {{0.41666666666666667, -0.08333333333333333},
                     {0.75, 0.25}};
                b = {0.75, 0.25};
                c = {0.3333333333333333, 1.0};
                break;
            }
            case 3: {
                A = {{0.19681547722366043, -0.06553542585019839, 0.023770974348220152},
                     {0.39442431473908727, 0.29207341166522846, -0.04154875212599793},
                     {0.37640306270046727, 0.5124858261884216, 0.11111111111111111}};
                b = {0.37640306270046727, 0.5124858261884216, 0.11111111111111111};
                c = {0.15505102572168219, 0.6449489742783178, 1.0};
                break;
            }
            case 4: {
                A = {{0.20689257393535890, 0.23438399574740026, -0.04785712804854072, 0.01604742280651627},
                     {-0.04030922072352221, 0.11299947932315619, 0.02580237742033639, -0.009904676507266424},
                     {0.4061232638673733, 0.21668178462325034, 0.18903651817005634, -0.02418210489983294},
                     {0.38819346884317188, 0.22046221117676838, 0.32884431998005974, 0.06250000000000000}};
                b = {0.38819346884317188, 0.22046221117676838, 0.32884431998005974, 0.0625};
                c = {0.4094668644407347, 0.08858795951270395, 0.7876594617608471, 1.0};
                break;
            }
            case 5: {
                A = {{0.14621486784749350, 0.15377523147918247, -0.03644456890512809, 0.02123306311930472, -0.007935579902728778},
                     {-0.02673533110794557, 0.07299886431790332, 0.01867692976398435, -0.01287910609330644, 0.005042839233882015},
                     {0.29896712949128348, 0.14006304568480987, 0.16758507013524896, -0.03396910168661775, 0.010944288744192252},
                     {0.27650006876015923, 0.14489430810953476, 0.32579792291042103, 0.12875675325490976, -0.015708917378805328},
                     {0.28135601514946206, 0.14371356079122594, 0.31182652297574125, 0.22310390108357074, 0.04000000000000000}};
                b = {0.28135601514946206, 0.14371356079122594, 0.31182652297574125, 0.22310390108357074, 0.04};
                c = {0.27684301363812383, 0.05710419611451768, 0.5835904323689168, 0.8602401356562194, 1.0};
                break;
            }
            default:
                throw std::runtime_error("Radau IIA only implemented for s = 2...5");
        }
    }
    else if (family_lower.substr(0, 7) == "lobatto") {
        switch (s) {
            case 2: {
                A = {{0.0, 0.0},
                     {0.5, 0.5}};
                b = {0.5, 0.5};
                c = {0.0, 1.0};
                break;
            }
            case 3: {
                A = {{0.0, 0.0, 0.0},
                     {0.20833333333333333, 0.3333333333333333, -0.041666666666666664},
                     {0.16666666666666667, 0.6666666666666666, 0.16666666666666666}};
                b = {0.16666666666666666, 0.6666666666666666, 0.16666666666666666};
                c = {0.0, 0.5, 1.0};
                break;
            }
            case 4: {
                A = {{0.0, 0.0, 0.0, 0.0},
                     {0.07303276685416842, 0.22696723314583158, 0.45057403089581055, -0.02696723314583158},
                     {0.11030056647916491, -0.03390736422914388, 0.18969943352083508, 0.01030056647916491},
                     {0.08333333333333333, 0.4166666666666667, 0.4166666666666667, 0.08333333333333333}};
                b = {0.08333333333333333, 0.4166666666666667, 0.4166666666666667, 0.08333333333333333};
                c = {0.0, 0.7236067977499790, 0.2763932022500210, 1.0};
                break;
            }
            case 5: {
                A = {{0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.05370013924241453, 0.15247745287881054, 0.37729127742211367, 0.26158639799680673, -0.017728432186156897},
                     {0.040625, -0.030961961100820556, 0.17777777777777778, 0.30318418332304278, 0.009375},
                     {0.06772843218615690, 0.01063582422541549, -0.021735721866558114, 0.11974476934341168, -0.003700139242414531},
                     {0.05, 0.2722222222222222, 0.35555555555555557, 0.2722222222222222, 0.05}};
                b = {0.05, 0.2722222222222222, 0.35555555555555557, 0.2722222222222222, 0.05};
                c = {0.0, 0.8273268353539886, 0.5, 0.1726731646460114, 1.0};
                break;
            }
            default:
                throw std::runtime_error("Lobatto IIIC only implemented for s = 2...5");
        }
    }
    else {
        throw std::runtime_error("Unknown IRK family '" + family + "'. Must be 'gauss', 'radau', or 'lobatto'.");
    }
}

//defines the IRK step with Gauss-Seidel relaxation
void step_collocation(ODE_Func f, double t, const std::vector<double>& y, double h,
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

//main solver for any collocation IRK method using Gauss-Seidel relaxation
void solve_collocation(ODE_Func f, double t0, double tf, const std::vector<double>& y0,
                       double h, int n, std::vector<std::vector<double>>& Yout,
                       std::vector<double>& tgrid, const std::string& family, int s,
                       int sweeps, double tol) {
    
    std::vector<std::vector<double>> A;
    std::vector<double> b, c;
    
    //gets the Butcher tableau
    get_tableau(family, s, A, b, c);
    
    int nsteps = static_cast<int>(std::ceil((tf - t0)/h));
    
    tgrid.resize(nsteps + 1);
    for (int k = 0; k <= nsteps; ++k) {
        tgrid[k] = t0 + k * h;
    }
    
    Yout.resize(nsteps + 1, std::vector<double>(n));
    Yout[0] = y0;
    std::vector<double> y = y0;
    
    for (int k = 0; k < nsteps; ++k) {
        step_collocation(f, tgrid[k], y, h, A, b, c, s, n, Yout[k+1], sweeps, tol);
        y = Yout[k+1];
    }
}
}