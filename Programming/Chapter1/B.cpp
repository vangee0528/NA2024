#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>

const double Pi = acos(-1.);

class F1 : public Function {
public:
    double operator() (double x) const {
        return 1.0/x-tan(x);
    }
};

void solve_f1() {
    std::cout << "Solving x^{-1} - \\tan x on [0, \\pi/2]" << std::endl;
    Bisection_Method solver_f1(F1(), 0, Pi/2);
    double x = solver_f1.solve();
    std::cout << "A root is: " << x << std::endl;
}

/* Type your code here */

int main() {
    solve_f1();
    /* Type your code here */
    return 0;
}