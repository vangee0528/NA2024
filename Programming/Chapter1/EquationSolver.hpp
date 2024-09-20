#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"

class EquationSolver{
protected:
    const Function & F;
public:
    EquationSolver(const Function& F) : F(F) {}
    virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver {
private:
    double a, b;
    double eps, delta;
    int Maxiter;
public:
    Bisection_Method(const Function &F, double a, double b, 
        double eps = 1e-7, double delta = 1e-6, int Maxiter = 50) :
        EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}
    
    virtual double solve() {
        /* Type your code here. */
        return 0;   // to be replaced.
    }
};

class Newton_Method : public EquationSolver {
private:
    double x0;
    double eps;
    int Maxiter;
public:
    Newton_Method(const Function &F, double x0, 
        double eps = 1e-7, double Maxiter = 8) :
        EquationSolver(F), x0(x0), Maxiter(Maxiter), eps(eps) {}
    
    virtual double solve() {
        /* Type your code here. */
        return 0;   // to be replaced.
    }
};

/* Type your code here. */

#endif