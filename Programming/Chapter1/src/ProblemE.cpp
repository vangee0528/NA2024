#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Function.hpp"          
#include "EquationSolver.hpp"    

// 定义函数 f(h) = -asin(h) + 0.5 * pi - h * sqrt(1 - h * h) - 1.24
class WaterDepthFunction : public Function {
public:
    virtual double operator() (double h) const override {
        return -asin(h) + 0.5 * M_PI - h * sqrt(1 - h * h) - 1.24;
    }

    virtual double derivative(double h) const override {
        double sqrt_term = std::sqrt(1 - h * h);
        return -1 / sqrt_term - sqrt_term - h * h / sqrt_term;
    }
};

int main() {
    try {
        WaterDepthFunction function;
        std::cout << "-----------PROBLEM E---------------" << std::endl << std::endl;
        
        // 1. 二分法求解
        std::cout << "Solving the equation with \033[1;32m(Bisection Method)\033[0m" << std::endl;
        Bisection_Method bisection_solver(function, 0, 1);
        double root_bisection = bisection_solver.solve();
        std::cout << "\033[32mRoot: " << root_bisection << " ft\033[0m" << std::endl << std::endl;

        // 2. 牛顿法求解
        std::cout << "--------------------------------" << std::endl << std::endl;
        std::cout << "Solving the equation with \033[1;33m(Newton Method)\033[0m" << std::endl;
        Newton_Method newton_solver(function, 0.5);
        double root_newton = newton_solver.solve();
        std::cout << "\033[33mRoot: " << root_newton << " ft\033[0m" << std::endl << std::endl;

        // 3. 割线法求解
        std::cout << "--------------------------------" << std::endl << std::endl;
        std::cout << "Solving the equation with \033[1;34m(Secant Method)\033[0m" << std::endl;
        Secant_Method secant_solver(function, 0.3, 0.7);
        double root_secant = secant_solver.solve();
        std::cout << "\033[34mRoot: " << root_secant << " ft\033[0m" << std::endl<< std::endl;

        std::cout << "--------------------------------" << std::endl << std::endl;


    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}