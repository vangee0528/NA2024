#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Function.hpp"          // 引入 Function 虚类
#include "EquationSolver.hpp"    // 引入 EquationSolver 虚类

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

        // 1. 牛顿法求解
        Newton_Method newton_solver(function, 0.5);
        double root_newton = newton_solver.solve();
        std::cout << "Root (Newton's Method): " << root_newton << " ft" << std::endl;

        // 2. 割线法求解
        Secant_Method secant_solver(function, 0.3, 0.7);
        double root_secant = secant_solver.solve();
        std::cout << "Root (Secant Method): " << root_secant << " ft" << std::endl;

        // 3. 二分法求解
        Bisection_Method bisection_solver(function, 0, 1);
        double root_bisection = bisection_solver.solve();
        std::cout << "Root (Bisection Method): " << root_bisection << " ft" << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}