#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Function.hpp"          // 引入 Function 虚类
#include "EquationSolver.hpp"    // 引入 EquationSolver 虚类

// 定义函数 f1
class f1 : public Function {
public:
    virtual double operator() (double x) const override {
        return sin(x / 2) - 1; // f(x) = sin(x/2) - 1
    }
};

// 定义函数 f2
class f2 : public Function {
public:
    virtual double operator() (double x) const override {
        return exp(x) - tan(x); // f(x) = e^x - tan(x)
    }
};

// 定义函数 f3
class f3 : public Function {
public:
    virtual double operator() (double x) const override {
        return x * x * x - 12 * x * x + 3 * x + 1; // f(x) = x^3 - 12x^2 + 3x + 1
    }
};


int main() {
    try {
        // 1. f(x) = sin(x/2) - 1
        f1 function1;
        Secant_Method solver1(function1, 0, M_PI / 2);
        double root1 = solver1.solve();
        std::cout << "Root of sin(x/2) - 1: " << root1 << std::endl;
        std::cout << "f(root) = " << function1(root1) << std::endl;

        // 2. f(x) = e^x - tan(x)
        f2 function2;
        Secant_Method solver2(function2, 1, 1.4);
        double root2 = solver2.solve();
        std::cout << "Root of e^x - tan(x): " << root2 << std::endl;
        std::cout << "f(root) = " << function2(root2) << std::endl;

        // 3. f(x) = x^3 - 12x^2 + 3x + 1
        f3 function3;
        Secant_Method solver3(function3, 0, -0.5);
        double root3 = solver3.solve();
        std::cout << "Root of x^3 - 12x^2 + 3x + 1: " << root3 << std::endl;
        std::cout << "f(root) = " << function3(root3) << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
