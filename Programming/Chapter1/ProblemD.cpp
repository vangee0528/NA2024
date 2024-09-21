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

void checkRoot(const Function &f, double root) {
    if(fabs(f(root)) < 1e-7) {
        std::cout << "\033[32mRoot is correct.\033[0m" << std::endl << std::endl;
        std::cout << "--------------------------------" << std::endl << std::endl;
    } else {
        std::cout << "\033[31mRoot is incorrect.\033[0m" << std::endl << std::endl;
        std::cout << "--------------------------------" << std::endl << std::endl;
    }
}


int main() {
    try {
        std::cout << "--------------------------------" << std::endl << std::endl;

        // 1. f(x) = sin(x/2) - 1
        {
            std::cout << "Testing f(x) = sin(x/2) - 1:" << std::endl;
            f1 function1;
            Secant_Method solver1(function1, 0, M_PI / 2);
            double root1 = solver1.solve();
            std::cout << "Root of sin(x/2) - 1: " << root1 << std::endl;
            std::cout << "f(root) = " << function1(root1) << std::endl;
            checkRoot(function1, root1);
        }

        // 2. f(x) = e^x - tan(x)
        {
            std::cout << "Testing f(x) = e^x - tan(x):" << std::endl;
            f2 function2;
            Secant_Method solver2(function2, 1, 1.4);
            double root2 = solver2.solve();
            std::cout << "Root of e^x - tan(x): " << root2 << std::endl;
            std::cout << "f(root) = " << function2(root2) << std::endl;
            checkRoot(function2, root2);
        }

        // 3. f(x) = x^3 - 12x^2 + 3x + 1
        {
            std::cout << "Testing f(x) = x^3 - 12x^2 + 3x + 1:" << std::endl;
            f3 function3;
            Secant_Method solver3(function3, 0, -0.5);
            double root3 = solver3.solve();
            std::cout << "Root of x^3 - 12x^2 + 3x + 1: " << root3 << std::endl;
            std::cout << "f(root) = " << function3(root3) << std::endl;
            checkRoot(function3, root3);
        }


    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
