#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>


// 定义函数 f(x) = 1/x - tan(x)
class Function1 : public Function {
public:
    virtual double operator() (double x) const override {
        if (x == 0) return INFINITY;  // 防止 1/x 在 x=0 处出现除零
        return 1.0 / x - std::tan(x);
    }
};


// 函数 f(x) = 1/x - 2^x
class Function2 : public Function {
public:
    virtual double operator() (double x) const override {
        if (x == 0) return INFINITY;  // 防止除零
        return 1.0 / x - std::pow(2, x);
    }
};

// 函数 f(x) = 2^{-x} + e^x + 2cos(x) - 6
class Function3 : public Function {
public:
    virtual double operator() (double x) const override {
        return std::pow(2, -x) + std::exp(x) + 2 * std::cos(x) - 6;
    }
};

// 函数 f(x) = (x^3 - 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2)
class Function4 : public Function {
public:
    virtual double operator() (double x) const override {
        double numerator = std::pow(x, 3) + 4 * std::pow(x, 2) + 3 * x + 5;
        double denominator = 2 * std::pow(x, 3) - 9 * std::pow(x, 2) + 18 * x - 2;
        if (denominator == 0) return INFINITY;  // 防止分母为零
        return numerator / denominator;
    }
};

int main() {

        // 测试函数1: f(x) = 1/x - tan(x) 在区间 [0, pi/2]
    {
        std::cout << "Testing f(x) = 1/x - tan(x) in [0, pi/2]:" << std::endl;
        Function1 f1;
        double a = 1e-6;
        double b = M_PI / 2 - 1e-6;  // 避免 tan(pi/2) 不定义
        Bisection_Method solver1(f1, a, b);
        double root1 = solver1.solve();
        std::cout << "Root: " << root1 << std::endl;
        std::cout << "f(root) = " << f1(root1) << std::endl << std::endl;
    }

    
    // 测试函数2: f(x) = 1/x - 2^x 在区间 [0, 1]
    {
        std::cout << "Testing f(x) = 1/x - 2^x in [0, 1]:" << std::endl;
        Function2 f2;
        double a = 1e-6;
        double b = 1;
        Bisection_Method solver2(f2, a, b);
        double root2 = solver2.solve();
        std::cout << "Root: " << root2 << std::endl;
        std::cout << "f(root) = " << f2(root2) << std::endl << std::endl;
    }

    // 测试函数3: f(x) = 2^{-x} + e^x + 2cos(x) - 6 在区间 [1, 3]
    {
        std::cout << "Testing f(x) = 2^{-x} + e^x + 2cos(x) - 6 in [1, 3]:" << std::endl;
        Function3 f3;
        double a = 1;
        double b = 3;
        Bisection_Method solver3(f3, a, b);
        double root3 = solver3.solve();
        std::cout << "Root: " << root3 << std::endl;
        std::cout << "f(root) = " << f3(root3) << std::endl << std::endl;
    }

    // 测试函数4: f(x) = (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2) 在区间 [0, 4]
    {
        std::cout << "Testing f(x) = (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2) in [0, 4]:" << std::endl;
        Function4 f4;
        double a = 0;
        double b = 4;
        Bisection_Method solver4(f4, a, b);
        std::cout << f4(0)<<f4(4) << std::endl;
        double root4 = solver4.solve();
        std::cout << "Root: " << root4 << std::endl;
        std::cout << "f(root) = " << f4(root4) << std::endl << std::endl;
    }

    return 0;
}
