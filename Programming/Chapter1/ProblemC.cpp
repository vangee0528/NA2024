#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Function.hpp"          // 引入 Function 虚类
#include "EquationSolver.hpp"    // 引入 EquationSolver 虚类

class Function2 : public Function {
public:
    virtual double operator() (double x) const override {
        return x - tan(x); // f(x) = x - tan(x)
    }

    virtual double derivative(double x) const override {
        return 1 - 1 / cos(x) / cos(x); // f'(x) = 1 - sec^2(x)
    }
};


int main() {
    try {
        Function2 f2;

        // 测试接近 4.5 的根
        double x0_1 = 4.5;
        Newton_Method solver1(f2, x0_1);
        double root1 = solver1.solve();
        std::cout << "Root near 4.5: " << root1 << std::endl;
        std::cout << "f(root) = " << f2(root1) << std::endl;

        // 测试接近 7.7 的根
        double x0_2 = 7.7;
        Newton_Method solver2(f2, x0_2);
        double root2 = solver2.solve();
        std::cout << "Root near 7.7: " << root2 << std::endl;
        std::cout << "f(root) = " << f2(root2) << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
